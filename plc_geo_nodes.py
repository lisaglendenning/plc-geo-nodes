#!/bin/env python

"""Allocates geographically diverse nodes to a PL slice using PLCAPI.

Uses a heuristic of dividing the latitude/longitude plane into a grid of cells.
Then for each cell in the grid, picks a close site, then picks a random node
from that site.
"""

##############################################################################
##############################################################################

import sys, optparse, getpass, math, random
from sets import Set # assuming python can be less than 2.6

import PLC.Shell

##############################################################################
##############################################################################

PLCAPI_URL = 'https://www.planet-lab.org/PLCAPI/'

USAGE = "%prog [options] USER"
DESCRIPTION = """Allocates geographically dispersed nodes to a using PLCAPI.
Required argument USER: PLC user or file containing PLC user.
Prompts for the password if not given on the command line.
Will choose some slice if no slice is specified.
"""

DEFAULT_NUMBER = 16 # default number of nodes to allocate

##############################################################################
##############################################################################

class Filter(object):
    """Useful filters for node and site objects."""
    
    @staticmethod
    def apply(filters, obj):
        for f in filters:
            if not f(obj):
                return False
        return True

    @staticmethod
    def filter_booted(node):
        # filter nodes that are likely offline
        boot_state = node['boot_state']
        if boot_state not in ('boot',):
            return False
        run_level = node['run_level']
        if run_level not in ('boot',):
            return False
        return True
    
    @staticmethod
    def filter_member(nodes):
        # filter nodes that are part of some set
        def inner(node):
            if node in nodes:
                return False
            return True
        return inner

    @staticmethod
    def filter_locatable(site):
        # filter sites without a valid long/lat
        p = (site['longitude'], site['latitude'])
        return None not in p
        
    @staticmethod
    def filter_available(nodes):
        # filter sites with no available nodes
        def inner(site):
            available = [n for n in site['node_ids'] if n in nodes]
            return len(available) > 0
        return inner

##############################################################################
##############################################################################

def connect(user, password, url=PLCAPI_URL):
    """Initializes PLC connection."""
    shell = PLC.Shell.Shell(globals = globals(),
                        url = url, 
                        xmlrpc = True,
                        method = 'password',
                        user = user, 
                        password = password)
    return shell

##############################################################################
##############################################################################

def fetch(shell, slice_name, node_fields=(), site_fields=()):
    """Returns a slice object, all nodes, and all sites."""
    return_fields = ['name', 'slice_id', 'node_ids']
    # keyword arguments don't seem to work ?
    slices = GetSlices(slice_name, return_fields)
    if len(slices) < 1:
        raise RuntimeError("No slices found")
    # if slice is not specified, use the first slice returned
    slice = slices[0]
    
    return_fields = ['site_id', 'abbreviated_name', 'name', 
                     'latitude', 'longitude', 'node_ids',]
    for f in site_fields:
        if f not in return_fields:
            return_fields.append(f)
    sites = GetSites(None, return_fields)
    return_fields = ['node_id', 'hostname', 'site_id',]
    for f in node_fields:
        if f not in return_fields:
            return_fields.append(f)
    nodes = GetNodes(None, return_fields)
    return slice, sites, nodes
    
##############################################################################
##############################################################################

def update(shell, slice, nodes):
    """Adds nodes to slice."""
    AddSliceToNodes(slice['slice_id'], list(nodes))
    
##############################################################################
##############################################################################

# namedtuple would be perfect
class Point(tuple):
    """2D point."""

    def __new__(cls, x, y):
        return super(Point, cls).__new__(cls, (x, y))
    
    x = property(lambda self: self[0])
    y = property(lambda self: self[1])

##############################################################################
##############################################################################

class Axis(tuple):
    
    def __new__(cls, max, min=0):
        return super(Axis, cls).__new__(cls, (min, max))
    
    min = property(lambda self: self[0])
    max = property(lambda self: self[-1])
    length = property(lambda self: self.max - self.min)

    def wrap(self, p):
        min, max = self
        # shift to make the math easier
        shift = 0
        if min < 0:
            shift = -min
            min, max = (0, max + shift)
            p += shift
        total = max - min
        if p < min:
            diff = (min - p) % total
            p = max - diff
        elif p > max:
            diff = (p - max) % total
            p = min + diff
        if shift != 0:
            p -= shift
        return p
    
##############################################################################
##############################################################################

class Plane(object):
    """2D plane."""
    x = None
    y = None
    
    def __init__(self, xmax, ymax, xmin=0, ymin=0):
        self.x = Axis(xmax, xmin)
        self.y = Axis(ymax, ymin)

    length = property(lambda self: self.x.size)
    height = property(lambda self: self.y.size)

##############################################################################
##############################################################################

class LonLatPlane(Plane):
    """Longitude-latitude plane."""
    
    LAT_MIN = -90
    LAT_MAX = 90
    LON_MIN = -180
    LON_MAX = 180

    def __init__(self):
        super(LonLatPlane, self).__init__(self.LON_MAX, self.LAT_MAX, 
                                          self.LON_MIN, self.LAT_MIN)
    
##############################################################################
##############################################################################

class GridAxis(object):
    
    axis = None
    tics = None
    
    @classmethod
    def create_tics(cls, axis, ntics):
        tics = []
        min, max = [int(i) for i in axis]
        total = max - min
        inc = total / ntics
        tics = [min + i*inc for i in xrange(0, ntics)]
        return tics
    
    def __init__(self, axis, ntics):
        self.axis = axis
        self.tics = self.create_tics(axis, ntics)
        assert len(self.tics) == ntics
    
    ntics = property(lambda self: len(self.tics))
    inc = property(lambda self: int(self.axis.length / self.ntics))
    
    def bin(self, i):
        # shift to zero-base to make the math easier
        min = self.axis.min
        if min != 0:
            i -= min
        bin = int(i / self.inc)
        ntics = self.ntics
        if bin > ntics - 1:
            bin = ntics - 1
        assert bin >= 0 and bin < ntics
        return bin
    
    def map(self, i):
        return self.tics[self.bin(i)]
    
    def next(self, i):
        current = self.bin(i)
        if current == self.ntics - 1:
            next = 0
        else:
            next = current + 1
        return self.tics[next]
    
    def prev(self, i):
        current = self.bin(i)
        if current == 0:
            prev = self.ntics - 1
        else:
            prev = current - 1
        return self.tics[prev]
    
##############################################################################
##############################################################################

class Grid(object):
    """Overlays a grid of cells on a 2D plane. Cell dimensions are integral."""
    
    DIRECTIONS = xrange(4)
    SOUTH, NORTH, WEST, EAST = DIRECTIONS
        
    @classmethod
    def create_cells(cls, xtics, ytics):
        """Returns a list of the bottom-left points in all cells."""
        cells = []
        for x in xtics:
            for y in ytics:
                cells.append(Point(x, y))
        return cells
    
    def __init__(self, plane, nrows, ncols=None):
        self.plane = plane
        if ncols is None:
            ncols = nrows
        self.rows = GridAxis(plane.y, nrows)
        self.cols = GridAxis(plane.x, ncols)
        self.cells = Set(self.create_cells(self.cols.tics, self.rows.tics))
    
    nrows = property(lambda self: self.rows.ntics)
    ncols = property(lambda self: self.cols.ntics)
    ncells = property(lambda self: self.nrows * self.ncols)
    col_inc = property(lambda self: self.cols.inc)
    row_inc = property(lambda self: self.rows.inc)
    
    def map(self, p):
        """Return cell containing point p."""
        px, py = p
        x = self.cols.map(px)
        y = self.rows.map(py)
        cell = Point(x, y)
        assert cell in self.cells
        return cell
    
    def neighbor(self, p, direction):
        x, y = self.cols.map(p.x), self.rows.map(p.y)
        if direction in (self.NORTH, self.SOUTH):
            if (direction == self.NORTH):
                y = self.rows.next(y)
            else:
                y = self.rows.prev(y)
        else:
            if (direction == self.EAST):
                x = self.cols.next(x)
            else:
                x = self.cols.prev(x)
        neighbor = Point(x, y)
        return neighbor

##############################################################################
##############################################################################

def site_to_point(site):
    return Point(site['longitude'], site['latitude'])

##############################################################################
##############################################################################

def bin_sites(sites_by_id, grid):
    """Assign each site to a cell of the grid."""
    bins = {}
    for site_id, site in sites_by_id.iteritems():
        p = site_to_point(site)
        bin = grid.map(p)
        if bin not in bins:
            bins[bin] = Set()
        bins[bin].add(site_id)
    return bins

##############################################################################
##############################################################################

def select_site(cell, grid, bins):
    """Select a close site for the given cell."""
    site = None
    distance = 0
    seen = Set()
    cells = []
    # Iteratively explore the neighboring cells an increasing distance away
    # from the given cell
    while len(seen) < grid.ncells:
        if distance == 0:
            cells.append(cell)
        else:
            xs = []
            ys = []
            for coord, axis, coords in ((cell.x, grid.cols, xs), (cell.y, grid.rows, ys)):
                c = coord
                coords.append(c)
                for i in xrange(1, distance+1):
                    c = axis.prev(c)
                    if c == coords[-1]:
                        break
                    coords.append(c)
                coords.reverse()
                c = coord
                for i in xrange(1, distance+1):
                    c = axis.next(c)
                    if c == coords[0]:
                        break
                    coords.append(c)
            if len(xs) == 0 or len(ys) == 0:
                break
            for xvals, yvals in ((xs, (ys[0], ys[-1])), ((xs[0], xs[-1]), ys)):
                for x in xvals:
                    for y in yvals:
                        point = Point(x, y)
                        if point not in cells and point not in seen:
                            cells.append(point)
            random.shuffle(cells)
        for next in cells:
            seen.add(next)
            if next in bins:
                site = random.choice(list(bins[next]))
                return site
        distance += 1
        del cells[:]
    return site
            
##############################################################################
##############################################################################

def allocate(sites_by_id, nodes_by_id, num_nodes=DEFAULT_NUMBER, 
             nrows=None, ncols=None):
    """Returns a set of num_node nodes from nodes_by_id."""
    # initialize grid
    if nrows is None:
        if ncols is not None:
            nrows = ncols
        else:
            # grid has the same number of rows and columns
            # round up so that the number of grid cells >= num_nodes
            nrows = int(math.ceil(math.sqrt(num_nodes)))
    plane = LonLatPlane()
    grid = Grid(plane, nrows, ncols)
    bins = bin_sites(sites_by_id, grid)
    cells = grid.cells
    if len(cells) < num_nodes:
        raise RuntimeError("Grid too coarse")
    elif len(cells) > num_nodes:
        cells = random.sample(cells, num_nodes)
    # maybe else: randomly permute the cells

    # for each grid cell, pick a node from a close site
    selection = Set()
    for cell in cells:
        site_id = select_site(cell, grid, bins)
        assert site_id is not None
        site = sites_by_id[site_id]
        # pick random available node from site
        nodes = [n for n in site['node_ids'] if (n in nodes_by_id and n not in selection)]
        assert len(nodes) > 0
        node = random.choice(nodes)
        selection.add(node)
        if len(nodes) == 1:
            # no more available nodes from site, so remove it
            bin = grid.map(site_to_point(site))
            bins[bin].remove(site_id)
            if len(bins[bin]) == 0:
                del bins[bin]
    return selection

##############################################################################
##############################################################################

def parse_options(argv):
    parser = optparse.OptionParser(usage=USAGE, description=DESCRIPTION)
    parser.add_option("-u", "--url", help="PLCAPI URL (default: %s)" % PLCAPI_URL, 
                      default=PLCAPI_URL)
    parser.add_option("-p", "--password", 
                      help="PLC password or file containing PLC password (default: prompt)")
    parser.add_option("-s", "--slice", 
                      help="Slice name (default: some slice)")
    parser.add_option("-n", "--nodes",
                      type=int, default=DEFAULT_NUMBER,
                      help="Number of nodes to add (default=%d)" % DEFAULT_NUMBER)
    parser.add_option("-r", "--rows",
                      type=int,
                      help="Number of grid rows")
    parser.add_option("-c", "--cols",
                      type=int,
                      help="Number of grid columns")
        
    opts, args = parser.parse_args(argv)
    if len(args) < 1:
        parser.error("Missing required argument (USER)")
    opts.user = args[0]
    try:
        opts.user = open(opts.user).read().strip()
    except IOError:
        pass
    if opts.password is None:
        try:
            opts.password = getpass.getpass()
        except (EOFError, KeyboardInterrupt):
            return 0
    else:
        try:
            opts.password = open(opts.password).read().strip()
        except IOError:
            pass
    return opts

##############################################################################
##############################################################################

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    opts = parse_options(argv)
    
    shell = connect(opts.user, opts.password, opts.url)
    slice, all_sites, all_nodes = fetch(shell, opts.slice, ('run_level', 'boot_state',))

    # filter nodes and index by id
    node_filters = [Filter.filter_booted, Filter.filter_member(slice['node_ids'])]
    nodes_by_id = {}
    for n in all_nodes:
        if Filter.apply(node_filters, n):
            nodes_by_id[n['node_id']] = n
    if len(nodes_by_id) < opts.nodes:
        raise RuntimeError("Only %d nodes available" % len(nodes_by_id))

    # filter sites and index by id
    site_filters = [Filter.filter_locatable, Filter.filter_available(nodes_by_id)]
    sites_by_id = {}
    for s in all_sites:
        if Filter.apply(site_filters, s):
            sites_by_id[s['site_id']] = s
    
    selection = allocate(sites_by_id, nodes_by_id, opts.nodes, opts.rows, opts.cols)
    update(shell, slice, selection)
    
##############################################################################
##############################################################################

if __name__ == '__main__':
    sys.exit(main())
    
##############################################################################
##############################################################################
