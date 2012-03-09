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

class Plane(object):
    """2D plane."""
    xrange = None
    yrange = None
    
    def __init__(self, xmax, ymax, xmin=0, ymin=0):
        self.xrange = (xmin, xmax)
        self.yrange = (ymin, ymax)

    length = property(lambda self: self.xrange[1] - self.xrange[0])
    height = property(lambda self: self.yrange[1] - self.yrange[0])

    def wrap(self, p, bounds):
        # shift to make the math easier
        shift = 0
        if bounds[0] < 0:
            shift = -bounds[0]
            bounds = (0, bounds[1] + shift)
            p += shift
        total = bounds[1] - bounds[0]
        if p < bounds[0]:
            diff = (bounds[0] - p) % total
            p = bounds[1] - diff
        elif p > bounds[1]:
            diff = (p - bounds[1]) % total
            p = bounds[0] + diff
        if shift != 0:
            p -= shift
        return p
    
    def wrapx(self, x):
        return self.wrap(x, self.xrange)
            
    def wrapy(self, y):
        return self.wrap(y, self.yrange)

##############################################################################
##############################################################################

class LonLatPlane(Plane):
    """Longitude-latitude plane."""
    
    LAT_MIN = -90.0
    LAT_MAX = 90.0
    LON_MIN = -180.0
    LON_MAX = 180.0

    def __init__(self):
        super(LonLatPlane, self).__init__(self.LON_MAX, self.LAT_MAX, 
                                          self.LON_MIN, self.LAT_MIN)
    
##############################################################################
##############################################################################
    
class Grid(object):
    """Overlays a grid of fixed-size cells on a 2D plane."""
    nrows = None
    ncols = None
    
    def __init__(self, plane, nrows, ncols=None):
        self.plane = plane
        self.nrows = nrows
        if ncols is None:
            ncols = nrows
        self.ncols = ncols
    
    ncells = property(lambda self: self.nrows * self.ncols)
    cell_height = property(lambda self: self.plane.height / self.nrows)
    cell_width = property(lambda self: self.plane.length / self.ncols)
    
    def cells(self):
        """Returns a list of the bottom-left points in all cells."""
        cells = []
        xinc = self.cell_width
        yinc = self.cell_height
        x = self.plane.xrange[0]
        while x < self.plane.xrange[1]:
            y = self.plane.yrange[0]
            while y < self.plane.yrange[1]:
                cells.append(Point(x, y))
                y += yinc
            x += xinc
        return cells
    
    def bin(self, p):
        """Assign point p to a cell."""
        xinc = self.cell_width
        yinc = self.cell_height
        # shift to make the math easier
        px, py = p
        xmin = self.plane.xrange[0]
        ymin = self.plane.yrange[0]
        if xmin < 0:
            px -= xmin
        if ymin < 0:
            py -= ymin 
        x = math.floor(px / xinc) * xinc
        y = math.floor(py / yinc) * yinc
        if xmin < 0:
            x += xmin
        if ymin < 0:
            y += ymin 
        return Point(x, y)

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
        bin = grid.bin(p)
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
    xinc = grid.cell_width
    yinc = grid.cell_height
    seen = Set()
    cells = []
    # Iteratively explore the neighboring cells an increasing distance away
    # from the given cell
    while len(seen) < grid.ncells:
        if distance == 0:
            cells.append(cell)
        else:
            xmin = grid.plane.wrapx(cell.x - xinc * distance)
            ymin = grid.plane.wrapy(cell.y - yinc * distance)
            xs = []
            ys = []
            for min, inc, wrap, bounds, result in ((xmin, xinc, grid.plane.wrapx, grid.plane.xrange, xs),
                                                   (ymin, yinc, grid.plane.wrapy, grid.plane.yrange, ys)):
                coord = min
                for i in xrange(distance*2+1):
                    result.append(coord)
                    coord += inc
                    coord = wrap(coord)
                    if coord == bounds[1]:
                        coord = bounds[0]
                    if coord == min and i > 0:
                        break
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
    cells = grid.cells()
    if len(cells) < num_nodes:
        raise RuntimeError("Grid too coarse")
    if len(cells) > num_nodes:
        cells = random.sample(cells, num_nodes)
    # maybe else: randomly permute the cells

    # for each grid cell, pick a node from a close site
    selection = Set()
    for cell in cells:
        site_id = select_site(cell, grid, bins)
        site = sites_by_id[site_id]
        # pick random available node from site
        nodes = [n for n in site['node_ids'] if (n in nodes_by_id and n not in selection)]
        assert len(nodes) > 0
        node = random.choice(nodes)
        selection.add(node)
        if len(nodes) == 1:
            # no more available nodes from site, so remove it
            bin = grid.bin(site_to_point(site))
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
        
    opts, args = parser.parse_args()
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
        argv = sys.argv
    opts = parse_options(argv)
    
    shell = connect(opts.user, opts.password, opts.url)
    slice, all_sites, all_nodes = fetch(shell, opts.slice, ('run_level', 'boot_state',))

    # filter nodes and index by id
    node_filters = [Filter.filter_booted, Filter.filter_member(slice['node_ids'])]
    nodes_by_id = {}
    for n in all_nodes:
        if Filter.apply(node_filters, n):
            nodes_by_id[n['node_id']] = n

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
    sys.exit(main(sys.argv))
    
##############################################################################
##############################################################################
