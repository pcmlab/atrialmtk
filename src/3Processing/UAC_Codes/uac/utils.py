#!/usr/bin/python3
#
import os
import re
import shutil
import struct
import numpy as np
import collections
from scipy.spatial import distance
import vtk
from vtk.util import numpy_support
from uac.element import Element
from uac.line import Line


# Number of vertices per element
eleminfo = {"Ln": 2, "Tr": 3, "Qd": 4, "Tt": 4, "Py": 5, "Pr": 6, "Hx": 8}
# Map number of vertices to CARP surface element
carpsurf = {2: "Ln", 3: "Tr", 4: "Qd"}
# Map CARP to VTK element types
vtkugrid_elemtypes = {"Ln": 3, "Tr": 5, "Qd": 9, "Tt": 10, "Py": 14, "Pr": 13, "Hx": 11}

verbose = 0


def find_index(list_to_find, elem_selection_func):
    elem = elem_selection_func(list_to_find)
    for index in range(len(list_to_find)):
        if list_to_find[index] == elem:
            break
    return index


def closest_node(node, nodes):
    deltas = nodes - node
    dist_2 = np.einsum("ij,ij->i", deltas, deltas)
    return dist_2


def load_landmarks(filepath, scale_factor):
    if filepath.endswith(".vtk"):
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(filepath)
        reader.Update()
        polydata = reader.GetOutput()
        size = polydata.GetNumberOfPoints()
        points = np.empty((size, 3), float)
        for i in range(size):
            points[i] = np.array(polydata.GetPoint(i))
    else:
        points = np.loadtxt(filepath, delimiter=",")

    return points * scale_factor


def read_vtk_surface(surface):
    surface_filter = vtk.vtkDataSetSurfaceFilter()
    surface_filter.SetInputData(surface)
    surface_filter.Update()
    polydata = surface_filter.GetOutput()

    size = polydata.GetNumberOfPoints()
    tot_pts = np.empty((size, 3), float)
    for i in range(size):
        tot_pts[i] = np.array(polydata.GetPoint(i))

    return tot_pts, surface_filter


def poly2carp(polydata):
    pts = numpy_support.vtk_to_numpy(polydata.GetPoints().GetData())
    ele = numpy_support.vtk_to_numpy(polydata.GetPolys().GetData())
    i = 0
    elems = []
    while i < len(ele):
        elems.append(ele[i + 1: i + ele[i] + 1])
        i += ele[i] + 1
    elems = [[carpsurf[len(e)]] + list(map(str, e)) for e in elems]
    return pts, elems


def poly2carp_with_labels(polydata):
    pts, elems = poly2carp(polydata)
    if polydata.GetCellData().GetNumberOfArrays():
        labels = numpy_support.vtk_to_numpy(polydata.GetCellData().GetScalars())
        for i in range(len(elems)):
            elems[i].append(str(labels[i]))
    return pts, elems


def convert_unstructureddata_to_polydata(surface):
    """Converts unstructured grid to polydata"""
    ugg = vtk.vtkGeometryFilter()
    if vtk.VTK_MAJOR_VERSION < 6:
        ugg.SetInput(surface)
    else:
        ugg.SetInputData(surface)
    ugg.Update()
    return ugg.GetOutput()


def clean_polydata(polydata):
    """removes floating nodes etc"""
    cleaner = vtk.vtkCleanPolyData()
    if vtk.VTK_MAJOR_VERSION < 6:
        cleaner.SetInput(polydata)
    else:
        cleaner.SetInputData(polydata)
    cleaner.Update()
    return cleaner.GetOutput()


def get_vtk_from_file(filename):
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllScalarsOn()
    reader.ReadAllVectorsOn()
    reader.ReadAllTensorsOn()
    reader.Update()
    return reader.GetOutput()


def find_edges(elems, pts):
    """
    Find border edges, and make a dictionary of elements across each edge also return the set of nodes on the border
    """
    if verbose >= 1:
        print("Finding border edges")

    bedges = dict()
    for e in elems:
        ees = e.edges()
        for edge in ees:
            if edge in bedges:
                del bedges[edge]
            else:
                bedges[edge] = e.index

    edgenode = set()
    for b in bedges:
        edgenode.add(b.n[0])
        edgenode.add(b.n[1])

    if verbose >= 2:
        bdat = [0] * len(pts)
        for n in edgenode:
            bdat[n] = 1
        print(bdat)

    return bedges, edgenode


# Join boundaries into connected lists
def coalesce(edges):
    """Join disparate edges into contiguous lines"""

    ##############################################
    def coalesce_loop(neighbours, line):
        """Function to spread to edge neighbours"""
        victims = set()

        for n in neighbours:
            edgelist = start.get(n)
            if edgelist:
                for e in edgelist:
                    line.add_edge(e)
                    victims.add(e.n[1])
                start.pop(n, None)
            edgelist = end.get(n)
            if edgelist:
                for e in edgelist:
                    line.add_edge(e)
                    victims.add(e.n[0])
                end.pop(n, None)
        return victims

    ##############################################

    if verbose >= 1:
        print(("Coalescing lines from", str(len(edges)), "edges"))

    start = collections.defaultdict(set)
    end = collections.defaultdict(set)
    for e in edges:
        start[e.n[0]].add(e)
        end[e.n[1]].add(e)

    lines = []
    while start:
        lines.append(Line())

        # start at an end if possible
        start_n = None
        for n in start:
            if n not in end:
                start_n = n
                break
        if not start_n:
            if n not in start:
                start_n = n
                break
        if not start_n:
            start_n = list(start.keys())[0]

        neighbours = set((start_n,))

        while neighbours:
            neighbours = coalesce_loop(neighbours, lines[-1])

        if verbose >= 2:
            print((len(lines[-1].edge), lines[-1].nodes))

    return lines


def read_array_fiber(lon_file):
    """Read a .lon file from CARP simulator with no indicator in the beginning"""
    with open(lon_file) as file:
        first_line = file.readline()

    skip_count = 1 if len(first_line.split()) == 1 else 0

    return np.loadtxt(lon_file, skiprows=skip_count)


def read_array_igb(igbfile):
    """Read a .igb file"""
    with open(igbfile, mode="rb") as file:
        header = file.read(1024)
        words = header.split()
        word = []
        for i in range(4):
            word.append(int([re.split(r"(\d+)", s.decode("utf-8")) for s in [words[i]]][0][1]))

        nnode = word[0] * word[1] * word[2]

        size = os.path.getsize(igbfile) // 4 // nnode
        data = np.empty((size, nnode), float)
        for i in range(size):
            data[i] = struct.unpack("f" * nnode, file.read(4 * nnode))

    return data


def write_array_igb(igbfile, data):
    """Write data to a .igb file"""
    with open(igbfile, mode="wb") as file:
        npts = len(data) if isinstance(data[0], float) else len(data[0])
        ntime = len(data)
        header = "x:" + str(npts) + " y:1 z:1 t:" + str(ntime) + " type:float"
        for i in range(1022 - len(header)):
            header = header + " "
        header = header + "\n\f"
        header = "".join("\n" if i % 80 == 0 else char for i, char in enumerate(header, 1))
        file.write(header.encode("utf8"))

        if ntime != npts:
            for n in range(ntime):
                for i in range(npts):
                    file.write(struct.pack("f", data[n][i]))
        else:
            for i in range(npts):
                file.write(struct.pack("f", data[i]))


def read_pts(modname, n=3, vtx=False, item_type=float):
    """Read pts/vertex file"""
    with open(modname + (".vtx" if vtx else ".pts")) as file:
        count = int(file.readline().split()[0])
        if vtx:
            file.readline()

        pts = np.empty((count, n), item_type)
        for i in range(count):
            pts[i] = [item_type(val) for val in file.readline().split()[0:n]]

    return pts if n > 1 else pts.flat


def read_element_file(filename):
    """Read element file"""
    with open(filename) as file:
        size = int(file.readline().split()[0])
        elems = [None] * size
        for i in range(size):
            e = file.readline().split()
            np = eleminfo[e[0]]
            elems[i] = Element(e[0], e[1: np + 1], e[np + 1] if len(e) > np + 1 else 0)

    return elems


def read_elem(modname):
    """Read element file"""
    return read_element_file(modname + ".elem")


def read_surf(modname):
    """Read surf file"""
    return read_element_file(modname + ".surf")


def read_carp(base_dir, mesh_name, return_surface=True):
    carp_mesh = base_dir + mesh_name
    lon_file = carp_mesh + ".lon"
    igb_file = carp_mesh + ".igb"
    pts = read_pts(carp_mesh)
    elems = read_elem(carp_mesh)
    fiber = read_array_fiber(lon_file) if os.path.isfile(lon_file) else None
    data = read_array_igb(igb_file) if os.path.isfile(igb_file) else None
    if not return_surface:
        return pts, elems, fiber, data

    tmp_file = base_dir + "tmp.vtk"
    write_vtk(pts, elems, fiber, data, tmp_file)
    surface = get_vtk_from_file(tmp_file)
    os.remove(tmp_file)
    surface = convert_unstructureddata_to_polydata(surface)

    return pts, elems, fiber, data, surface


def to_ele(elems, data):
    """Transfer nodal data to elements"""
    elems_len = len(elems)
    col_count = data.shape[1] if len(data.shape) > 1 else 1
    ls = np.empty((elems_len, col_count), float)
    for i in range(elems_len):
        for j in range(col_count):
            ls[i][j] = elems[i].node2ctr(data[:, j] if col_count > 1 else data)

    return ls


def write_dat_simple(data, filename):
    with open(filename, "w") as file:
        for d in data:
            file.write("".join(str(d)) + "\n")


def add_celldata_array(polydata, arrayname, initval=0):
    """Add a new point data array to polydata with name arrayname"""
    newarray = vtk.vtkDoubleArray()
    newarray.SetName(arrayname)
    newarray.SetNumberOfTuples(polydata.GetNumberOfCells())
    newarray.FillComponent(0, initval)
    polydata.GetCellData().AddArray(newarray)
    return polydata


def add_labels(surface, reg_name, src_regs):
    """
    Adds a new array to a surface and marks with marker (default 1.0).
    All the point within a circular region around a centre and with a specified radius.
    Unmarked points are marked with unmarked (default 0.0)
    """
    surface = add_celldata_array(surface, reg_name)
    regs = surface.GetCellData().GetArray(reg_name)
    for point_id in range(surface.GetNumberOfCells()):
        regs.SetValue(point_id, src_regs[point_id])
    return surface


def write_surface_to_vtk(surface, filename, write_ascii=False):
    writer = vtk.vtkPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION < 6:
        writer.SetInput(surface)
    else:
        writer.SetInputData(surface)
    if write_ascii:
        writer.SetFileTypeToASCII()
    else:
        writer.SetFileTypeToBinary()
    writer.SetFileName(filename)
    writer.Write()


def write_surface_to_carp_mesh(surface, outfile, labels=[]):
    """Writes the surface in the carp mesh format (.pts and .elem)"""
    write_carp_nodes(surface, outfile)

    def generate_line(elem, label=0):
        e = surface.GetCell(elem)
        line = "Tr {} {} {} {}\n".format(e.GetPointId(0), e.GetPointId(1), e.GetPointId(2), label)
        return line.encode()

    with open(outfile + ".elem", "wb") as file:
        line = str(surface.GetNumberOfCells()) + "\n"
        file.write(line.encode())
        if len(labels) == 0:
            # no region label
            for elem in range(surface.GetNumberOfCells()):
                file.write(generate_line(elem))
        else:
            for elem, elem_label in zip(list(range(surface.GetNumberOfCells())), labels):
                file.write(generate_line(elem, elem_label))


def write_carp_nodes(surface, outfile):
    """Writes the .pts file for carp meshes"""
    with open(outfile + ".pts", "wb") as file:
        line = str(surface.GetNumberOfPoints()) + "\n"
        file.write(line.encode())
        for i in range(surface.GetNumberOfPoints()):
            point = surface.GetPoint(i)
            line = "{} {} {}\n".format(point[0], point[1], point[2])
            file.write(line.encode())


def write_carp_nodes_uacs_3d(coords, outfile):
    """Writes the .pts file for carp meshes"""
    with open(outfile + ".pts", "wb") as file:
        line = str(len(coords)) + "\n"
        file.write(line.encode())
        for point in coords:
            line = "{} {} {}\n".format(point[0], point[1], point[2])
            file.write(line.encode())


def write_vtx_strength(ids, strength, filename, intra):
    with open(filename + ".vtx", "wb") as fp:
        line = str(len(ids)) + "\n"
        fp.write(line.encode())
        fp.write("intra\n".encode() if intra else "extra\n".encode())
        for i in range(len(ids)):
            line = "{} {}\n".format(ids[i], strength[i])
            fp.write(line.encode())


def write_vtx_strength_intra(ids, strength, filename):
    write_vtx_strength(ids, strength, filename, True)


def write_vtx_strength_extra(ids, strength, filename):
    write_vtx_strength(ids, strength, filename, False)


def write_vtk(pts, elems, fibers, data, outmesh):
    with open(outmesh, "w") as outfile:
        outfile.write("# vtk DataFile Version 2.0\n")
        outfile.write("converted from a CARP mesh\n")
        outfile.write("ASCII\n")
        outfile.write("DATASET UNSTRUCTURED_GRID\n")
        outfile.write("POINTS %d float\n" % len(pts))
        for p in pts:
            outfile.write("{0[0]} {0[1]} {0[2]}\n".format(p))
        outfile.write("\n")

        numcells = sum(eleminfo[e.type] + 1 for e in elems)
        outfile.write("CELLS {0} {1}\n".format(len(elems), numcells))
        for e in elems:
            outfile.write(str(eleminfo[e.type]) + " " + " ".join(map(str, e.n)) + "\n")
        outfile.write("\n")

        outfile.write("CELL_TYPES {0}\n".format(len(elems)))
        for e in elems:
            outfile.write(str(vtkugrid_elemtypes[e.type]) + "\n")
        outfile.write("\n")

    if data is not None and len(data[0]) != len(elems):
        ntime = len(data)
        if ntime == 1:
            with open(outmesh, "a") as outfile:
                outfile.write("POINT_DATA {0}\n".format(len(pts)))
                outfile.write("SCALARS data double\n")
                outfile.write("LOOKUP_TABLE default\n")
                for d in data:
                    outfile.write("{0}\n".format(d))
        else:
            if os.path.exists(outmesh[:-4]):
                shutil.rmtree(outmesh[:-4])
            os.makedirs(outmesh[:-4])
            for n in range(ntime):
                shutil.copyfile(outmesh, outmesh[:-4] + "/time{0}.vtk".format(str(n).zfill(4)))
                with open(outmesh[:-4] + "/time{0}.vtk".format(str(n).zfill(4)), "a") as f:
                    f.write("POINT_DATA {0}\n".format(len(pts)))
                    f.write("SCALARS data double\n")
                    f.write("LOOKUP_TABLE default\n")
                    for d in data[n]:
                        f.write("{0}\n".format(d))

    if fibers is not None:
        with open(outmesh, "a") as outfile:
            outfile.write("CELL_DATA {0}\n".format(len(elems)))
            outfile.write("VECTORS fibers float\n")
            for f in fibers:
                outfile.write("{0[0]} {0[1]} {0[2]}\n".format(f))
            outfile.write("\n")

    if data is None and elems[0].reg > 0:
        with open(outmesh, "a") as outfile:
            if fibers is None:
                outfile.write("CELL_DATA {0}\n".format(len(elems)))
            outfile.write("SCALARS data double\n")
            outfile.write("LOOKUP_TABLE default\n")
            for e in elems:
                outfile.write("{0}\n".format(e.reg))
            outfile.write("\n")


def write_carp(pts, elems, fiber, data, filename):
    with open(filename + ".pts", "w") as ptout:
        ptout.write(str(len(pts)) + "\n")
        for p in pts:
            ptout.write("{0[0]} {0[1]} {0[2]}\n".format(p))

    with open(filename + ".elem", "w") as eout:
        eout.write(str(len(elems)) + "\n")
        if isinstance(elems[0], Element):
            for e in elems:
                eout.write(e.type + " " + " ".join(map(str, e.n)) + " " + str(e.reg) + "\n")
        else:
            for e in elems:
                eout.write(" ".join(e) + "\n")

    with open(filename + ".lon", "w") as lonout:
        if fiber is None:
            for i in range(len(elems)):
                lonout.write("1 0 0\n")
        else:
            for f in fiber:
                lonout.write("{0[0]} {0[1]} {0[2]}\n".format(f))

    if data is not None:
        write_array_igb(filename + ".igb", data)


def find_junc_nodes(elems, reg1, reg2):
    """Find nodes at junction of the given regions"""
    nodes = [set(), set()]
    for e in elems:
        # Find list of nodes in the regions
        if e.reg == reg1:
            nodes[0].update(e.n)
        elif e.reg == reg2:
            nodes[1].update(e.n)

    tot_nodes = nodes[0].intersection(nodes[1])
    return np.fromiter(tot_nodes, int, len(tot_nodes))


def find_region_nodes(elems, reg):
    """Find nodes in the given region"""
    nodes = set()
    for e in elems:
        # Find list of nodes in the region
        if e.reg == reg:
            nodes.update(e.n)
    return np.fromiter(nodes, int, len(nodes))


def find_nodes(elems, pts, no=-1, ascending=False):
    """Finds MV nodes as those in the largest connected list of edges"""
    border, _ = find_edges(elems, pts)
    lines = coalesce(border)

    size = len(lines)
    line_size = np.empty(size, int)
    for i in range(size):
        line_size[i] = len(lines[i].nodes)

    if no == -1:
        # the largest line corresponds to the MV
        index = np.argmax(line_size)
    else:
        # sort and select which one
        ls = np.argsort(line_size)
        if not ascending:
            ls = ls[::-1]
        index = ls[no]

    nodes = np.fromiter(lines[index].nodes, int, len(lines[index].nodes))
    return index, lines, nodes


def boundary_markers(index1, index2, tot_pts):
    """Calculates markers 1-4 in UAC paper"""
    indices = [index1, index2]
    points_t = [None] * 2
    for k in range(2):
        size = len(indices[k])
        points = np.empty((size, 3), float)
        for i in range(size):
            for j in range(3):
                points[i][j] = tot_pts[indices[k][i]][j]

        points_t[k] = points

    y = distance.cdist(points_t[0], points_t[1], "euclidean")
    indices = np.unravel_index(np.argmin(y, axis=None), y.shape)

    return indices[0], indices[1]


def element2delete(elems, node_del, elems_del_src=None, excluded=False):
    """
    Returns elements (not) to delete
        elems_del_src: defines the deletion source; it uses this list to create the list instead of elems
        excluded: if it's True, returns elements not to delete
    """
    if elems_del_src is None:
        elems_del_src = elems

    node_del = set(node_del.flat)
    elems_len = len(elems)
    elems_del = np.empty(elems_len, object)
    edi = 0
    for i in range(elems_len):
        if bool(node_del.intersection(elems[i].n)) == (not excluded):
            elems_del[edi] = elems_del_src[i]
            edi += 1
    elems_del.resize(edi)

    return elems_del


def geod_lspv_rspv(surface_filter, lspv_pt, rspv_pt, tot_pts):
    """Path between LSPV and RSPV markers"""
    pv_pts = [lspv_pt, rspv_pt]
    pts = [None] * 2
    for ind in range(len(pts)):
        pt = ((tot_pts - pv_pts[ind]) ** 2).sum(1)
        pts[ind] = find_index(pt, min)

    dijkstra = vtk.vtkDijkstraGraphGeodesicPath()
    dijkstra.SetInputData(surface_filter.GetOutput())
    dijkstra.SetEndVertex(pts[0])
    dijkstra.SetStartVertex(pts[1])
    dijkstra.StopWhenEndReachedOn()
    dijkstra.UseScalarWeightsOff()
    dijkstra.Update()

    id_len = dijkstra.GetIdList().GetNumberOfIds()
    indices = np.empty(id_len, int)
    mpts = np.empty((id_len, 3), float)
    for i in range(id_len):
        indices[i] = dijkstra.GetIdList().GetId(id_len - i - 1)
        mpts[i] = tot_pts[indices[i]]

    return mpts, indices


def geod_midpoint(surface_filter, mid_point, nodes_la_pv, increment, tot_pts):
    pointd = ((tot_pts - mid_point) ** 2).sum(1)
    index = find_index(pointd, min)

    dijkstra = vtk.vtkDijkstraGraphGeodesicPath()
    dijkstra.SetInputData(surface_filter.GetOutput())
    dijkstra.StopWhenEndReachedOn()
    dijkstra.UseScalarWeightsOff()
    dijkstra.SetEndVertex(index)

    len_store = []
    for ep in range(0, len(nodes_la_pv), increment):
        dijkstra.SetStartVertex(nodes_la_pv[ep])
        dijkstra.Update()
        len_store.append(dijkstra.GetIdList().GetNumberOfIds())

    index = find_index(len_store, max)
    index = list(range(0, len(nodes_la_pv), increment))[index]
    found_index = nodes_la_pv[index]

    return found_index


def geod_marker_mv(surface_filter, marker, increment, tot_pts, nodes_mv):
    """Path from MV to LAA or FO marker"""
    pointd = ((tot_pts - marker) ** 2).sum(1)
    marker_node = find_index(pointd, min)

    dijkstra = vtk.vtkDijkstraGraphGeodesicPath()
    dijkstra.SetInputData(surface_filter.GetOutput())
    dijkstra.StopWhenEndReachedOn()
    dijkstra.UseScalarWeightsOff()
    dijkstra.SetStartVertex(marker_node)

    len_sp = []
    for sp in range(0, len(nodes_mv), increment):
        dijkstra.SetEndVertex(nodes_mv[sp])
        dijkstra.Update()
        len_sp.append(dijkstra.GetIdList().GetNumberOfIds())

    # now choose the min of these
    index = find_index(len_sp, min)
    index = list(range(0, len(nodes_mv), increment))[index]
    dijkstra.SetEndVertex(nodes_mv[index])
    dijkstra.Update()

    id_len = dijkstra.GetIdList().GetNumberOfIds()
    indices = np.empty(id_len, int)
    for i in range(id_len):
        indices[i] = dijkstra.GetIdList().GetId(id_len - i - 1)

    return indices


def write_vtx(ids, filename, intra):
    with open(filename, "w") as file:
        file.write(str(len(ids)) + "\n")
        file.write("intra\n" if intra else "extra\n")
        for d in ids:
            file.write("".join(str(int(d))) + "\n")


def write_vtx_intra(ids, filename):
    write_vtx(ids, filename, True)


def write_vtx_extra(ids, filename):
    write_vtx(ids, filename, False)


def extract_region_by_marker(polydata, marker):
    """Extracts the largest connected region"""
    surfer = vtk.vtkDataSetSurfaceFilter()
    if vtk.VTK_MAJOR_VERSION < 6:
        surfer.SetInput(polydata)
    else:
        surfer.SetInputData(polydata)
    surfer.Update()
    # NOTE: preventive measures: clean before connectivity filter
    # to avoid artificial regionIds
    # It slices the surface down the middle
    cleaner1 = vtk.vtkCleanPolyData()
    cleaner1.SetInputConnection(surfer.GetOutputPort())
    cleaner1.Update()
    # connectivity filter
    connect = vtk.vtkPolyDataConnectivityFilter()
    connect.SetInputConnection(cleaner1.GetOutputPort())
    connect.SetExtractionModeToClosestPointRegion()
    connect.SetClosestPoint(marker)
    connect.Update()
    # leaves phantom points ....
    cleaner2 = vtk.vtkCleanPolyData()
    cleaner2.SetInputConnection(connect.GetOutputPort())
    cleaner2.Update()
    return cleaner2.GetOutput()


def pv_index_remap(pts_src, pts_pv, np_func):
    dmap = [None] * len(pts_pv)
    for i in range(len(pts_pv)):
        pointd = ((pts_pv[i] - pts_src) ** 2).sum(1)
        dmap[i] = np_func(pointd)

    found_index = find_index(dmap, min)
    return found_index


def pv_boundary(pv_index, elems, pts):
    *_, nodes_pv = find_nodes(elems, pts, pv_index + 1)  # +1 since largest is MV
    pts_pv = pts[nodes_pv]
    return pts_pv, nodes_pv


def vein2circle(nodes, pa_ud, ls_lr, tot_pts, surface, xshift, yshift, r1, r2):
    ud_pv = pa_ud[nodes]
    lr_pv = ls_lr[nodes]

    max_ud_pv_ind = nodes[find_index(ud_pv, max)]
    max_lr_pv_ind = nodes[find_index(lr_pv, max)]
    min_ud_pv_ind = nodes[find_index(ud_pv, min)]
    min_lr_pv_ind = nodes[find_index(lr_pv, min)]

    p = [tot_pts[max_lr_pv_ind], tot_pts[max_ud_pv_ind], tot_pts[min_lr_pv_ind], tot_pts[min_ud_pv_ind]]
    size = len(p)
    nodes_found = [None] * size
    strengths_ls = [None] * size
    strengths_pa = [None] * size
    multiplier = 0.0
    increment = 0.5

    for i in range(size):
        _, nodes_found[i] = geod_lspv_rspv(surface, p[i], p[(i + 1) % size], tot_pts)
        nodes_found[i] = np.flip(nodes_found[i], 0)

        # split as four regions
        theta = np.linspace(np.pi * multiplier, np.pi * (multiplier + increment), len(nodes_found[i]))
        strengths_ls[i] = xshift + r1 * np.cos(theta)
        strengths_pa[i] = yshift + r2 * np.sin(theta)
        multiplier += increment

    nodes_found = np.concatenate(nodes_found)
    strengths_ls = np.concatenate(strengths_ls)
    strengths_pa = np.concatenate(strengths_pa)

    return nodes_found, strengths_pa, strengths_ls, max_ud_pv_ind, min_ud_pv_ind, max_lr_pv_ind, min_lr_pv_ind


def get_iso_points(data, value, data_lr, ll, delta, ylist2d, decision_func):
    ind = 0
    size = len(data)
    nl = np.empty(size, int)
    diff = np.empty(size, float)
    for i in range(size):
        if decision_func(i, data, value, data_lr, ll, delta, ylist2d):
            nl[ind] = i
            diff[ind] = data[i] - value
            ind += 1

    nl.resize(ind)
    diff.resize(ind)
    return nl, diff


def get_iso_points_gt_post_gl(data, value, data_lr, ll, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] < (value + delta) and data[i] > value and ylist2d[i] < 0.5 and data_lr[i] > ll

    return get_iso_points(data, value, data_lr, ll, delta, ylist2d, decision_func)


def get_iso_points_lt_post_gl(data, value, data_lr, ll, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] > (value - delta) and data[i] < value and ylist2d[i] < 0.5 and data_lr[i] > ll

    return get_iso_points(data, value, data_lr, ll, delta, ylist2d, decision_func)


def get_iso_points_gt_post_ll(data, value, data_lr, ll, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] < (value + delta) and data[i] > value and ylist2d[i] < 0.5 and data_lr[i] < ll

    return get_iso_points(data, value, data_lr, ll, delta, ylist2d, decision_func)


def get_iso_points_lt_post_ll(data, value, data_lr, ll, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] > (value - delta) and data[i] < value and ylist2d[i] < 0.5 and data_lr[i] < ll

    return get_iso_points(data, value, data_lr, ll, delta, ylist2d, decision_func)


def get_iso_points_gt_post(data, value, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] < (value + delta) and data[i] > value and ylist2d[i] < 0.5

    return get_iso_points(data, value, None, None, delta, ylist2d, decision_func)


def get_iso_points_lt_post(data, value, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] > (value - delta) and data[i] < value and ylist2d[i] < 0.5

    return get_iso_points(data, value, None, None, delta, ylist2d, decision_func)


def get_iso_points_gt_ant(data, value, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] < (value + delta) and data[i] > value and ylist2d[i] > 0.5

    return get_iso_points(data, value, None, None, delta, ylist2d, decision_func)


def get_iso_points_lt_ant(data, value, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] > (value - delta) and data[i] < value and ylist2d[i] > 0.5

    return get_iso_points(data, value, None, None, delta, ylist2d, decision_func)


def get_iso_points_gt_ant_gl(data, value, data_lr, ll, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] < (value + delta) and data[i] > value and ylist2d[i] > 0.5 and data_lr[i] > ll

    return get_iso_points(data, value, None, None, delta, ylist2d, decision_func)


def get_iso_points_lt_ant_gl(data, value, data_lr, ll, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] > (value - delta) and data[i] < value and ylist2d[i] > 0.5 and data_lr[i] > ll

    return get_iso_points(data, value, None, None, delta, ylist2d, decision_func)


def get_iso_points_gt_ant_ll(data, value, data_lr, ll, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] < (value + delta) and data[i] > value and ylist2d[i] > 0.5 and data_lr[i] < ll

    return get_iso_points(data, value, data_lr, ll, delta, ylist2d, decision_func)


def get_iso_points_lt_ant_ll(data, value, data_lr, ll, delta, ylist2d):
    def decision_func(i, data, value, data_lr, ll, delta, ylist2d):
        return data[i] > (value - delta) and data[i] < value and ylist2d[i] > 0.5 and data_lr[i] < ll

    return get_iso_points(data, value, data_lr, ll, delta, ylist2d, decision_func)


def lspv2line(surface, tot_pts, mid_pt, end_pt, xshift, start_seed, r1):
    """Distribute LS along straight line of superior PV"""
    nodes = [None] * 2
    _, nodes[0] = geod_lspv_rspv(surface, mid_pt, tot_pts[start_seed], tot_pts)
    _, nodes[1] = geod_lspv_rspv(surface, end_pt, mid_pt, tot_pts)

    strengths = [None] * 2
    strengths[0] = np.linspace(-r1, 0, len(nodes[0])) + xshift
    strengths[1] = np.linspace(0, r1, len(nodes[1])) + xshift

    nodes = np.concatenate(nodes)
    strengths = np.concatenate(strengths)

    return nodes, strengths


def rspv2line(surface, xshift, tot_pts, start_seed, r1, mid_pt, end_pt):
    """Distribute LS along straight line of superior PV"""
    nodes, strengths = lspv2line(surface, tot_pts, mid_pt, end_pt, xshift, start_seed, r1)
    nodes = np.flip(nodes, 0)
    return nodes, strengths


def rspv2line_ra(xshift, path, r1):
    strength = np.linspace(r1, -r1, len(path)) + xshift
    return path, strength


def lspv2line_ra(xshift, r1, path):
    nodes, strengths = rspv2line_ra(xshift, path, r1)
    strengths = np.flip(strengths, 0)
    return nodes, strengths


def lspv2circle(surface, tot_pts, mid_pt, end_pt, xshift, start_seed, r1, r2):
    _, nodes, strengths_pa, strengths_ls = rspv2circle(surface, xshift, tot_pts, start_seed, r1, r2, end_pt, mid_pt)
    node_ff = nodes[find_index(abs(strengths_ls - 0.05), min)]
    return nodes, strengths_pa, strengths_ls, node_ff


def rspv2circle(surface, xshift, tot_pts, start_seed, r1, r2, end_pt, mid_pt):
    return rspv2circle_ra(surface, xshift, tot_pts, r1, r2, mid_pt, end_pt, tot_pts[start_seed], -1)


def rspv2circle_ra(surface, xshift, tot_pts, r1, r2, mid_pt, end_pt, start_pt, mult=1):
    nodes_pa = [None] * 2
    _, nodes_pa[0] = geod_lspv_rspv(surface, mid_pt, start_pt, tot_pts)
    _, nodes_pa[1] = geod_lspv_rspv(surface, end_pt, mid_pt, tot_pts)

    # split as four regions
    strengths_ls = [None] * 2
    strengths_pa = [None] * 2

    theta = np.linspace(0, -np.pi / 2 * mult, len(nodes_pa[0]))
    strengths_ls[0] = xshift - r1 * np.cos(theta)
    strengths_pa[0] = 1 + r2 * np.sin(theta) * mult

    theta = np.linspace(-np.pi / 2 * mult, -np.pi * mult, len(nodes_pa[1]))
    strengths_ls[1] = xshift - r1 * np.cos(theta)
    strengths_pa[1] = 1 + r2 * np.sin(theta) * mult

    nodes_pa = np.concatenate(nodes_pa)
    nodes_ls = np.flip(nodes_pa, 0)
    strengths_pa = np.concatenate(strengths_pa)
    strengths_ls = np.concatenate(strengths_ls)

    return nodes_ls, nodes_pa, strengths_pa, strengths_ls


def lspv2circle_ra(surface, tot_pts, end_pt, mid_pt, start_pt, xshift, r1, r2):
    _, nodes, strengths_pa, strengths_ls = rspv2circle_ra(surface, xshift, tot_pts, r1, r2, mid_pt, end_pt, start_pt)
    node_ff = nodes[find_index(abs(strengths_ls - 0.05), min)]
    return nodes, strengths_pa, strengths_ls, node_ff


def triangle_calc(pts, elems, no):
    mpts = np.zeros((len(elems), 3), float)
    for i in range(len(elems)):
        mpts[i] = pts[elems[i].n[no]]

    return mpts


def mp_calc_ele(pts, elems, n):
    """Calculate midpoints"""
    mpts = np.zeros((len(elems), 3), float)
    for i in range(len(elems)):
        for j in range(n):
            mpts[i] += pts[elems[i].n[j]]
        mpts[i] /= n

    return mpts


def mp_calc_ele_3(pts, elems):
    return mp_calc_ele(pts, elems, 3)


def mp_calc_ele_4(pts, elems):
    return mp_calc_ele(pts, elems, 4)


def barycentric_coord_v1(p, a, b, c):
    v0 = c - a
    v1 = b - a
    v2 = p - a
    d00 = np.dot(v0, v0)
    d01 = np.dot(v0, v1)
    d02 = np.dot(v0, v2)
    d11 = np.dot(v1, v1)
    d12 = np.dot(v1, v2)
    denom = d00 * d11 - d01 * d01
    if denom == 0:
        u = 1
        v = 0
        w = 0
    else:
        u = (d11 * d02 - d01 * d12) / denom
        v = (d00 * d12 - d01 * d02) / denom
        w = 1.0 - u - v

    return u, v, w


def barycentric_coord_v2(p, a, b, c, d):
    vap = p - a
    vbp = p - b
    vab = b - a
    vac = c - a
    vad = d - a
    vbc = c - b
    vbd = d - b

    # <a, b x c>

    va6 = np.linalg.det(np.dstack([vbp, vbd, vbc]))
    vb6 = np.linalg.det(np.dstack([vap, vac, vad]))
    vc6 = np.linalg.det(np.dstack([vap, vad, vab]))
    vd6 = np.linalg.det(np.dstack([vap, vab, vac]))
    denom = np.linalg.det(np.dstack([vab, vac, vad]))

    if denom == 0:
        u = 1
        v = 0
        w = 0
        x = 0
    else:
        u = va6 / denom
        v = vb6 / denom
        w = vc6 / denom
        x = vd6 / denom

    return u, v, w, x


def relabel_elems_in_anteriormesh(elems, ylist2d, svc_label, ivc_label):
    """Relabel SVC & IVC if they are in AnteriorMesh i.e. >0.5"""
    ylist2d_elements = to_ele(elems, ylist2d)
    labels = np.empty(len(elems), int)
    for i in range(len(elems)):
        labels[i] = elems[i].reg
        if labels[i] == svc_label or labels[i] == ivc_label:
            if ylist2d_elements[i][0] >= 0.5:
                labels[i] = 1

    return labels


def get_sub_surface(base_dir, pts, elems, node_del, elems_del_src=None, return_pts=False):
    elems_del = element2delete(elems, node_del, elems_del_src)
    tmp_file = base_dir + "tmp.vtk"
    write_vtk(pts, elems_del, None, None, tmp_file)
    surface = get_vtk_from_file(tmp_file)
    os.remove(tmp_file)
    surface = convert_unstructureddata_to_polydata(surface)
    tot_pts, surface = read_vtk_surface(surface)
    return tot_pts if return_pts else surface
