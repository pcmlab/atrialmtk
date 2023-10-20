#!/usr/bin/env python

import numpy as np
from vtk.util import numpy_support
from xml.etree.ElementTree import ElementTree
from collections import namedtuple


class vtuWriter:
    """
    USAGE:
    vtk_writer = vtuWriter()
    vtk_writer.writeVTU("filename.vtu", pts, all other optional arguments...)
    vtk_writer.writePVD("filename.pvd")
    """

    def __init__(self):
        self.fileNames = []

    def toString(self, a):
        if len(np.shape(a)) > 1:
            return " ".join(repr(num) for row in a for num in row)
        else:
            return " ".join(repr(num) for num in a)

    def writeVTU(self, fileName, nodes, cells=[], celltype=[], celldata={}, pointdata={}, format=[]):
        """
        ARGUMENTS:
        fileName        file name and/or path/filename
        nodes           array of coordinates of points
        cells           optional array of particle connectivity
        celltype        optional array of type of connectivity for cells
        .
        """
        import xml.dom.minidom

        # import xml.dom.ext # python 2.5 and later

        # Document and root element
        doc = xml.dom.minidom.Document()
        root_element = doc.createElementNS("VTK", "VTKFile")
        root_element.setAttribute("type", "UnstructuredGrid")
        root_element.setAttribute("version", "0.1")
        root_element.setAttribute("byte_order", "LittleEndian")
        doc.appendChild(root_element)

        if len(format) == 0:
            format = "ascii"

        # Unstructured grid element
        unstructuredGrid = doc.createElementNS("VTK", "UnstructuredGrid")
        root_element.appendChild(unstructuredGrid)

        # Piece 0 (only one)
        piece = doc.createElementNS("VTK", "Piece")
        piece.setAttribute("NumberOfPoints", str(len(nodes)))
        piece.setAttribute("NumberOfCells", str(len(cells)))
        unstructuredGrid.appendChild(piece)

        # Points
        points = doc.createElementNS("VTK", "Points")
        piece.appendChild(points)

        # Point location data
        point_coords = doc.createElementNS("VTK", "DataArray")
        point_coords.setAttribute("type", "Float32")
        point_coords.setAttribute("format", format)
        point_coords.setAttribute("NumberOfComponents", "3")
        points.appendChild(point_coords)

        string = self.toString(np.array(nodes))
        point_coords_data = doc.createTextNode(string)
        point_coords.appendChild(point_coords_data)

        # Cells
        elems = doc.createElementNS("VTK", "Cells")
        piece.appendChild(elems)

        # Cell connectivity
        cell_connectivity = doc.createElementNS("VTK", "DataArray")
        cell_connectivity.setAttribute("type", "Int32")
        cell_connectivity.setAttribute("Name", "connectivity")
        cell_connectivity.setAttribute("format", format)
        elems.appendChild(cell_connectivity)
        if len(cells) > 0:
            string = self.toString(np.array(np.hstack(cells)))
            cells_connectivity_data = doc.createTextNode(string)
        else:
            cells_connectivity_data = doc.createTextNode("")
        cell_connectivity.appendChild(cells_connectivity_data)

        # Cell offsets
        cell_offsets = doc.createElementNS("VTK", "DataArray")
        cell_offsets.setAttribute("type", "Int32")
        cell_offsets.setAttribute("Name", "offsets")
        cell_offsets.setAttribute("format", format)
        elems.appendChild(cell_offsets)

        if len(cells) > 0:
            count = 0
            offset = []
            for e in cells:
                count = count + len(e)
                offset.append(count)
            string = self.toString(offset)
            offsets = doc.createTextNode(string)
        else:
            offsets = doc.createTextNode("0")
        cell_offsets.appendChild(offsets)

        # Cell type
        cell_types = doc.createElementNS("VTK", "DataArray")
        cell_types.setAttribute("type", "UInt8")
        cell_types.setAttribute("Name", "types")
        cell_types.setAttribute("format", format)
        elems.appendChild(cell_types)

        if len(celltype) > 0:
            string = self.toString(celltype)
            cell_types_data = doc.createTextNode(string)
        else:
            if len(cells) > 0:
                ctype = []
                for e in cells:
                    if len(e) == 1:
                        ctype.append(1)
                    if len(e) == 2:
                        ctype.append(3)
                    if len(e) == 3:
                        ctype.append(5)
                    if len(e) == 4:
                        ctype.append(10)
                    if len(e) == 8:
                        ctype.append(12)
                string = self.toString(ctype)
                cell_types_data = doc.createTextNode(string)
            else:
                cell_types_data = doc.createTextNode("")
        cell_types.appendChild(cell_types_data)

        # PointData
        point_data = doc.createElementNS("VTK", "PointData")
        piece.appendChild(point_data)
        for (name, value) in list(pointdata.items()):
            data = doc.createElementNS("VTK", "DataArray")
            data.setAttribute("Name", name)
            if len(np.shape(value)) == 2:
                data.setAttribute("NumberOfComponents", "3")
            data.setAttribute("type", "Float32")
            data.setAttribute("format", format)
            point_data.appendChild(data)
            string = self.toString(np.array(value))
            xdata = doc.createTextNode(string)
            data.appendChild(xdata)

        # CellData
        cell_data = doc.createElementNS("VTK", "CellData")
        piece.appendChild(cell_data)
        for (name, value) in list(celldata.items()):
            cdata = doc.createElementNS("VTK", "DataArray")
            cdata.setAttribute("Name", name)
            if len(np.shape(value)) == 2:
                cdata.setAttribute("NumberOfComponents", "3")
            cdata.setAttribute("type", "Float32")
            cdata.setAttribute("format", format)
            cell_data.appendChild(cdata)
            string = self.toString(np.array(value))
            xdata = doc.createTextNode(string)
            cdata.appendChild(xdata)

        # Write to file and exit
        outFile = open(fileName, "w")
        #        xml.dom.ext.PrettyPrint(doc, file)
        doc.writexml(outFile, newl="\n")
        outFile.close()
        self.fileNames.append(fileName)

    def writePVD(self, fileName):
        outFile = open(fileName, "w")
        import xml.dom.minidom

        pvd = xml.dom.minidom.Document()
        pvd_root = pvd.createElementNS("VTK", "VTKFile")
        pvd_root.setAttribute("type", "Collection")
        pvd_root.setAttribute("version", "0.1")
        pvd_root.setAttribute("byte_order", "LittleEndian")
        pvd.appendChild(pvd_root)

        collection = pvd.createElementNS("VTK", "Collection")
        pvd_root.appendChild(collection)

        for i in range(len(self.fileNames)):
            dataSet = pvd.createElementNS("VTK", "DataSet")
            dataSet.setAttribute("timestep", str(i))
            dataSet.setAttribute("group", "")
            dataSet.setAttribute("part", "0")
            dataSet.setAttribute("file", str(self.fileNames[i]))
            collection.appendChild(dataSet)

        outFile = open(fileName, "w")
        pvd.writexml(outFile, newl="\n")
        outFile.close()


DataArray = namedtuple("DataArray", ["name", "data"])
Piece = namedtuple("Piece", ["points", "cells", "point_data", "cell_data"])
Cell = namedtuple("Cell", ["connectivity", "offsets", "types"])


def parse_data_array(data_array):
    vtk_to_np_type = {
        "UInt8": "uint8",
        "UInt16": "uint16",
        "Int32": "int32",
        "Int64": "int64",
        "Float32": "float32",
        "Float64": "float64",
    }
    assert data_array.tag == "DataArray"
    name = data_array.get("Name")
    noc = int(data_array.get("NumberOfComponents", default="1"))
    type = data_array.get("type")
    # format = data_array.get("format", default="ascii")
    # assert (format == 'ascii')
    data = [float(val) for val in data_array.text.split()]
    array = np.array(data, dtype=vtk_to_np_type[type])
    components = []
    for i in range(noc):
        components.append(array[i::noc])
    return DataArray(name, components)


def parse_all_data_array(element):
    d = {}
    if element.findall("DataArray"):
        for data_array in element.findall("DataArray"):
            data = parse_data_array(data_array)
            d[data.name] = data.data
    return d


def parse_points(points):
    pts = None
    if points:
        data = parse_data_array(points.find("DataArray"))
        pts = [list(i) for i in zip(*data.data)]
    return pts


def parse_cells(cells):
    elems = None
    if cells:
        d = parse_all_data_array(cells)
        for (key, value) in list(d.items()):
            d[key] = value[0]
        elems = Cell(**d)
    return elems


def parse_piece(piece):
    assert piece.tag == "Piece"
    points = parse_points(piece.find("Points"))
    cells = parse_cells(piece.find("Cells"))
    point_data = parse_all_data_array(piece.find("PointData"))
    cell_data = parse_all_data_array(piece.find("CellData"))
    #    cell_count  = int (piece.get ('NumberOfCells'))
    #    point_count = int (piece.get ('NumberOfPoints'))
    #    assert (cell_count == len (cells.offsets))
    #    assert (cell_count == len (cells.types))
    return Piece(points=points, cells=cells, point_data=point_data, cell_data=cell_data)


def vtuReader(fileName):
    """
    USAGE:
    pts, elem, types, fibers, data = vtuReader(vtufname)
    """
    tree = ElementTree()
    tree.parse(fileName)
    root = tree.getroot()

    if root.tag != "VTKFile":
        raise Exception("Root element must be 'VTKFile'")

    type = root.get("type")
    version = root.get("version")
    byte_order = root.get("byte_order")
    if type is None or version is None or byte_order is None:
        raise Exception("Missing attribute!")

    grid = root.find(type)
    piece = parse_piece(grid.find("Piece"))

    elems = []
    elems.append(piece.cells.connectivity[0: piece.cells.offsets[0]])
    for i in range(1, len(piece.cells.offsets)):
        elems.append(piece.cells.connectivity[piece.cells.offsets[i - 1]: piece.cells.offsets[i]])

    data = None
    for (name, value) in list(piece.point_data.items()):
        if len(value) == 1:
            data = value

    fibers = None
    for (name, value) in list(piece.cell_data.items()):
        if name == "fibers" and len(value) == 3:
            fibers = [list(i) for i in zip(*value)]

    return piece.points, elems, piece.cells.types, fibers, data


def vtiReader(fileName):
    """
    USAGE:
    data, extent, origin, spacing = vtiReader(vtufname)
    """
    tree = ElementTree()
    tree.parse(fileName)
    root = tree.getroot()

    if root.tag != "VTKFile":
        raise Exception("Root element must be 'VTKFile'")

    type = root.get("type")
    version = root.get("version")
    byte_order = root.get("byte_order")
    if type != "ImageData" or version is None or byte_order is None:
        raise Exception("Missing attribute!")

    grid = root.find(type)
    extent = list(map(int, grid.get("WholeExtent").split()))
    origin = list(map(float, grid.get("Origin").split()))
    spacing = list(map(float, grid.get("Spacing").split()))
    piece = parse_piece(grid.find("Piece"))

    return piece.point_data, extent, origin, spacing


def vtkImageToNumPy(image, pixelDims):
    pointData = image.GetPointData()
    arrayData = pointData.GetArray(0)
    ArrayDicom = numpy_support.vtk_to_numpy(arrayData)
    ArrayDicom = ArrayDicom.reshape(pixelDims, order="F")
    return ArrayDicom
