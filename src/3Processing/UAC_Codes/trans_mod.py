#!/usr/bin/env python3
#
# This script converts between 2 model formats
# Formats are guessed from the suffix.
# All models are first converted to CARP format which is 0-index based

import sys, gzip, pdb, time, shutil,re
import os.path
import argparse
import numpy as np
import scipy.io as io
import struct
import meshio
from collections import defaultdict
from argparse import RawTextHelpFormatter

# custom formats - formats provided in this file
custom_formats = [ 'carp', 'carpbin', 'CARTO' ]

# number of vertices per element
eleminfo    = { 'Ln':2, 'Tr':3, 'Qd':4, 'Tt':4,  'Py':5,  'Pr':6, 'Hx':8, 'Pt':1 }
nnodes_elem = { 'line':2, 'triangle':3, 'quad':4, 'tetra':4, 'pyramid':5, 'wedge':6, 'hexahedron':8 }

#extensions of known formats
format_exts={ 'ply':'ply', 'vtk':'vtk', 'vtu':'vtu', 'obj':'obj', 'off':'off', 'msh':'gmsh22', 'mat':'matlab',
        'pts_t':'carp_aux', 'mesh':'medit', 'meshb':'medit', 'pts':'carp',
        'inp':'abaqus', 'avs':'avsucd', 'cgns':'cgns', 'xml':'dolphin-xml', 'exo':'exodus', 
        'f3grid':'flac3d', 'h5m':'h5m', 'mdpa':'mdpa', 'vol':'netgen', 
        'vol.gz':'netgen', 'stl':'stl', 'su2':'su2', 'ugrid':'ugrid',
             'wkt':'wkt', 'bpts':'carpbin', 'bwmesh':'CARTO' }

carp_predef = [ 'lon', 'sheet', 'region' ]

class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

#
# class to read/write binary data
#
class  BinaryFile :

    def __init__ (self, filename, mode ) :
        self.file = open( filename, mode+'b' )
        self.mode = mode

    def __enter__ (self ) :
        return self
    
    def __exit__ (self,a,b,c):
        self.file.close()

    def get( self, fmt ) :
        sz = struct.calcsize( fmt )
        buf = self.file.read( sz )
        return struct.unpack_from( fmt, buf )

    def set_position( self, pos ) :
        self.file.seek( pos )

    def put( self, fmt, *data ) :
        d = struct.pack( fmt, *data )
        self.file.write( d )

    # output a string
    def write( self, s ):
        self.file.write( s.encode() )

    def position( self ) :
        return self.file.ftell()

# Class to read in one word at a time ignoring line breaks
# A single character can be specified for beginning of comments
# which last until the end of the line. Use None to disable this. '#' is the default
class WordReader :

    # constructor
    # comstart is the character signifying start of a comment, can be None
    def  __init__(self, file, comstart='#'):
        self.file= file
        self.cword = 0
        self.line = []
        self.comstart = comstart
    
    # get next word as a string
    def get_next(self) :
        self.cword += 1
        if self.cword >= len(self.line) :
            while True :
                l = self.file.readline()
                if not l :
                    self.close()
                    return l
                if ',' in l :
                    l = re.sub(',',' ',l) # ?????
                if self.comstart :
                    self.line = l.split(self.comstart)[0].split()  #filter comments
                else :
                    self.line = l.split()
                if len(self.line) > 0 :
                    break
            self.cword = 0
        return self.line[self.cword]

    # get last word returned
    def last( self ) : 
        if self.closed() :
            return None
        else:
            return self.line[self.cword]

    # get next n words as strings, list for n>1
    def __next__( self, n=1 ):
        if n==1 :
            return self.get_next()
        else :
            ret = []
            for i in range(n) :
                ret.append( self.get_next() )
            return ret

    # get next n words as strings, list for n>1
    def next( self, n=1 ):
        return self.__next__(n)

    # get next n words as integers, list for n>1
    def nextI (self,n=1) :
        if n==1 :
            return int(next(self))
        else :
            ret=[0]*n;
            for i in range(n) :
                ret[i] = int(self.get_next())
            return ret

    # get next n words as floats, list for n>1
    def nextF (self,n=1) :
        if n==1 :
            return float(next(self))
        else :
            ret=[0.]*n;
            for i in range(n) :
                ret[i] = float(self.get_next())
            return ret

    # ignore rest of line
    def clear_line( self ) :
        self.cword = len(self.line)

    # read next line possibly ignoring commented and blank lines
    def next_line( self, ignoreComments=False ) :
        if ignoreComments :
            self.get_next_info_line()
        else :
            return self.file.readline()

    # skip the next n words
    def skip( self, n=1 ) :
        for i in range(n) :
            next(self)

   # skip the next n lines
    def skipline( self,n=1, ignoreComments=False ) :
        if ignoreComments:
            for i in range(n) :
                self.get_next_info_line()
        else :
            for i in range(n) :
                self.file.readline() 

    # get next line which is not blank or commented out
    def get_next_info_line( self ) :
        while True :
            line = self.file.readline() 
            if self.comstart : line = line.split(self.comstart)[0]
            if len( line.split() ) : return line

    # read file until the word matches target
    def until( self, target ):
        while self.next() != target :
            if self.last() is None : return False
        return True

    # close file
    def close( self ) :
        self.file.close()

    # is closed?
    def closed( self ):
        return self.file.closed


def format_pt( p, prec ):
    return f'{p[0]:.{prec}} {p[1]:.{prec}} {p[2]:.{prec}}'
    

        
################################################
# mat file format (Matlab file)
################################################

# map CARP <-> mat types
mat_elemtypes = { 'line':2, 'triangle':3, 'quad':8, 'tetra':4}
mat2mesh_etype = dict([(v,k) for (k,v) in list(mat_elemtypes.items())])

def read_mat( matfname ) :
    b    = io.loadmat(matfname)
    
    pts  = b['pts'].tolist() if 'pts' in b else None

    elems   = defaultdict(list)
    cdata   = defaultdict(list)
    lon     = defaultdict(list)
    
    if 'elems' in b:
        fibs = b['fibers'].tolist() if 'fibers' in b else None
        data = b['data'].tolist()   if 'data' in b else None
        for i,e in enumerate(b['elems']) :
          et = mat2carp_etype[len(e)]
          elems[et].append([x-1 for x in e])
          if fibs is not None:
              lon[et].append(fibers[i])
          if data is not None: 
              cdata[et].append(data[i])
    
    ptdata  = b['pointdata'].tolist() if 'pointdata' in b else None

    cells    = []
    cellData = defaultdict(list)

    for k in met_elemtypes:
        if k in elems :
            cells.append( (k, elems[k]) )
            if fibs:
                cellData['lon'].append( lon[k] )
            if data:
                cellData['data'].append( cdata[k] )

    return meshio.Mesh( pts, elems, cell_data=cellData, point_data=ptdata )


def write_mat (mesh, outputMesh):
    a = {}
    a['pts'] = mesh.points
    ele = []
    lon = []
    fibres = mesh.cell_data['lon']   if 'lon'  in mesh.cell_data else None
    data   = mesh.point_data['data'] if 'data' in mesh.cell_data else None
    for cb in mesh.cells:
        for i,e in enumerate(cb.data):
            ele.append([x+1 for x in e])
            if fibres: lon.append( fibres[i] )
    a['elems'] = ele
    if (fibers != None):
      a['fibers'] = lon
    if (data != None):
      a['pointdata'] = data
    io.savemat(outputMesh,a)

################################################
# CARP format
################################################
carp2meshio = { 'Ln':'line', 'Tr':'triangle', 'Qd':'quad', 'Tt':'tetra',
        'Pr':'wedge', 'Py':'pyramid', 'Hx':'hexahedron', 'Pt':'vertex'         }

meshio2carp = dict([(v,k) for (k,v) in list(carp2meshio.items())])

def read_carp( carpMesh ) :
    ptsFile  = carpMesh + '.pts'
    elemFile = carpMesh + '.elem'
    lonFile  = carpMesh + '.lon'

    with open( ptsFile, 'r' ) as ptin :
        ptfile = ptin.readlines()
    pts = []
    for l in ptfile[1:]:
        l = l.strip()
        if l: pts.append( [float(a) for a in l.split()[0:3]] )

    elems   = defaultdict(list)
    fibres  = defaultdict(list)
    regions = defaultdict(list)
    sheet   = defaultdict(list)
    with open( elemFile, 'r' ) as elemin :
        elemf = elemin.readlines()
    with open( lonFile, 'r' ) as lonin :
        lonf = lonin.readlines()
    nr = len(lonf[0].strip().split())
    if nr < 3:
        ndir = lonf[0].strip().split()[0]
        del lonf[0]
    elif nr == 3 :
        ndir = 1
    else:
        ndir = 2
    if elemf[-1].strip() == '': del elemf[-1]
    if lonf[-1].strip()  == '': del lonf[-1]

    if abs(len(elemf)-len(lonf)) > 1 :
        print( "Elem and lon file mismatch" )
        sys.exit(1)

    for e,l in zip(elemf[1:],lonf):
        words = e.strip().split()
        etype = carp2meshio[words[0]]
        indx  = np.asarray(words[1:eleminfo[words[0]]+1], int)
        elems[etype].append(indx)
        axes =  np.asarray(l.split(), float)
        fibres[etype].append( axes[:3] )
        if ndir == 2: sheet[etype].append( axes[3:6] )
        if len(words) > len(indx)+1 :
            regions[etype].append( int(words[-1]) )

    if 'hexahedron' in elems:
        for ind in elems['hexahedron']:
            ind[5],ind[7] = ind[7],ind[5]
    if 'wedge' in elems:
        for ind in elems['wedge']:
            ind[3],ind[4] = ind[4],ind[3]

    cells      = []
    cellData   = defaultdict(list)
    region_dat = np.any([ len(v)>0 for k,v in regions.items() ])
    
    for k in meshio2carp :
        if k in elems :
            cells.append( (k, elems[k]) )
            cellData['lon'].append( fibres[k] )
            if sheet: cellData['sheet'].append( sheet[k] )
            if region_dat :
                cellData['region'].append( regions[k] )

    return meshio.Mesh( pts, cells, cell_data=cellData )


def write_carp( mesh, base ) :

    ext = os.path.splitext(base)[1]
    if ext=='.pts'  : base=base[:-4]

    #point filehttps://www.programiz.com/python-programming/methods/built-in/eval
    with open( base+'.pts', 'w' ) as ptout:
        ptout.write( f'{np.shape(mesh.points)[0]}\n' )
        for p in mesh.points :
            ptout.write( f'{p[0]} {p[1]} {p[2]}\n' )

    # point data
    for dt in mesh.point_data :
        if dt not in carp_predef :
            dshape = np.shape(mesh.point_data[dt])
            ndat   = 1 if len(dshape)==1 else dshape[1]
            ext    = 'dat' if ndat == 1 else 'vec'
            fformat = '{0}' if ndat == 1 else '{0[0]}'
            for i in range(1,ndat):
                fformat += f' {{0[{i}]}}'
            fformat += '\n'
            print(f'Writing data file: {base}-{dt}.{ext}')
            with open( f'{base}-{dt}.{ext}', 'w' ) as datout:
                for d in mesh.point_data[dt]:
                    datout.write( fformat.format(d) )

    # count elements
    numele = 0
    for cellblock in mesh.cells:
        if cellblock.type != 'vertex' : 
            numele += np.shape(cellblock.data)[0]
    regdat = mesh.cell_data.get('region')

    # write elem file
    with open( base+'.elem', 'w' ) as eout:
        eout.write( f'{numele}\n' )
        for cb, cellblock in enumerate(mesh.cells):
            if cellblock.type == 'vertex': continue
            fmtstr = meshio2carp[cellblock.type]
            nodes = [ *range(eleminfo[meshio2carp[cellblock.type]]) ]
            if   cellblock.type == 'hexahedron' : nodes[5],nodes[7]=nodes[7],nodes[5]
            elif cellblock.type == 'wedge'      : nodes[3],nodes[4]=nodes[4],nodes[3]
            for i in nodes:
                fmtstr += ' {0['+str(i)+']}'
            if regdat : cb_reg = regdat[cb]
            for i,e in enumerate(cellblock.data):
                eout.write(fmtstr.format(e))
                if regdat:
                    eout.write(f' {cb_reg[i]}')
                eout.write('\n')

    # write lon file
    with open( base+'.lon', 'w' ) as lonout:
        londata = mesh.cell_data.get('lon')
        if londata:
            for cb in londata:
                for f in cb:
                    lonout.write( '{0[0]} {0[1]} {0[2]}\n'.format(f) )
        else:
            for i in range(numele) :
                lonout.write('1 0 0\n')
 

################################################
# CARP binary format
################################################
carpb2meshio = { 0:'tetra', 1:'hexahedron', 2:'octahedron', 3:'pyramid', 4:'wedge',
                 5:'quad', 6:'triangle', 7:'line' }
meshio2carpb = dict([(v,k) for (k,v) in list(carpb2meshio.items())])

def read_carpb( base, data=None ) :
    pts = []
    with BinaryFile( base+'.bpts', "r" ) as ptsf:
        hdr = ptsf.get('1024s')
        npts = int(hdr[0].strip().split()[0])
        for i in range(npts) :
            x,y,z = ptsf.get('fff')
            pts.append( [x,y,z] )

    elems   = defaultdict(list)
    regions = defaultdict(list) 
    lon     = defaultdict(list) 
    sheet   = defaultdict(list)
    with BinaryFile( base+'.belem', "r" ) as elef, BinaryFile( base+'.blon', 'r') as lonf:
        hdr   = elef.get('1024s')
        nelem = int(hdr[0].strip().split()[0])
        hdr   = lonf.get('1024s')
        ndir  = int(hdr[0].strip().split()[0])
        for i in range(nelem) :
            x, = elef.get('i')
            et = carpb2meshio[x]
            a  = elef.get('i'*(nnodes_elem[et]+1))
            elems[et].append( a[0:nnodes_elem[et]] )
            regions[et].append( a[-1] )
            b  = lonf.get('f'*3*ndir)
            lon[et].append(b[0:3])
            if ndir == 2 : sheet[et].append(b[3:])
        
    if 'hexahedron' in elems:
        for ind in elems['hexahedron']:
            ind[5],ind[7] = ind[7],ind[5]
    if 'wedge' in elems:
        for ind in elems['wedge']:
            ind[3],ind[4] = ind[4],ind[3]

    cells    = []
    cellData = defaultdict(list)
    
    for k in meshio2carp :
        if k in elems :
            cells.append( (k, elems[k]) )
            cellData['lon'].append( lon[k] )
            cellData['region'].append( regions[k] )
            if sheet : cellData['sheet'].append(sheet[k])

    mesh = meshio.Mesh( pts, cells, cell_data=cellData )
    return mesh


def write_carpb( mesh, base ) :
    
    ext = os.path.splitext(base)[1]
    if ext=='.pts'  : base=base[:-4]
    if ext=='.bpts' : base=base[:-5]

    with BinaryFile( base+'.bpts', 'w' ) as ptout:
        ptout.write( f'{np.shape(mesh.points)[0]:<1024}' )
        for p in mesh.points :
            ptout.put( 'fff', *p[0:3] )

    numele = 0
    for cellblock in mesh.cells:
        numele += np.shape(cellblock.data)[0]
    regdat = mesh.cell_data.get('region')

    with BinaryFile( base+'.belem', 'w' ) as eout:
        eout.write( f'{numele:<1024}' )
        for cb, cellblock in enumerate(mesh.cells):
            cet    = meshio2carpb[cellblock.type]
            fmtstr = 'i'*(nnodes_elem[cellblock.type]+2)
            if regdat : cb_reg = regdat[cb]
            for i,e in enumerate(cellblock.data):
                if   cellblock.type == 'hexahedron' : e[5],e[7]=e[7],e[5]
                elif cellblock.type == 'wedge' :      e[3],e[4]=e[4],e[3]
                eout.put(fmtstr, cet, *e, cb_reg[i] if regdat else 1)

    with BinaryFile( base+'.blon', 'w' ) as lonout:
        londata = mesh.cell_data.get('lon')
        numfib = 1
        lonout.write(f'{numfib:<1024}')
        if londata:
            for cb in londata:
                for f in cb:
                    lonout.put( 'fff', *f[0:3] )
        else:
            for i in range(numele) :
                lonout.put( 'fff', 1., 0, 0 )
 
    # point data
    for dt in mesh.point_data :
        if dt not in carp_predef :
            write_array_igb( f'{base}-{dt}.igb', mesh.point_data[dt] )


################################################
# aux format -- suffix .pts_t
################################################

def read_aux( auxfname ) :
    ptstFile  = auxfname + '.pts_t'
    dattFile  = auxfname + '.dat_t'
    row,ptst  = read_array_pts_t (ptstFile)
    row,elemt = read_array_dat_t (dattFile)
    return ptst, elemt, row


def read_pts_t (ptstFile):
    """
    Purpose: Function to read a .pts_t file from CARP, where the first
    line contains header, second line number of time frames and the remaining pts list
    for each time frame.
    """
    try:
        f = open(ptstFile)
        line = f.readline()
    except IOError:
        print(("IOError: File %s not found." % ptstFile))
        sys.exit(-1)
    
    if line[:10] != '## Format ':
        raise Exception(ptstFile+' not recognized as a pts_t file')

    wread = WordReader( f, None )
    nt = wread.nextI()
    row = []
    col = []
    row.append(0)
    for i in range(nt) :
        if wread.get_next()[:1] != '#':
          raise Exception(ptstFile+' at time # '+str(i+1)+' not recognized')
        wread.clear_line()
        n = wread.nextI()
        row.append(row[-1] + n)
        for j in range(n):
            col.append(wread.nextF(3))
            wread.clear_line()
    return row, col

# end of read_array_pts_t

def read_array_dat_t (dattFile):
    try:
        f = open(dattFile)
        line = f.readline()
    except IOError:
        print(("IOError: File %s not found." % dattFile))
        sys.exit(-1)
    
    if line[:10] != '## Format ':
        raise Exception(dattFile+' not recognized as a dat_t file')

    wread = WordReader( f, None )
    nt = wread.nextI()
    row = []
    col = []
    row.append(0)
    for i in range(nt) :
        if wread.get_next()[:1] != '#':
          raise Exception(dattFile+'time'+str(i)+' # not recognized')
        n = wread.nextI()
        row.append(row[-1] + n)
        for j in range(n):
            col.append(wread.nextI())
    return row, col

def write_pts_t ( pts, elems, fibers, data, base ) :
    base = base[:-6]
    with open( base+'.pts_t', 'w' ) as ptout:
        ptout.write( f'1\n{len(pts)}\n' )
        for p in pts :
            ptout.write( '{0[0]} {0[1]} {0[2]}\n'.format(p) )

    with open( base+'.elem_t', 'w' ) as eout:
        eout.write( f'1\n{len(elems)}\n' )
        for e in elems :
            eout.write( ' '.join(e)+'\n')

    if data :
        print( "Ignoring IGB data for now", out=stderr )

################################################
# CARTO format -- suffix .bwmesh
################################################

def read_CARTO( file, data=None ):

    with open(file,'r',errors='ignore') as infile:
        meshin = WordReader( infile, ';' )

        meshin.until('NumVertex')
        meshin.skip()
        nVtx = meshin.nextI();

        meshin.until('NumTriangle')
        meshin.skip()
        nTri = meshin.nextI()

        meshin.until('[VerticesSection]')
        pts = []
        for p in range(nVtx):
            meshin.skip(2)
            pts.append( meshin.nextF(3) )
            meshin.clear_line()
        
        meshin.until('[TrianglesSection]')
        elems   = defaultdict(list)
        fibres  = defaultdict(list)
        regions = defaultdict(list)
        axes    = np.asarray( [1., 0., 0.] )
        for e in range(nTri):
            meshin.skip(2)
            indx = np.asarray(meshin.nextI(3))
            elems['triangle'].append(indx)
            fibres['triangle'].append(axes)
            meshin.skip(3)
            regions['triangle'].append(meshin.nextI())
            meshin.clear_line()

    cells    = []
    cellData = defaultdict(list)
    for k in meshio2carp :
        if k in elems :
            cells.append( (k, elems[k]) )
            cellData['lon'].append( fibres[k] )
            cellData['region'].append( regions[k] )
                                   
    mesh = meshio.Mesh( pts, cells, cell_data=cellData )
    return mesh

###########################
# End Formats
###########################

# remove unused points and renumber all elements
def rm_unused_pts( mesh ):
    # identify used nodes
    nmap    = [-1]*len(mesh.points)
    for cb in mesh.cells :
        for elem in cb.data:
            for n in elem: nmap[n]=1

    # construct map of old to new indices
    red_pts = []
    for i,p in enumerate(mesh.points):
        if nmap[i] == 1 :
            nmap[i] = len(red_pts)
            red_pts.append( p )
    
    # output list of kept nodes
    with open( "kept_pts.dat", "w" ) as kpout:
        for i,n in enumerate(nmap):
            if n != -1 : print( i, file=kpout )

    if len(red_pts) == len(mesh.points): return mesh

    # map element node numbers
    new_cells = []
    for cb in mesh.cells:
        new_cells.append( (cb.type,list()) )
        for elem in cb.data:
            new_cells[-1][1].append([ nmap[n] for n in elem ])

    # remove unused nodal data
    for data in mesh.point_data :
        mesh.point_data[data] = [d for d,m in zip(mesh.point_data[data],nmap) if m!=-1]

    return meshio.Mesh( red_pts, new_cells, cell_data=mesh.cell_data, point_data=mesh.point_data )


# remove all regions not in the provided regs
def select_regions( mesh, regs ):
    if 'region' not in mesh.cell_data:
        print( 'No regions specified in model: cannot select by region' )
        return mesh

    dellist = []
    elim    = []
    for ereg in mesh.cell_data['region']:
        dellist.append(list())
        for i,r in enumerate(ereg):
            if r not in regs : dellist[-1].append(i)
        elim.append( len(dellist[-1]) == len(ereg) )

    cells = list()
    for i,cb in enumerate(mesh.cells):
        if not elim[i] :
            cells.append( (cb.type, np.delete(cb.data,dellist[i],0)) )

    cellData = defaultdict(list)
    for dt, datalst in mesh.cell_data.items():
        for i,dlst in enumerate(dellist):
            if not elim[i] :
                cellData[dt].append(np.delete(datalst[i], dlst, 0))

    return meshio.Mesh( mesh.points, cells, cell_data=cellData, point_data=mesh.point_data )

# write an IGB file
def write_array_igb( filename, data ):
        npts = np.shape(data)[0]
        ndat = 1 if len(np.shape(data))==1 else np.shape(data)[1]
        if data.dtype == float:
            datype = 'float' if ndat == 1 else f'vec{ndat}f'
            fformat = 'f'*ndat
        else:
            datype = 'int'   if ndat == 1 else f'vec{ndat}d'
            fformat = 'i'*ndat
        hdr = f'x:{npts} y:1 z:1 t:{ndat} type:{datype}'
        hdr += '\n'*(1023-len(hdr))+'\f'
        with BinaryFile( filename, 'w' ) as datout:
            datout.write(hdr)
            for d in data:
                datout.put( fformat, *d if ndat>1 else d )


def find_format( filename, specform=None ) :
    if specform is None:
        ext = os.path.splitext(filename)[1][1:]
        if ext in format_exts :
            fformat = format_exts[ext]
        else:
            fformat = 'carp'
    else:
        if specform in format_exts.values() :
            fformat = specform
        else:
            print('Illegal input format specified: '+specform)
            sys.exit(1)

    return fformat


#add scalar data field from ASCII file
def add_data( mesh, label, datafile ):
    with open( datafile, 'r' ) as datin:
        lines = datin.readlines()

    data = []
    for l in lines:
        data.append(l.strip().split())

    if len(data) != np.shape(mesh.points)[0] : 
        print( "Data file does not match number of mesh points" )
        sys.exit(1)

    mesh.point_data[label] = np.asarray( data, float )
    
    #mesh = meshio.Mesh( mesh.points, mesh.cells, cell_data=mesh.cell_data, point_data=mesh.point_data )

    print(f'Added data field {label} to mesh')
    return mesh


# simplistic mesh merge, do not look for coincident points
def merge_mesh( m0, m1, reginc):

    if m0 is None : return m1
    if m1 is None : return m0

    npts0  = np.shape(m0.points)[0]
    points = np.concatenate( (m0.points, m1.points) )

    # point data
    for k in m0.point_data:
        if k in m1.points_data :
            m0.point_data[k] = np.concatenate( (m0.point_data[k], m1.point_data.pop(k)) )
        else:
            m0.point_data.pop(k)
            print( f'Ignoring point data: {k}' )
    for k in m1.point_data:
        print( f'Ignoring point data: {k}' )

    # elements
    for cb in m1.cells:
        m0.cells.append( (cb.type, cb.data+npts0) )

    # element data
    for dtype, data in m0.cell_data.items():
        if dtype not in m1.cell_data:
            m0.cell_data.pop(dtype)
            print( f'Ignoring cell data: {dtype}' )
        else:
            if dtype == 'region' : 
                for arr in m1.cell_data['region']: arr+=reginc
            m0.cell_data[dtype] += m1.cell_data.pop(dtype)
    for dtype in m1.cell_data:
        print( f'Ignoring cell data: {dtype}' )

    mesh = meshio.Mesh( points, m0.cells, cell_data=m0.cell_data )
    return mesh


# rename data fields
def map_data( data, mappings ):
    for m in mappings:
        if m[0] in data:
            data[m[1]] = data.pop(m[0])


# make extensions and formats pretty
def prettify_ext():
    e = list(format_exts.keys())
    ml_list=''
    max_line = 6
    while len(e) :
        if len(ml_list) : ml_list += '\t\t'
        np = min(max_line,len(e))
        for i in range(np):
            ml_list += f'{format_exts[e[i]]} (.{e[i]}), '
        e = e[np:]
        if len(e):
            ml_list += '\n'
        else:
            return ml_list


###i###################
# begin main
###i###################
if __name__ == "__main__":

    e = prettify_ext()
    parser = argparse.ArgumentParser(description=f'convert between mesh formats\nknown formats={e}',
                                                            formatter_class=RawTextHelpFormatter)
    parser.add_argument('fromf', nargs='+', help='source mesh(es)' )
    parser.add_argument('tof', nargs='?', default=None, help='target mesh' )
    parser.add_argument('--from-format', default=None, help='specify source mesh format' )
    parser.add_argument('--to-format',   default=None, help='specify target mesh format' )
    parser.add_argument('--pointdata', nargs=2, metavar=('label','file'), default=None, 
                         help='point scalar data' )
#    parser.add_argument('--element-vector-data',default=None, help='file of element vector data' )
#    parser.add_argument('--element-scalar-data',default=None, help='file of element scalar data' )
    parser.add_argument('--precision',default=6, help='number of decimal places for points' )
    parser.add_argument('--binary',action='store_true', help='store in binary format' )
    parser.add_argument('--rm-lower-dim',action='store_true', help='only keep highest dimensional elements' )
    parser.add_argument('--map-cell-data', default=[], action='append',metavar=('old','new'), 
                        nargs=2,help='map old cell data name to new name' )
    parser.add_argument('--map-pt-data',default=[],action='append',metavar=('old','new'),
                       nargs=2,help='map old point data name to new name' )
    parser.add_argument('--del-data', default=[], nargs='+', help='data fields to ignore' )
    parser.add_argument('--sel-regions', default=[], nargs='+', type=int,help='save only these regions' )
    parser.add_argument('--region-inc', type=int, default=100, help='increment element regions by this much for each model')
    parser.add_argument('--keep-unused-pts', action='store_true', help='keep points not connected to any region')

    opts   = parser.parse_args()
    binary = opts.binary

    if opts.tof is None and len(opts.fromf)>1 :
        opts.tof = opts.fromf.pop()

    mesh = None
    for i,fromfs in enumerate(opts.fromf):
        from_format = find_format( fromfs, opts.from_format )
        print(f'Reading {fromfs} in {from_format} format')
        try :
            if from_format in custom_formats:
                Mesh = eval( f"read_{from_format}( fromfs )" )
            else:    
                Mesh = meshio.read( fromfs )
        except Exception as e: 
            print(('\n'+str(e)+'\n'))
            sys.exit(1)

        print( "Finished reading "+fromfs )
        print( Mesh, '\n' )

        mesh = merge_mesh( mesh, Mesh, opts.region_inc*i )

    if opts.tof is None : sys.exit(0)
    to_format = find_format( opts.tof  )

    for d in opts.del_data:
        print(f'Removing data field {d}')
        mesh.point_data.pop( d, None )
        mesh.cell_data.pop( d, None )

    if opts.pointdata: mesh = add_data( mesh, *opts.pointdata )

    map_data( mesh.point_data, opts.map_pt_data )
    map_data( mesh.cell_data,  opts.map_cell_data )

    if opts.sel_regions: mesh = select_regions( mesh, opts.sel_regions )

    if opts.rm_lower_dim : mesh.remove_lower_dimensional_cells()

    if not opts.keep_unused_pts : mesh = rm_unused_pts( mesh )

    print(f"Writing {opts.tof} in {to_format} {'binary' if binary else 'ASCII'} format")
    if to_format  == 'carp' and not binary:
        write_carp( mesh, opts.tof )
    elif to_format in ['carp', 'carpbin']: 
        write_carpb( mesh, opts.tof  )
    else:
        try:
            mesh.write( opts.tof, file_format=to_format, binary=binary )
        except TypeError:
            print('binary not an argument')
            mesh.write( opts.tof, file_format=to_format )

# end of main

