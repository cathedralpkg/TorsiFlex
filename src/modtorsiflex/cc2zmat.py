'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: TorsiFlex
Version     : 2022.1
License     : MIT/x11

Copyright (c) 2022, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
---------------------------

*----------------------------------*
| Module     :  modtorsiflex       |
| Sub-module :  cc2zmat            |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#==============================================================#
import numpy as np
import os
import sys
#--------------------------------------------------------------#
sys.setrecursionlimit(10000) # change maximum recursion depth  #
#--------------------------------------------------------------#
import common.files    as ff
import common.fncs     as fncs
import common.internal as intl
import common.MolGraph as mg
import common.physcons as pc
from   common.criteria import CONNECTSCAL,EPS_LINEAR
#--------------------------------------------------------------#
from   modtorsiflex.tfrw import read_adjacency, write_adjacency
#==============================================================#




#==============================================================#
def key_for_dihedral(nodeD,nodeC,nodeB,nodeA,ugdict,symbols,incycle,brpairs):
    ''' ugdict must include dummy atoms'''
    bond = set([nodeB,nodeC])
    if not intl.isproper(nodeA,nodeB,nodeC,nodeD,ugdict): key = "itor"
    elif   nodeB not in incycle                    : key = "ptor"
    elif   nodeC not in incycle                    : key = "ptor"
    elif   bond      in brpairs                    : key = "ptor"
    else                                           : key = "rtor"
    # Exclude N=N, C=N and C=C proper torsions
    if key == "ptor":
       nodeB_C = symbols[nodeB].upper() == "C" and len(ugdict[nodeB]) < 4
       nodeB_N = symbols[nodeB].upper() == "N" and len(ugdict[nodeB]) < 3
       nodeC_C = symbols[nodeC].upper() == "C" and len(ugdict[nodeC]) < 4
       nodeC_N = symbols[nodeC].upper() == "N" and len(ugdict[nodeC]) < 3
       if (nodeB_C or nodeB_N) and (nodeC_C or nodeC_N): key = "etor"
    return key
#--------------------------------------------------------------#
def string_zmat(CACHE):
    ''' returns the string with the z-matrix '''

    # Data of interest from cache
    xcc     = CACHE["XCC"]
    symbols = CACHE["SYMBS"]
    zmatrix = CACHE["ZMAT" ]
    ugdict  = CACHE["ALIST"]
    incycle = CACHE["INCYCLE"]
    brpairs = CACHE["BRPAIRS"]

    # some formatting
    nn = fncs.num_of_digits(len(zmatrix.keys()))
    iformat = "%"+"0%i"%nn+"i"

    # Create a copy of ugdict including dummy atoms
    cugdict = ugdict.copy()
    for node,conns in zmatrix.items():
        if type(node) != int: cugdict[node] = [conns[0]]

    # New labelling
    old2new = {node:idx+1 for idx,node in enumerate(zmatrix.keys())}

    # Generate string
    string,keys,proper,idx = "",{},[],1
    for nodeD,(nodeC,nodeB,nodeA) in zmatrix.items():
        # current node
        string += " %-2s"%(symbols[nodeD])
        # distance D-C
        xD = xcc.get(nodeD,None)
        xC = xcc.get(nodeC,None)
        xB = xcc.get(nodeB,None)
        xA = xcc.get(nodeA,None)
        if nodeC is not None:
           key = "dist"+iformat%idx
           string += "  %2i  %7s"%(old2new[nodeC],key)
           keys[key] = fncs.distance(xD,xC) * pc.ANGSTROM
        # bond angle D-C-B
        if nodeB is not None:
           key = "angl"+iformat%idx
           string += "  %2i  %7s"%(old2new[nodeB],key)
           keys[key] = np.rad2deg( fncs.angle(xD,xC,xB))
        # dihedral D-C-B-A
        if nodeA is not None:
           args = (nodeD,nodeC,nodeB,nodeA,cugdict,symbols,incycle,brpairs)
           key  = key_for_dihedral(*args)+iformat%idx
           # save proper torsion
           if "ptor" in key: proper.append(key)
           # add to string
           string += "  %2i  %7s"%(old2new[nodeA],key)
           keys[key] = np.rad2deg( fncs.dihedral(xD,xC,xB,xA))
        # end of line
        idx += 1
        string += "\n"
    string += "\n"
    # Add values for each key
    for key,value in keys.items():
        string += "%-6s   %10.5f\n"%(key,value)
    string += "\n"
    return string
#==============================================================#

#==============================================================#
def same_cycle(node1,node2,cycle,ugdict):
    visited = set([node2])
    tovisit = ugdict[node2].intersection(cycle).difference(visited.union([node1]))
    did_i_visit_node_1 = False
    while len(tovisit) != 0:
        # visited node
        node_i = tovisit.pop()
        visited.add(node_i)
        # update nodes to visit
        neighs  = ugdict[node_i].intersection(cycle).difference(visited)
        tovisit = tovisit.union(neighs)
        # node1 in neighs???
        if node1 in neighs:
           did_i_visit_node_1 = True
           break
    return did_i_visit_node_1
#--------------------------------------------------------------#
def separate_cycles(cycle,ugdict):
    '''
    For a complex cycle, check if any pair of bonded atoms(*)
    belong to the same ring;

    Returns a list of pairs where one node belongs to one ring
    and the other to another ring

    This is useful to notice that Ph-Ph is actually made by two
    separated rings with one proper torsion of interest

    (*) both with more than 2 connections each
    '''
    broken_rpairs = []
    nodes = [node for node in cycle if len(ugdict[node].intersection(cycle)) > 2]
    for idxi,ni in enumerate(nodes):
        for nj in nodes[idxi+1:]:
            if nj in ugdict[ni]:
                # can we go ni --> nj in a cycle?
                did_i_visit_node_1 = same_cycle(ni,nj,cycle,ugdict)
                if not did_i_visit_node_1:
                   broken_rpairs.append( set([ni,nj]) )
    return broken_rpairs
#==============================================================#

#==============================================================#
def require_dummy(node_i,node_j,node_k,xcc):
    '''
    for i-j-k, calculates bond angle and
    decides if dummy atoms is required

    xcc is a dictionary xcc[node_j] --> xcc for node_j
    '''
    if None in (node_i,node_j,node_k): return False
    # Cartesian Coordinates
    if type(xcc) == dict:
       x_i = list(xcc[node_i])
       x_j = list(xcc[node_j])
       x_k = list(xcc[node_k])
    else:
       x_i = list(xcc[3*node_i:3*node_i+3])
       x_j = list(xcc[3*node_j:3*node_j+3])
       x_k = list(xcc[3*node_k:3*node_k+3])
    # Calculate angle g-h-i
    angle = np.rad2deg(fncs.angle(x_i,x_j,x_k))
    # Check if dummy is required
    return abs(180-angle) < EPS_LINEAR or angle < EPS_LINEAR
#--------------------------------------------------------------#
def add_dummy(ni,nj,nk,nl,CACHE):
    '''
    new node: nl              ni
                               \
    dummy is required           nj-nk-[nl]
                                   |
                                   X
    '''
    nX = "dummy%i"%nk
    xk = CACHE["XCC"][nk]
    xj = CACHE["XCC"][nj]
    xi = CACHE["XCC"].get(ni,None)
    # coordinates for dummy atom
    if xi is None: xX = intl.zmat_thirdatom(xj,xk   ,1.0,np.pi/2    )
    else         : xX = intl.zmat_nextpoint(xi,xj,xk,1.0,np.pi/2,0.0)
    # update xcc
    CACHE["XCC"][nX] = np.array(xX)
    # update symbols
    CACHE["SYMBS"][nX] = "XX"
    # update zmatrix
    CACHE["ZMAT"][nX] = (nk,nj,ni)
    # update CONNS for nk
    if None not in (nj,nk,nX): CACHE["CONNS"][nk] = (nj,nX)
    # return
    return CACHE,nX
#--------------------------------------------------------------#
def update_zmatrix(paths,CACHE,iscycle=False):

    for m0,path in paths:
        # first node of fragment
        n0 = path[0]
        # Nodes of previous fragment for connections
        if m0 is None:
           m1,m2 = None,None
        elif m0 in CACHE["CONNS"]:
           m1,m2 = CACHE["CONNS"][m0]
           if require_dummy(m1,m0,n0,CACHE["XCC"]):
              m1,m2 = m2,m1
              assert not require_dummy(m1,m0,n0,CACHE["XCC"])
        else:
           m1,m2 = CACHE["ENDS"][m0]
        # Update zmatrix
        ni,nj,nk = m2,m1,m0
        for idx_l,nl in enumerate(path):
            # (1) check if dummy is required
            if require_dummy(nj,nk,nl,CACHE["XCC"]):
               CACHE,nX = add_dummy(ni,nj,nk,nl,CACHE)
               # update reference nodes
               nj,ni = nX,nj
            # (2) Update CONNS
            elif idx_l != 0 and None not in (nj,nk,nl):
               CACHE["CONNS"][nk] = (nj,nl)
            # (3) LAST NODE IN CHAIN
            if idx_l+1 == len(path) and not iscycle and None not in (nl,nk,nj):
                  CACHE["ENDS"][nl] = (nk,nj)
            # (4) LAST NODE IN CYCLE
            if idx_l+1 == len(path) and iscycle:
               # connection for last node
               n_0 = CACHE["ALIST"][nl].intersection(CACHE["ZMAT"].keys()).difference(path[1:])
               if len(n_0) == 1: n_0 = list(n_0)[0]
               else            : n_0 = list(n_0.difference(path[0:1]))[0]
               if None not in (nk,nl,n_0): CACHE["CONNS"][nl] = (nk,n_0)
               # connection for first node if FIRST CYCLE
               if m0 is None:
                  CACHE["CONNS"][path[0]] = (path[1],nl)
               # rename nodes for the definition of nl
               else:
                  nk = n_0
                  nj,ni = CACHE["CONNS"][nk]
                  if require_dummy(nl,nk,nj,CACHE["XCC"]): nj,ni = ni,nj
                  assert not require_dummy(nl,nk,nj,CACHE["XCC"])
            # (4) add to z-matrix
            CACHE["ZMAT"][nl] = (nk,nj,ni)
            # (5) update nodes
            ni,nj,nk = nj,nk,nl
    # return dictionary
    return CACHE
#==============================================================#




#==============================================================#
def fragment_it_belongs(node,noncycles,cycles):
    fragment = [(fragment,False) for fragment in noncycles if node in fragment] + \
               [(fragment,True ) for fragment in    cycles if node in fragment] 

    if   len(fragment) == 0: return None,None
    elif len(fragment) == 1: return fragment[0]
    raise Exception
#--------------------------------------------------------------#
def prepare_system(xcc,symbols,fugdict,write_conn=False):
    '''
    Prepares system from XYZ file
    '''

    #------------------#
    # (1) Create graph #
    #------------------#

    mgraph = mg.MolGraph(xcc, symbols, cscal=CONNECTSCAL, epslin=EPS_LINEAR)

    if os.path.exists(fugdict):
       ugdict = read_adjacency(fugdict)
       mgraph.set_from_alist(ugdict)
    else:
       mgraph.calculate_connectivity()
       ugdict = mgraph._ugdict.copy()
       if write_conn: write_adjacency(fugdict,ugdict)

    #-------------------#
    # (2) Detect cycles #
    #-------------------#

    mgraph.detect_cycles()
    num_in_cycle = len(mgraph._incycle)

    #------------------------#
    # (3) Fragment the graph #
    #------------------------#

    cycles,noncycles,only1,Hs = mgraph.split_graph()

    #------------------#
    # (4) Create CACHE #
    #------------------#

    CACHE = {}
    CACHE["XCC"      ] = {n:np.array(xcc[3*n:3*n+3]) for n   in range(len(symbols))}
    CACHE["SYMBS"    ] = {n:str(s)                   for n,s in enumerate(symbols) }
    CACHE["ALIST"    ] = ugdict
    CACHE["CYCLES"   ] = list(cycles)
    CACHE["NONCYCLES"] = list(noncycles)
    CACHE["ONLY1"    ] = list(only1)
    CACHE["HS"       ] = list(Hs)
    CACHE["INCYCLE"  ] = [node for cycle in cycles for node in cycle]
    CACHE["BRPAIRS"  ] = [pair for cycle in cycles for pair in separate_cycles(cycle,ugdict)]
    CACHE["ZMAT"     ] = {}
    CACHE["CONNS"    ] = {}
    CACHE["ENDS"     ] = {}

    return mgraph,CACHE
#--------------------------------------------------------------#
def cartesian_to_zmatrix(xcc,symbols,fugdict,ofile=None):

    # prepare system
    mgraph,CACHE = prepare_system(xcc,symbols,fugdict)

    # (1) FIRST PATH
    cycles    = CACHE["CYCLES"]
    noncycles = CACHE["NONCYCLES"]
    # (a) there are no cycles
    if   len(cycles) == 0                   : usecycle = False
    # (b) there are only cycles
    elif len(noncycles) == 0                : usecycle = True
    # (c) the biggest cycle contains more atoms than biggest non-cyclic
    elif len(cycles[0]) >= len(noncycles[0]): usecycle = True
    # (d) the biggest cycle contains less atoms than biggest non-cyclic
    else                                    : usecycle = False


    # Deal with first fragment
    if usecycle:
       fragm0  = CACHE["CYCLES"].pop(0)
       paths   = mgraph.paths_from_cycles(fragm0,m0=None)
       FIRST_PATH = []
    else:
       fragm0  = CACHE["NONCYCLES"].pop(0)
       paths   = mgraph.paths_from_chains(fragm0,CACHE["SYMBS"],m0=None)
       FIRST_PATH = paths[0][1]
      #FIRST_PATH = list(CACHE["ZMAT"].keys())

    # Update z-matrix
    CACHE = update_zmatrix(paths,CACHE,usecycle)

    # (2) CONNECT FRAGMENTS TO PREVIOUS PATHS
    tovisit = list(fragm0)
    while len(tovisit) != 0:
        # Get node
        m0 = tovisit.pop(0)
        # Get neighbors of node (that weren't visited yet)
        neighbors = set(CACHE["ALIST"][m0]).difference(CACHE["ZMAT"].keys())
        # Add new fragments
        for n0 in neighbors:
            # determine fragment:
            fragment,iscycle = fragment_it_belongs(n0,CACHE["NONCYCLES"],CACHE["CYCLES"])
            # Fragment is a single atom that will be added later
            if fragment is None:
               tovisit += [n0]
               continue
            elif iscycle: paths = mgraph.paths_from_cycles(fragment,m0=m0)
            else        : paths = mgraph.paths_from_chains(fragment,CACHE["SYMBS"],m0=m0)
            CACHE = update_zmatrix(paths,CACHE,iscycle)
            tovisit += fragment

    # (3) SINGLE ATOMS (DIFFERENT FROM H) CONNECTED TO CYCLES
    for n0 in CACHE["ONLY1"]:
        m0 = list(CACHE["ALIST"][n0])[0]
        paths = [(m0,[n0])]
        CACHE = update_zmatrix(paths,CACHE,False)

    # (4) TERMINAL H ATOM AT BEGINNING/END OF FIRST PATH
    Hs = list(CACHE["HS"])
    for atH in Hs:
        # H-A-B-C-...
        atA = list(CACHE["ALIST"][atH])[0]
        if atA not in FIRST_PATH: continue
        # assert atA not in CONNS
        if atA in CACHE["CONNS"]: continue
        # (a) at the beggining of FIRST_PATH
        if atA == FIRST_PATH[0]:
           try   : atB = FIRST_PATH[1]
           except: atB = None
           try   : atC = FIRST_PATH[2]
           except: atC = None
           FIRST_PATH = [atH]+FIRST_PATH
        # (b) at the end of FIRST_PATH
        else:
           try   : atB = FIRST_PATH[-2]
           except: atB = None
           try   : atC = FIRST_PATH[-3]
           except: atC = None
           FIRST_PATH = FIRST_PATH+[atH]
        # atB is None?
        if atB is None:
           CACHE["ZMAT" ][atH] = (atA,atB,atC)
        # Dummy required?
        elif require_dummy(atB,atA,atH,CACHE["XCC"]):
           CACHE,atX = add_dummy(atC,atB,atA,atH,CACHE)
           # update connection & ZMATRIX
           CACHE["CONNS"][atA] = (atB,atX)
           CACHE["ZMAT" ][atH] = (atA,atX,atB)
        else:
           CACHE["CONNS"][atA] = (atB,atH)
           CACHE["ZMAT" ][atH] = (atA,atB,atC)
    # update list of H atoms
    Hs = [H for H in Hs if H not in CACHE["ZMAT"]]

    # (5) TERMINAL H ATOMS INVOLVED AT THE END OF A TORSION (C-B-A-H)
    for atH in Hs:
        # H-A-B-C
        atA = list(CACHE["ALIST"][atH])[0]
        if atA not in CACHE["ENDS"]: continue
        if atA     in CACHE["CONNS"]: continue
        # at the END of any path
        atB,atC = CACHE["ENDS"][atA]
        # Dummy required?
        if require_dummy(atB,atA,atH,CACHE["XCC"]):
           CACHE,atX = add_dummy(atC,atB,atA,atH,CACHE)
           # update connection & ZMATRIX
           CACHE["CONNS"][atA] = (atB,atX)
           CACHE["ZMAT" ][atH] = (atA,atX,atB)
        else:
           CACHE["CONNS"][atA] = (atB,atH)
           CACHE["ZMAT" ][atH] = (atA,atB,atC)
    Hs = [H for H in Hs if H not in CACHE["ZMAT"]]

    # (6) H ATOMS NOT INVOLVED IN TORSIONS
    for n0 in Hs:
        m0      = list(CACHE["ALIST"][n0])[0]
        paths   = mgraph.paths_from_chains([n0],CACHE["SYMBS"],m0)
        CACHE = update_zmatrix(paths,CACHE,False)

    # Write Z-matrix file
    if ofile is not None:
       with open(ofile,'w') as asdf: asdf.write(string_zmat(CACHE))
    else:
       return CACHE
#--------------------------------------------------------------#
def xyz_to_zmatrix(fxyz,ofile):

    if fxyz.endswith(".xyz"): fugdict = fxyz[:-4]+".conn"
    else                    : fugdict = fxyz     +".conn"

    # Read xyz and prepare cache variable Create CACHE #
    xcc, symbols, masses = ff.read_xyz(fxyz)
    if len(symbols) < 3: raise Exception

    # Write z-matrix file
    cartesian_to_zmatrix(xcc,symbols,fugdict,ofile)
#==============================================================#



