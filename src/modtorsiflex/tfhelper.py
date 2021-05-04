'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: TorsiFlex
Version     : 2021.2
License     : MIT/x11

Copyright (c) 2021, David Ferro Costas (david.ferro@usc.es) and
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
| Sub-module :  tfhelper           |
| Last Update:  2020/12/21 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different
functions used in TorsiFlex
'''

#==========================================#
import os
import numpy as np
import random
#------------------------------------------#
import common.fncs       as fncs
import common.internal   as intl
import common.Exceptions as exc
from   common.physcons   import ANGSTROM
from   common.files      import mkdir_recursive
#------------------------------------------#
import modtorsiflex.printing    as pp
from   modtorsiflex.tpespoint   import TorPESpoint
import modtorsiflex.enantor     as enan
#==========================================#


#==================================================#
def create_dir(directory):
    if os.path.exists(directory): return
    try: mkdir_recursive(directory)
    except:
        pp.print_createdir(directory,tvars.NIBS2)
        raise Exception
#==================================================#

#==================================================#
def deal_with_HNRH(HNRH_atoms,xcc):
    '''
    HNRH_atoms --> idxH1, idxN, idxR, idxH2
    idxH1: index for first H atom
    idxN : index for N atom
    idxR : index for R group
    idxH2: index for second H atom
    '''
    HNRH_value = fncs.dihedral(*(xcc[3*at:3*at+3] for at in HNRH_atoms))
    HNRH_bool  = np.rad2deg(HNRH_value)%360 < 180
    return HNRH_value,HNRH_bool
#--------------------------------------------------#
def lonepair_dihedrals_from_phis(phi1,phi2):
    # to degrees
    phi1,phi2 = np.rad2deg(phi1)%360,np.rad2deg(phi2)%360
    # sort
    phi1,phi2 = min(phi1,phi2), max(phi1,phi2)
    # dihedral for lone pair
    if (phi2-phi1) < 180: dihedral = ((360+phi1+phi2)/2)%360
    else                : dihedral = ((phi2+phi1)/2)%360
    return np.deg2rad(dihedral)
#--------------------------------------------------#
def lonepair_dihedral_from_xs(xN,xH1,xH2,xB,xA):
    '''
    A-B-NH2
    return dihedral A-B-N-lp, lp being the lone pair
    '''
    # dihedral A-B-N-H1
    phi1 = fncs.dihedral(xA,xB,xN,xH1)
    # dihedral A-B-N-H2
    phi2 = fncs.dihedral(xA,xB,xN,xH2)
    # dihedral A-B-N-lp
    return lonepair_dihedrals_from_phis(phi1,phi2)
#==================================================#

#==================================================#
def adjmatrix_from_zmatrix(lzmat,zmatvals,cfactor):
    # convert to cartesian coordinates
    xcc = intl.zmat2xcc(lzmat,zmatvals)
    # symbols
    symbols = [pack[0] for pack in lzmat]
    # Get connection matrix
    cmatrix = intl.get_adjmatrix(xcc,symbols,scale=cfactor,mode="bool")[0]
    return cmatrix
#--------------------------------------------------#
def random_name(nn=5):
    characters = "abcdefghijklmnopqrstuvwxyz"
    numbers    = "0123456789"
    allch      = characters+characters.upper()+numbers
    return "".join([random.choice(allch) for ii in range(nn)])
#--------------------------------------------------#
def random_angle(tdomain=[(0,360)]):
    for ii in range(10**6):
        # random angle as integer
        angle = int(round(random.random()*360))
        # assert in range
        if fncs.float_in_domain(angle,tdomain): return angle
    # unable to generate angle
    exception = exc.UnableGenRandAng(Exception)
    exception._domain = tdomain
    raise exception
#--------------------------------------------------#
def detect_nh2(symbols,cmatrix):
    NH2s = []
    for idx,symbol in enumerate(symbols):
        if symbol != "N": continue
        bonded = [index for index,boolen in enumerate(cmatrix[idx]) if boolen and symbols[index]=="H"]
        if len(bonded) != 2: continue
        NH2s.append( [idx]+sorted(bonded) )
    return NH2s
#--------------------------------------------------#
def bool2yesno(logical_variable):
    if logical_variable: return "yes"
    else               : return "no"
#--------------------------------------------------#
def str2interval(string,bint=False):
    intervals = [interval.replace("(","").replace(")","").split(",") \
                 for interval in string.split("U")]
    intervals = [(float(phi1),float(phi2)) for (phi1,phi2) in intervals]
    # exclude intervals where phi1>phi2
    intervals = [(phi1,phi2) for (phi1,phi2) in intervals if phi1<=phi2]
    # deal with negative values (for angles)
    final_intervals = []
    for phi1,phi2 in intervals:
        if   phi1 < 0.0 and phi2 < 0.0:
             final_intervals.append( (phi1%360,phi2%360) )
        elif phi1 < 0.0 and phi2 == 0.0:
             final_intervals.append( (    0.0,0.0) )
             final_intervals.append( (phi1%360,360.0) )
        elif phi1 < 0.0 and phi2 >  0.0:
             final_intervals.append( ( 0.0, phi2) )
             final_intervals.append( (phi1%360,360.0) )
        else:
             final_intervals.append( (phi1,phi2) )
    # data as integers?
    if bint: final_intervals = [(int(round(p1)),int(round(p2))) for p1,p2 in final_intervals] 
    return final_intervals
#--------------------------------------------------#
def name2vec(name):
    return tuple([int(value)%360 for value in name.split("_")])
#--------------------------------------------------#
def vec2name(vec):
    return "_".join(["%003.0f"%(value%360) for value in vec])
#--------------------------------------------------#
def vecinlimits(vec,tlimit):
    for vx,lx in zip(vec,tlimit):
        if vx > lx: return False
    return True
#--------------------------------------------------#
def zmat2name(zmatvals,torsions):
    vec = zmat2vec(zmatvals,torsions)
    return vec2name(vec)
#--------------------------------------------------#
def zmat2vec(zmatvals,torsions,negative=False):
    if negative:
       return [(-zmatvals[torsion])%360 for torsion in torsions]
    else:
       return [zmatvals[torsion]%360 for torsion in torsions]
#--------------------------------------------------#
def xcc2vec(xcc,inpvars):
    phis = []
    for X in inpvars._ttorsions:
        xs  = [xcc[3*at:3*at+3] for at in inpvars._tatoms[X]]
        phis.append( np.rad2deg(fncs.dihedral(*xs)) )
    return TorPESpoint(phis)
#--------------------------------------------------#
def folder_points(folder):
    if not os.path.exists(folder): return []
    # only those log files with freq calculations
    vecs = [fname.split(".")[1] for fname in os.listdir(folder) \
            if fname.endswith(".log")]
    # The vectors
    return [TorPESpoint(point) for point in vecs]
#--------------------------------------------------#
def random_point(tdomain,tlimits):
    vec = [random_angle(domain) for domain in tdomain]
    return TorPESpoint(vec,tlimits)
#--------------------------------------------------#
def yield_angles(inpvars):
    # angles generated according pre-defined cases
    if inpvars._prec:
       for vec in inpvars.yield_precond(): yield TorPESpoint(vec)
    # angles generated randomly
    else:
       for cycle in range(inpvars._ncycles):
           yield random_point(inpvars._tdomain,inpvars._tlimit)
#--------------------------------------------------#
def enantiovec(lzmat,zmatvals,dcorr,inpvars):
    if dcorr == "Cs":
       vec = zmat2vec(zmatvals,inpvars._tic,negative=True)
       return TorPESpoint(vec)
    else:
       xcc0        = intl.zmat2xcc(lzmat,zmatvals)
       xcc_enantio = enan.generate_enantio(xcc0)
       xcc_enantio = enan.reorder_enantio(xcc_enantio,dcorr)
       vec_enantio = xcc2vec(xcc_enantio,inpvars)
    return vec_enantio
#==================================================#


#==================================================#
def in_explored_domain(vec,ddomains,inpvars):
    # Get the points of the explored domain
    points  = [fpoint for domain in ddomains.keys() \
                      for ipoint,fpoint,statusFRQ in ddomains[domain] \
                      if fpoint is not None]
    points += [ipoint for domain in ddomains.keys() \
                      for ipoint,fpoint,statusFRQ in ddomains[domain] \
                      if ipoint is not None]
    if len(points) == 0: return False, None
    # (a) already explored?
    spoints = [str(point) for point in points]
    if str(vec) in spoints: return True, vec
    # (b) Get closest point
    closest = vec.closest(points)
    inside  = vec.is_inside(closest,inpvars._dist1d)
    if inside: return True , closest
    return False, closest
#--------------------------------------------------#
def test_similarity_redundancy(point,inpvars,ddomains):
    '''
    returns:
       -1 --> point in domain
        0 --> belongs to domain of other point
        1 --> guess not in visited domain
    '''
    # Check visited domain
    bool_indomain, closest = in_explored_domain(point,ddomains,inpvars)
    # (a) Point already evaluated
    if bool_indomain and str(closest) == str(point): return -1, closest
    # (b) Point in domain
    elif bool_indomain: return 0, closest
    # (c) guess not in visited domain
    else:
       points  = folder_points(inpvars._dirll)
       if len(points) == 0: return 1, None
       closest = point.closest(points)
       return 1, closest
#--------------------------------------------------#
def test_connectivity(lzmat,zmatvals,cfactor,cmatrix):
    # calculate connectivity of current geometry
    cmatrix_current = adjmatrix_from_zmatrix(lzmat,zmatvals,cfactor)
    # return if test is passed
    return cmatrix == cmatrix_current
#--------------------------------------------------#
def test_hessian():
    pass
#--------------------------------------------------#
def test_hsconstraints(lzmat,zmatvals,constr,which="hard"):
    if constr is None or len(constr) == 0: return True
    # the xcc
    xcc = None
    # check constraints
    for ic,icdomain in constr:
        if ic in zmatvals:
           icvalue = zmatvals[ic]
           if icvalue < 0.0: icvalue = icvalue%360
        else:
           # Get atoms
           icatoms  = [int(at)-1 for at in ic.split("-")]
           nicatoms = len(icatoms)
           # Get xcc for each atom
           if xcc is None: xcc = intl.zmat2xcc(lzmat,zmatvals)
           xyzs     = [fncs.xyz(xcc,at) for at in icatoms]
           # calculate value in xcc
           if   nicatoms == 2: icvalue =  fncs.distance(*xyzs)*ANGSTROM
           elif nicatoms == 3: icvalue = np.rad2deg(fncs.angle(*xyzs)   )%360
           elif nicatoms == 4: icvalue = np.rad2deg(fncs.dihedral(*xyzs))%360
           else: raise Exception
        # in domain?
        boolean = fncs.float_in_domain(icvalue,icdomain)
        if which=="hard" and boolean is False: return False
        if which=="soft" and boolean is True : return True
    # (hard) all of them are fulfilled
    if which=="hard": return True
    # (soft) all of them failed
    if which=="soft": return False
#--------------------------------------------------#
##  def test_hconstraints(xcc,hconstr):
##      if hconstr is None or len(hconstr) == 0: return True
##      # check hard constraints
##      for ic,icdomain in hconstr:
##          # Get atoms
##          icatoms  = [int(at)-1 for at in ic.split("-")]
##          nicatoms = len(icatoms)
##          xyzs     = [fncs.xyz(xcc,at) for at in icatoms]
##          # calculate value in xcc
##          if   nicatoms == 2: icvalue =  fncs.distance(*xyzs)*ANGSTROM
##          elif nicatoms == 3: icvalue = np.rad2deg(fncs.angle(*xyzs)   )%360
##          elif nicatoms == 4: icvalue = np.rad2deg(fncs.dihedral(*xyzs))%360
##          else: raise Exception
##          # in domain?
##          boolean = fncs.float_in_domain(icvalue,icdomain)
##          if boolean is False: return False
##      # all of them are fulfilled
##      return True
##  #--------------------------------------------------#
##  def test_sconstraints(xcc,sconstr):
##      if sconstr is None or len(sconstr) == 0: return True
##      # check soft constraints
##      for ic,icdomain in sconstr:
##          # Get atoms
##          icatoms  = [int(at)-1 for at in ic.split("-")]
##          nicatoms = len(icatoms)
##          xyzs     = (fncs.xyz(xcc,at) for at in icatoms)
##          # calculate value in xcc
##          if   nicatoms == 2: icvalue = fncs.distance(*xyzs)*ANGSTROM
##          elif nicatoms == 3: icvalue = np.rad2deg(fncs.angle(*xyzs))
##          elif nicatoms == 4: icvalue = np.rad2deg(fncs.dihedral(*xyzs))%360
##          else: raise Exception
##          # in domain?
##          boolean = fncs.float_in_domain(icvalue,icdomain)
##          if boolean is True: return True
##      # all of them failed
##      return False
#--------------------------------------------------#
def precheck_geom(lzmat,zmatvals,cmatrix,inpvars):
    bools = [None,None,None,None]
    # (a) connectivity test
    bools[0] = test_connectivity(lzmat,zmatvals,inpvars._cfactor,cmatrix)
    # (b) in domain test
    vec = TorPESpoint([zmatvals[ic] for ic in inpvars._tic])
    bools[1] = vec.is_in_domain(inpvars._tdomain)
    # (c) Hard constraints
    try   : bools[2] = test_hsconstraints(lzmat,zmatvals,inpvars._hconstr,"hard")
    except: raise exc.ErrorHConstraint
    # (d) Soft constraints
    try   : bools[3] = test_hsconstraints(lzmat,zmatvals,inpvars._sconstr,"soft")
    except: raise exc.ErrorSConstraint
    # return all comparisons
    return bools
#==================================================#

