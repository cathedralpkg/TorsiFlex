#!/usr/bin/python3
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
| Program    :  torsiflex          |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#----------------------------------------------------------#
#          Importing basic and personal libraries          #
#----------------------------------------------------------#
import sys
if sys.version_info.major < 3: exit("Python 3 required!!")
#----------------------------------------------------------#
import importlib.util
import traceback
import os
import time
import numpy as np
#----------------------------------------------------------#
import common.files              as ff
import common.internal           as intl
import common.fncs               as fncs
import common.Exceptions         as exc
import common.enantor            as enan
from   common.physcons           import ANGSTROM, AMU
from   common.pgs                import get_pgs
from   common.criteria           import CONNECTSCAL
import common.MolGraph           as mg
from   common.smiles             import smiles_to_cc
from   common.Molecule           import Molecule
from   common.torsions           import torsion_of_cx3
#----------------------------------------------------------#
import modtorsiflex.printing     as pp
import modtorsiflex.tfvars       as tvars
import modtorsiflex.version      as version
import modtorsiflex.tfgau        as itf
import modtorsiflex.tfhelper     as tfh
import modtorsiflex.tfrw         as rw
from   modtorsiflex.options      import get_options_from_prompt
from   modtorsiflex.inpvars      import InpVars
from   modtorsiflex.printing     import sprint
from   modtorsiflex.cc2zmat      import cartesian_to_zmatrix, xyz_to_zmatrix
#----------------------------------------------------------#
from   modtorsiflex.optSEARCH    import search_conformers
from   modtorsiflex.optMSHO      import classify_files
from   modtorsiflex.optHLOPT     import highlevel_reopt
from   modtorsiflex.optREGEN     import regen_from_tmp
from   modtorsiflex.optTORSIONS  import redefine_torsions
from   modtorsiflex.optMSTOR     import gen_mstor
#----------------------------------------------------------#


#==================================================#
def execute_code(function,args,string):
    pp.print_head(string,tvars.NIBS)
    t1 = time.time()
    function(*args)
    t2 = time.time()
    pp.print_elapsed(t2-t1,tvars.NIBS+3)
#--------------------------------------------------#
def print_error():
    message  = str(traceback.format_exc())
    lines    = message.split("\n")
    maxlen   = max([len(line) for line in lines])
    division = "#="+"="*maxlen+"=#"
    print(division)
    for line in lines:
        if line == "": continue
        line += " "*(maxlen-len(line))
        print("# "+line+" #")
    print(division)
#--------------------------------------------------#
def can_smiles_be_used(smiles):
    if smiles is None: return False
    # is rdkit is available??
    answer = True
    pkgs   = ["rdkit","rdkit.Chem"]
    for pkg in pkgs:
        try   : spec = importlib.util.find_spec(pkg)
        except: spec = None
        if spec is None:
           pp.print_not_installed(pkg)
           answer = False
    return answer
#==================================================#

#--------------------------------------------------#
def detect_ch3(cmatrix,symbols,inpvars):
    lCH3 = []
    for idxC,symbol in enumerate(symbols):
        if symbol != "C": continue
        # bonded to?
        bonded_C = [idx2 for idx2,bool2 in enumerate(cmatrix[idxC]) if bool2]
        if len(bonded_C) != 4: continue
        hatoms = sorted([idx for idx in bonded_C if symbols[idx] == "H"])
        if len(hatoms) != 3: continue
        # (A_n)-Y-C-(H_3)
        idxY = [idx for idx in bonded_C if idx not in hatoms][0]
        # bonded atoms to idxY?
        bonded_Y = [idx2 for idx2,bool2 in enumerate(cmatrix[idxY]) if bool2 and idx2 !=idxC]
        # in torsion?
        save = True
        for X,tatoms in inpvars._tatoms.items():
            if sorted(tatoms[1:3]) == sorted([idxY,idxC]):
               save = False
               break
        if save: lCH3.append( [min(bonded_Y),idxY,idxC,min(hatoms)] )
    return lCH3
#--------------------------------------------------#
def detect_nh2(cmatrix,symbols,lzmat,zmatvals,inpvars):
    lNH2 = []
    xcc  = None
    enantio = {}
    for idxN,symbol in enumerate(symbols):
        if symbol != "N": continue
        # check bonded atoms
        bonded_N = [idx2 for idx2,bool2 in enumerate(cmatrix[idxN]) if bool2]
        if len(bonded_N) != 3: continue
        hatoms   = [idx2 for idx2 in bonded_N if symbols[idx2] == "H"]
        ratoms   = [idx2 for idx2 in bonded_N if symbols[idx2] != "H"]
        if len(hatoms) != 2: continue
        # Get Cartesians
        if xcc is None: xcc = intl.zmat2xcc(lzmat,zmatvals)
        # Define dihedral for reference configuration
        idxH1,idxH2 = sorted(hatoms)
        idxR        = ratoms[0]
        HNRH_atoms  = (idxH1,idxN,idxR,idxH2)
        HNRH_value,HNRH_bool = tfh.deal_with_HNRH(HNRH_atoms,xcc)
        lNH2.append( (HNRH_atoms,HNRH_value,HNRH_bool) )
        # in target torsion?
        for X,tatoms in inpvars._tatoms.items():
            atA,atB,atC,atD = tatoms
            # see case and define dihedral for the other H
            if   (atB == idxN and atA == idxH1) or (atC == idxN and atD == idxH1):
               enantio[X] = +HNRH_value
            elif (atB == idxN and atA == idxH2) or (atC == idxN and atD == idxH2):
               enantio[X] = -HNRH_value
            else: continue
    return lNH2

#==================================================#
def torsions_from_zmatfile(zmatfile):
    pp.print_reading_file(zmatfile)
    # read zmatrix
    (lzmat,zmatvals,zmatatoms), symbols, masses = ff.read_zmat(zmatfile)
    # Cartesian coordinates
    xcc = intl.zmat2xcc(lzmat,zmatvals)
    # Molecular Graph
    cmatrix = intl.get_adjmatrix(xcc,symbols,scale=CONNECTSCAL,mode="bool")[0]
    #cmatrix = intl.link_fragments(xcc,cmatrix,nfrags=1)[0]
    # Check torsions
    torsions = []
    for ic,atoms in zmatatoms.items():
        if len(atoms) != 4: continue
        if tfh.is_proper_torsion(atoms,cmatrix): torsions.append(ic)
    # Exclude according to name
    torsions = [torsion for torsion in torsions if not torsion.startswith("itor")]
    torsions = [torsion for torsion in torsions if not torsion.startswith("etor")]
    # Exclude CX3 torsions
    cx3 = [torsion for torsion in torsions if torsion_of_cx3(zmatatoms[torsion],symbols,cmatrix)]
    torsions = [torsion for torsion in torsions if torsion not in cx3]
    # Print info
    pp.print_torsions(torsions,cx3,zmatatoms,symbols)
    # Return
    return torsions,cx3
#--------------------------------------------------#
def readfile_xyz(fname):
    string = ""
    zmat, symbols, masses = ff.read_zmat(fname)
    xcc = intl.zmat2xcc(zmat[0],zmat[1])
    # Print more information
    ndummy = symbols.count("XX")
    natoms = len(symbols)-ndummy
    # without dummies
    symbols_wo , xcc_wo    = fncs.clean_dummies(symbols,xcc=xcc)
    symbols_wo , masses_wo = fncs.clean_dummies(symbols,masses=masses)
    molecule = Molecule()
    molecule.setvar(xcc=xcc_wo,symbols=symbols_wo,masses=masses_wo,ch=0,mtp=1)
    molecule.prepare()
    molecule.setup()
    string += "Molecular formula     : %s\n"%molecule._mform
    string += "Number of atoms       : %i\n"%natoms
    string += "Number of dummy atoms : %i\n"%ndummy
    string += "Vibrational d.o.f.    : %i\n"%molecule._nvdof
    string += "\n"
    # Cartesian Coordinates
    string += "Cartesian coordinates (generated from Z-matrix; in Angstrom):\n"
    sidx = "%%%si"%len(str(len(symbols)))
    at_wo = -1
    dwithout = {}
    for at,symbol in enumerate(symbols):
        xi,yi,zi = [value*ANGSTROM for value in xcc[3*at : 3*at+3]]
        mass     = masses[at]*AMU
        datainline = (sidx%(at+1),symbol,xi,yi,zi,mass)
        if symbol != "XX": at_wo += 1; dwithout[at] = at_wo
        else: dwithout[at] = None
        string += "[%s]  %-2s  %+12.7f  %+12.7f  %+12.7f  (%7.3f amu)\n"%datainline
    string += "\n"
    return xcc,zmat,symbols,masses,(dwithout,molecule,natoms,ndummy),string
#--------------------------------------------------#
def ask_for_enantio(enantio):
    '''Ask for torsional enantiomerism'''
    pp.print_question_enantio(enantio)
    answer  = input(tvars.NIBS*" "+">> ").strip().lower()
    pp.sprint("")
    return answer in ["y","ye","yes"]
#--------------------------------------------------#
def ask_for_ts(ts):
    pp.print_question_ts(ts)
    answer  = input(tvars.NIBS*" "+">> ").strip().lower()
    pp.sprint("")
    return answer == "1"
#--------------------------------------------------#
def ask_for_ch_mtp(ch,mtp):
    # ask for charge
    default_charge  = 0
    pp.print_question_ch(ch)
    ch  = input(tvars.NIBS*" "+">> ").strip().lower().replace(","," ")
    pp.sprint("")
    if ch != "": 
       try   : ch = int(ch)
       except: ch = default_charge
    else     : ch = default_charge
    # ask for multiplicity
    default_multipl = 1
    pp.print_question_mtp(mtp)
    mtp = input(tvars.NIBS*" "+">> ").strip().lower().replace(","," ")
    pp.sprint("")
    if mtp != "": 
       try   : mtp = int(mtp)
       except: mtp = default_multipl
    else     : mtp = default_multipl
    # return data
    return ch,mtp
#--------------------------------------------------#
def create_input():
    # initialize InpVars object
    inpvars = InpVars()
    # write default values
    status  = inpvars.write_default(tvars.IFILE,version.PROGNAME)
    if status == 1: pp.print_inpcreation(tvars.IFILE)
    if status == 0: pp.print_found(tvars.IFILE)
    # Generate Gaussian Templates
    num_files_generated  = itf.generate_templates()
    pp.print_msgGauTempl(num_files_generated,tvars.DIRTEMPL)
    return inpvars
#--------------------------------------------------#
def read_input():
    if not os.path.exists(tvars.IFILE):
       pp.print_notfound(tvars.IFILE)
       raise exc.END
    pp.print_found(tvars.IFILE)
    # initialize InpVars object
    inpvars = InpVars()
    # Read input file
    inpvars.read_input(tvars.IFILE)
    return inpvars
#--------------------------------------------------#
def lonepair(xcc,lNH2):
    lps = {}
    for X,torsion1,torsion2 in lNH2:
        phi1 = fncs.dihedral(*(xcc[3*at:3*at+3] for at in torsion1))
        phi2 = fncs.dihedral(*(xcc[3*at:3*at+3] for at in torsion2))
        phi  = tfh.lonepair_dihedrals_from_phis(phi1,phi2)
        lps[X] = phi
    return lps
#--------------------------------------------------#
def zmat_preparation(inpvars):
    # Read zmat file
    pp.print_found(inpvars._zmatfile)
    try:
        xcc,zmat,symbols,masses,other,string = readfile_xyz(inpvars._zmatfile)
        dwithout,molecule,natoms,ndummy = other
        for line in string.split("\n"): print("       "+line)
    except:
        pp.print_sthwrong(inpvars._zmatfile)
        raise Exception

    # Connectivity matrix
    sprint("   * Adjacency (connection) matrix",tvars.NIBS,1)
    fugdict = inpvars._zmatfile[:-5]+".conn"
    mgraph  = mg.MolGraph(xcc, symbols, cscal=inpvars._cfactor)
    if os.path.exists(fugdict):
       sprint("     - reading connectivity file: %s"%fugdict,tvars.NIBS,1)
       ugdict = rw.read_adjacency(fugdict)
       mgraph.set_from_alist(ugdict)
       cmatrix = mgraph.get_amatrix()
    else:
       sprint("     - connectivity scale factor: %.3f"%inpvars._cfactor,tvars.NIBS)
       sprint("     - writing connectivity file: %s"%fugdict,tvars.NIBS,1)
       mgraph.calculate_connectivity(link_fragments=False)
       ugdict = mgraph._ugdict.copy()
       rw.write_adjacency(fugdict,ugdict)
       cmatrix = np.asarray(mgraph._cmatrix)

    # print bonds
    nrows, ncols = len(cmatrix),len(cmatrix)
    bonds = []
    for idx1 in range(nrows):
        for idx2 in range(idx1+1,ncols):
            if not cmatrix[idx1][idx2]: continue
            bonds.append("%s%i-%s%i"%(symbols[idx1],idx1+1,symbols[idx2],idx2+1))
    sprint("     - number of bonds: %i"%len(bonds),tvars.NIBS)
    ml = max([len(bond) for bond in bonds])
    for idx in range(0,len(bonds),5):
        line = "  ".join(["%%-%is"%ml%bond for bond in bonds[idx:idx+5]])
        sprint(line,tvars.NIBS+7)
    sprint()

    # fragments
    fragments = list(intl.get_fragments_from_adjmatrix(cmatrix))
    sprint("     - number of fragments: %i"%len(fragments),tvars.NIBS)
    for count,fragment in enumerate(fragments):
        frag_mformu = fncs.get_molformula([symbols[idx] for idx in fragment])
        sprint("       (%i) %s"%(count+1,frag_mformu),tvars.NIBS)
    sprint()

    # pairs to be skipped
    if inpvars._skipconn != []:
       sprint("     - Pair(s) of atoms to be skipped in the connectivity test",tvars.NIBS)
       for at1,at2 in inpvars._skipconn:
           sprint("       (%s%i,%s%i)"%(symbols[at1],at1+1,symbols[at2],at2+1),tvars.NIBS)
       sprint()

    lzmat, zmatvals, zmatatoms = zmat

    # Get atoms of target torsions
    sprint("   * Available torsions:",tvars.NIBS)
    for ic,atoms in zmatatoms.items():
        if len(atoms) != 4: continue
        if not intl.isproper(atoms[0],atoms[1],atoms[2],atoms[3],cmatrix=cmatrix): continue
        storsion = "-".join([symbols[atom]+"%i"%(atom+1) for atom in atoms])
        if torsion_of_cx3(atoms,symbols,cmatrix):
           sprint("  %s --> %s  [CX3 rotation]"%(ic,storsion),tvars.NIBS+3)
        else:
           sprint("  %s --> %s"%(ic,storsion),tvars.NIBS+3)
    sprint()

    # Get atoms of target torsions
    sprint("   * selected torsions:",tvars.NIBS)
    for X,ic in zip(inpvars._ttorsions,inpvars._tic):
        if ic not in zmatatoms.keys():
            sprint("   ERROR: Unable to find '%s' in z-matrix..."%ic,tvars.NIBS)
            raise Exception
        inpvars._tatoms[X] = list(zmatatoms[ic])
        storsion = "-".join([symbols[atom]+"%i"%(atom+1) for atom in inpvars._tatoms[X]])
        sprint("  torsion%s (%s) --> %s"%(X,ic,storsion),tvars.NIBS+3)
    sprint()

    # Check lists and dicts related to torsions
    sthwrong = False
    for idx,X in enumerate(inpvars._ttorsions):
        sprint("   * variables for torsion%s (angles in degrees):"%X,tvars.NIBS)
        sprint("tsigma  : %i"%(inpvars._tsigma[idx]),tvars.NIBS+5)
        sprint("maxvalue: %.0f"%(inpvars._tlimit[idx]),tvars.NIBS+5)
        values = "U".join("(%i,%i)"%(phi1,phi2) for phi1,phi2 in inpvars._tdomain[idx])
        sprint("tdomain : %s"%(values),tvars.NIBS+5)
        values = ",".join("%i"%val for val in inpvars._precond[idx])
        sprint("precond : %s"%(values),tvars.NIBS+5)
        sprint()
        if len(inpvars._tdomain[idx]) == 0:
           sprint("ERROR! Something wrong with this torsion!",tvars.NIBS+5)
           raise exc.END

    # Important variable not defined
    status, variable = inpvars.check()
    if status == -1:
       pp.print_readerror(variable)
       end_program = True
       raise Exception

    # Detect CH3 groups
    lCH3 = detect_ch3(cmatrix,symbols,inpvars)
    if len(lCH3) != 0:
        sprint("* CH3 group(s) detected!",tvars.NIBS+3)
        for CH3 in lCH3:
            sprint("dihedral: "+"-".join(["%s%i"%(symbols[idx],idx+1) for idx in CH3]),tvars.NIBS2+2)
        if len(lCH3) == 1: sprint("It will be considered when using --mstor",tvars.NIBS+5)
        else             : sprint("They will be considered when using --mstor",tvars.NIBS+5)
        sprint()

    # Detect NH2 groups
    lNH2 = detect_nh2(cmatrix,symbols,lzmat,zmatvals,inpvars)
    if len(lNH2) != 0:
        sprint("* NH2 group(s) detected!",tvars.NIBS+3)
        for HNRH_atoms,HNRH_value,HNRH_bool in lNH2:
            H1,N,R,H2 = HNRH_atoms
            nh2 = "H%i-N%i-H%i"%(H1+1,N+1,H2+1)
            angle_args = (H1+1,N+1,symbols[R],R+1,H2+1,np.rad2deg(HNRH_value))
            angle = "angle(H%i-N%i-%s%i...H%i) = %.0f degrees"%angle_args
            sprint("%s  --> %s"%(nh2,angle),tvars.NIBS2+2)
        sprint()

    return inpvars,zmat,symbols,masses,cmatrix,lCH3,lNH2,ndummy
#--------------------------------------------------#
def check_optmode(inpvars,ndummy):
    if ndummy != 0 and inpvars._optmode != 0:
       sprint("WARNING! Dummy atoms detected! Keyword 'optmode' will be set to 0!!",tvars.NIBS,1)
       inpvars._optmode = 0
    return inpvars
#==================================================#



#==================================================#
def prepare_system():
    # Read input file
    inpvars = read_input()
    inpvars.prepare_variables()

    # zmat file exists?
    if not os.path.exists(inpvars._zmatfile):
       pp.print_filenotfound(inpvars._zmatfile,tvars.NIBS,1)
       raise exc.END

    # Read zmat, check variables, generate templates and update current variables
    if len(inpvars._ttorsions) == 0:
       sprint("No torsions selected in input file!",tvars.NIBS)
       if len(options) == 0:
          sprint("No options were given! Use --help for more information",tvars.NIBS,1)
       raise exc.END

    inpvars,zmat,symbols,masses,cmatrix,lCH3,lNH2,ndummy = zmat_preparation(inpvars)

    # Enantiomers? Correlate numbering
    dcorr = None
    if inpvars._enantio:
       sprint("--> enantio keyword activated!",tvars.NIBS2)
       sprint("    testing on reference...",tvars.NIBS2)
       lzmat, zmatvals, zmatatoms = zmat
       # reference geom
       xcc0 = intl.zmat2xcc(lzmat,zmatvals)
       vec0 = tfh.xcc2vec(xcc0,inpvars)
       sprint("    reference : %s"%str(vec0),tvars.NIBS2)
       try:
          # (1) First, try checking Cs symmetry!
          # 1.1 - create geom with all torsions at 180
          zmatvals2 = zmatvals.copy()
          for ic in inpvars._tic: zmatvals2[ic] = 180.0
          xcc0 = intl.zmat2xcc(lzmat,zmatvals2)
          # 1.2 check is point group is Cs
          pg, rotsigma = get_pgs(symbols,masses,xcc0)
          if pg == "Cs":
             dcorr = "Cs"
             sprint("    enantiomer: phi --> -phi",tvars.NIBS2,1)
          # (2) Correlate atoms with enantiomer
          else:
             # enantiomer
             xcc_enantio = enan.generate_enantio(xcc0)
             dcorr       = enan.correlate_enantio(xcc0,xcc_enantio,symbols)
             xcc_enantio = enan.reorder_enantio(xcc_enantio,dcorr)
             vec_enantio = tfh.xcc2vec(xcc_enantio,inpvars)
             sprint("    enantiomer: %s"%str(vec_enantio),tvars.NIBS2,1)
       except:
          dcorr = None
          sprint("    something went wrong... enantio deactivated!",tvars.NIBS2,1)
          inpvars._enantio = False

    # dealing with TS?
    if not inpvars._ts:
       inpvars._ifqrangeLL = []
       inpvars._ifqrangeHL = []

    return inpvars,zmat,symbols,masses,cmatrix,lCH3,lNH2,ndummy,dcorr
#==================================================#



#==================================================#
def option_smiles(smiles,zmatfile):
    fugdict = zmatfile[:-4]+".conn"
    # Print info
    pp.print_smiles_geomfile(smiles,zmatfile)
    # Deal with smiles
    if can_smiles_be_used(smiles):
       if os.path.exists(zmatfile):
          pp.print_file_exists(zmatfile,"Warning")
          pp.print_no_smiles()
       else:
          pp.print_converting_smiles()
          # smiles --> zmat
          xcc, symbols = smiles_to_cc(smiles)
          cartesian_to_zmatrix(xcc,symbols,fugdict,zmatfile)
          pp.print_smiles_done()
    else:
       pp.print_no_smiles()
#--------------------------------------------------#
def option_cartesian(fxyz,fzmat):
    # Assert existence of files
    if not os.path.exists(fxyz) : pp.print_filenotfound(fxyz,tvars.NIBS,1) ; return
    if     os.path.exists(fzmat): pp.print_file_exists(fzmat,"Warning"   ) ; return
    # Convert
    pp.print_cartesian(fxyz,fzmat)
    xyz_to_zmatrix(fxyz,fzmat)
#--------------------------------------------------#
def option_input(zmatfile):

    # Read file (if exists)
    if not zmatfile.endswith(".zmat"):
       pp.print_geomfile_ext()
       raise exc.END

    #-----------------------------#
    # Create TorsiFlex input file #
    #-----------------------------#
    create = True
    if os.path.exists(tvars.IFILE):
       pp.print_file_exists(tvars.IFILE,"WARNING")
       answer = input(tvars.NIBS*" "+"re-write file [Y/n]? ").strip().lower()
       print("")
       if answer in ["y","yes",""]: create = True
       else                       : create = False

    if create:
       inpvars = InpVars()
       #- - - - - - - - - - - - - - - - - - -#
       # Ask for some info to complete input #
       #- - - - - - - - - - - - - - - - - - -#
       pp.print_question_defaults()
       # (a) min or ts
       inpvars._ts = ask_for_ts(inpvars._ts)
       # (b) Charge and multiplicity
       inpvars._charge, inpvars._multipl = ask_for_ch_mtp(inpvars._charge,inpvars._multipl)
       # (c) torsional enantiomerism
       inpvars._enantio = ask_for_enantio(inpvars._enantio)
       # print output
       pp.print_question_answers(inpvars._ts,inpvars._charge,inpvars._multipl,inpvars._enantio)
       #- - - - - - - - - - - - - - - - - - -#
       # Read torsions from z-matrix         #
       #- - - - - - - - - - - - - - - - - - -#
       if os.path.exists(zmatfile):
          torsions,cx3 = torsions_from_zmatfile(zmatfile)
       else:
          torsions = None
       #- - - - - - - - - - - - - - - - - - -#
       # Write inputfile using previous info #
       #- - - - - - - - - - - - - - - - - - -#
       pp.print_inpcreation(tvars.IFILE)
       status  = inpvars.write_default(tvars.IFILE,version.PROGNAME,torsions,True)

    # Create Gaussian templates
    nf_gen = itf.generate_templates()
    if nf_gen != 0: pp.print_msgGauTempl(tvars.GAUTEMPL)
#==================================================#



#==================================================#
def main():

    # Presentation and arguments
    options  = get_options_from_prompt()

    #============================#
    # OPTION == --help           #
    #============================#
    if "-h" in options or "--help" in options:
        pp.print_welcome(version.PROGNAME)
        sprint(pp.HSTRING)
        raise exc.END
    #----------------------------#

    #============================#
    # OPTION == --version        #
    #============================#
    if "-v" in options or "--version" in options:
        sprint(version.PROGNAME.split(".py")[0]+" "+version.PROGVER)
        raise exc.END
    #----------------------------#

    pp.print_welcome(version.PROGNAME)
    pp.print_user_info()

    #============================#
    # OPTION == --smiles         #
    #============================#
    if "--smiles" in options:
       execute_code(option_smiles,options["--smiles"],"Cartesian coordinates from SMILES")
       raise exc.END
    #============================#

    #============================#
    # OPTION == --cartesian      #
    #============================#
    if "--cartesian" in options:
       execute_code(option_cartesian,options["--cartesian"],"Cartesian coordinates to Z-matrix")
       raise exc.END
    #============================#

    #============================#
    # OPTION == --input          #
    #============================#
    if "--input" in options:
      execute_code(option_input,(options["--input"],),"Input creation")
      raise exc.END
    #----------------------------#

    #============================#
    # NO INPUT FILE              #
    #============================#
    if not os.path.exists(tvars.IFILE):
      pp.print_notfound(tvars.IFILE)
      raise exc.END
    #----------------------------#

    inpvars,zmat,symbols,masses,cmatrix,lCH3,lNH2,ndummy,dcorr = prepare_system()

    #============================#
    # NO OPTIONS                 #
    #============================#
    if len(options) == 0:
       pp.print_no_options()
       raise exc.END
    #----------------------------#

    #============================#
    # OPTION == --prec or --stoc #
    #============================#
    if "--prec" in options or "--stoc" in options:
        inpvars = check_optmode(inpvars,ndummy)

        # Preconditioned search? Stochastic search?
        inpvars._prec = options.get("--prec",False)
        inpvars._stoc = options.get("--stoc",[]   )

        args = (inpvars,zmat,symbols,cmatrix,dcorr,lNH2)
        execute_code(search_conformers,args,"Low-Level search")
    #----------------------------#

    #============================#
    # OPTION == --msho ll        #
    #============================#
    if options.get("--msho",None) in ["ll","all"]:
       argsLL = (inpvars,cmatrix,dcorr,"LL")
       execute_code(classify_files,argsLL,"Low-Level MSHO partition functions")
    #----------------------------#

    #============================#
    # OPTION == --mstor ll       #
    #============================#
    if options.get("--mstor",None) in ["ll","all"]:
       args = (inpvars,dcorr,lCH3,"LL")
       execute_code(gen_mstor,args,"Generating MsTor input (LL)")
    #----------------------------#

    #============================#
    # OPTION == --hlopt          #
    #============================#
    if "--hlopt" in options:
       inpvars = check_optmode(inpvars,ndummy)
       mode_hlopt,lowerconformer,upperconformer = options["--hlopt"]
       if   mode_hlopt == "0":
          args = (inpvars,cmatrix,lowerconformer,upperconformer,True,lNH2)
          execute_code(highlevel_reopt,args,"High-Level calculations")
       elif mode_hlopt == "1":
          args = (inpvars,cmatrix,lowerconformer,upperconformer,False,lNH2)
          execute_code(highlevel_reopt,args,"Generating High-Level input files")
          raise exc.END
    #----------------------------#

    #============================#
    # OPTION == --msho hl        #
    #============================#
    if options.get("--msho",None) in ["hl","all"]:
       argsHL = (inpvars,cmatrix,dcorr,"HL")
       execute_code(classify_files,argsHL,"High-Level MSHO partition functions")
    #----------------------------#

    #============================#
    # OPTION == --mstor hl       #
    #============================#
    if options.get("--mstor",None) in ["hl","all"]:
       args = (inpvars,dcorr,lCH3,"HL")
       execute_code(gen_mstor,args,"Generating MsTor input (HL)")
    #----------------------------#

    #============================#
    # OPTION == --regen          #
    #============================#
    if "--regen" in options:
       args = (inpvars,dcorr)
       execute_code(regen_from_tmp,args,"Regenerating domains")
       raise exc.END
    #----------------------------#

    #============================#
    # OPTION == --torsions       #
    #============================#
    if "--torsions" in options:
       args = (inpvars,dcorr,zmat,cmatrix,symbols)
       execute_code(redefine_torsions,args,"Re-defining torsions")
       raise exc.END
    #----------------------------#

#==================================================#



#==================================================#
if __name__ == '__main__':
   try           : main()
   except exc.END: pass
   except exc.ErrorTorsion1 as error:
       EXITMESS = "\n    First torsion must be named torsion1 instead of torsion%i!\n"%error._var
       print(EXITMESS)
   except exc.ErrorTorsionN as error:
       EXITMESS = "\n    There are %i target torsions, but last torsion is torsion%i!\n"%error._var
       print(EXITMESS)
   except exc.ErrorTorsionRepe as error:
       EXITMESS = "\n    Target torsion (torsion%s) is duplicated in input file!\n"%error._var
       print(EXITMESS)
       EXITMESS = "\n  ERROR while executing %s!\n"%version.PROGNAMEnopy
       print(EXITMESS)
       print_error()
   except exc.WrongDimension as exception:
       sprint("ERROR! Wrong dimension in file '%s'!"%inpvars._pcfile,tvars.NIBS2)
       sprint("First line of file contains %i torsion(s)!"%exception._ntor,tvars.NIBS2+7)
       sprint("--> %s "%exception._fline,tvars.NIBS2+7)
       sprint("Line which causes this error:",tvars.NIBS2+7)
       sprint("--> %s "%exception._line,tvars.NIBS2+7)
       sprint()
       raise exc.END
   except exc.StocSearchWrongDim as exception:
       sprint("ERROR! Wrong dimension in stochastic vector!",tvars.NIBS2)
       sprint("Input    dimension: %i (%s)"%(exception._dim[0],exception._vec),tvars.NIBS2+7)
       sprint("Expected dimension: %i"%exception._dim[1],tvars.NIBS2+7)
       sprint()
       raise exc.END
   except exc.UnableGenRandAng as exception:
        sprint("ERROR! Unable to generate valid random angle!",tvars.NIBS2)
        if len(exception._domain) == 0:
           sprint("Empty domain for one of the torsions!",tvars.NIBS2+7)
        else:
           the_domain = "U".join(["(%.0f,%.0f)"%(p1,p2) for p1,p2 in exception._domain])
           sprint("Domain for torsion: %s"%the_domain,tvars.NIBS2+7)
        sprint()
        raise exc.END
   except exc.ErrorHConstraint:
        sprint("ERROR! Problem(s) when checking hard constraint(s)!\n",tvars.NIBS2)
        raise exc.END
   except exc.ErrorSConstraint:
        sprint("ERROR! Problem(s) when checking soft constraint(s)!\n",tvars.NIBS2)
        raise exc.END
#==================================================#
    

