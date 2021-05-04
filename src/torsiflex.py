#!/usr/bin/python3
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
| Program    :  torsiflex          |
| Last Update:  2020/12/21 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#----------------------------------------------------------#
#          Importing basic and personal libraries          #
#----------------------------------------------------------#
import sys
if sys.version_info.major < 3: exit("Python 3 required!!")
#----------------------------------------------------------#
import traceback
import os
import time
import numpy as np
#----------------------------------------------------------#
import common.internal           as intl
import common.fncs               as fncs
import common.Exceptions         as exc
from   common.pgs                import get_pgs
#----------------------------------------------------------#
import modtorsiflex.printing     as pp
import modtorsiflex.tfvars       as tvars
import modtorsiflex.tfgau        as itf
import modtorsiflex.tfhelper     as tfh
import modtorsiflex.tfrw         as rw
import modtorsiflex.enantor      as enan
from   modtorsiflex.inpvars      import InpVars
from   modtorsiflex.printing     import sprint
#----------------------------------------------------------#
from   modtorsiflex.optSEARCH    import search_conformers
from   modtorsiflex.optMSHO      import classify_files
from   modtorsiflex.optHLOPT     import highlevel_reopt
from   modtorsiflex.optREGEN     import regen_from_tmp
from   modtorsiflex.optMSTOR     import gen_mstor
#----------------------------------------------------------#

OPTIONS__  = "help,version,inp,"
OPTIONS__ += "prec,stoc,hlopt,"
OPTIONS__ += "msho,mstor,"
OPTIONS__ += "regen,"
OPTIONS_   = "h,v"


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
#==================================================#


#==================================================#
def get_options_from_prompt():
    # Get user arguments (options)
    args = sys.argv[1:]
    # The valid options
    valid_args = ["--"+arg for arg in OPTIONS__.split(",")]+\
                 [ "-"+arg for arg in OPTIONS_.split(",")]
    # bools
    argbools = {option:False for option in OPTIONS__.split(",")}
    # Are they known?
    current   = None
    prec1     = None
    prec2     = None
    mode_hlopt = "0"
    mode_msho  = "all"
    mode_mstor = "all"
    for arg in args:
        if arg in valid_args: current = arg
        # arguments for option --prec
        elif current == "--prec" and not arg.startswith("-"):
             try:
                if   prec1 is None: prec1 = int(arg)
                elif prec2 is None: prec2 = int(arg)
             except:
                sprint("arguments for --prec must be integers!!")
                raise exc.END
             continue
        # argument for option --hlopt
        elif current == "--hlopt" and not arg.startswith("-"):
             if arg == "nocalc": mode_hlopt = "1"; continue
             else: sprint("argument for --hlopt can only be 'nocalc'!!"); raise exc.END
        # argument for option --msho
        elif current == "--msho" and not arg.startswith("-"):
             if   arg == "ll": mode_msho = "ll"; continue
             elif arg == "hl": mode_msho = "hl"; continue
             else: sprint("argument for --msho can only be 'll' ot 'hl'!!"); raise exc.END
        # argument for option --mstor
        elif current == "--mstor" and not arg.startswith("-"):
             if   arg == "ll": mode_mstor = "ll"; continue
             elif arg == "hl": mode_mstor = "hl"; continue
             else: sprint("argument for --mstor can only be 'll' ot 'hl'!!"); raise exc.END

        if not arg.startswith(""): continue
        elif arg not in valid_args:
             pp.print_unknown_option(arg,tvars.PROGNAME)
             raise exc.END
        argbools[arg[2:]] = True
    # value for prec
    if prec1 is None: prec1 = 1
    if prec2 is None: prec2 = 1
    if argbools["prec"]: argbools["prec"] = (max(prec1,prec2,1),max(min(prec1,prec2),1))
    # value for hlopt, msho and mstor
    if argbools["hlopt"] is not False: argbools["hlopt"] = mode_hlopt
    if argbools["msho" ] is not False: argbools["msho" ] = mode_msho
    if argbools["mstor"] is not False: argbools["mstor"] = mode_mstor
    # check if user asks for help
    if "-h" in args or "--help" in args:
        argbools["help"] = True
        pp.print_welcome(tvars.PROGNAME)
        sprint(pp.HSTRING)
        raise exc.END
    # check if user asks for version
    if "-v" in args or "--version" in args:
        argbools["version"] = True
        sprint(tvars.PROGNAME.split(".py")[0]+" "+tvars.PROGVER)
        raise exc.END
    pp.print_welcome(tvars.PROGNAME)
    pp.print_user_info()
    return argbools
#--------------------------------------------------#
def deal_with_input(case="create"):
    # initialize InpVars object
    inpvars = InpVars()
    # act according case
    if case == "create":
       status  = inpvars.write_default(tvars.IFILE,tvars.PROGNAME)
       if status == 1: pp.print_inpcreation(tvars.IFILE)
       if status == 0: pp.print_found(tvars.IFILE)
       # Generate Gaussian Templates
       num_files_generated  = itf.generate_templates()
       pp.print_msgGauTempl(num_files_generated,tvars.DIRTEMPL)
    elif case == "read":
       if os.path.exists(tvars.IFILE): pp.print_found(tvars.IFILE)
       else: pp.print_notfound(tvars.IFILE); raise exc.END
       inpvars.read_input(tvars.IFILE)
    # return InpVars
    return inpvars
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
        xcc,zmat,symbols,masses,other,string = rw.readfile_xyz(inpvars._zmatfile)
        dwithout,molecule,natoms,ndummy = other
        for line in string.split("\n"): print("       "+line)
    except:
        pp.print_sthwrong(inpvars._zmatfile)
        raise Exception

    # Get connectivity matrix
    sprint("   * calculating adjacency (connection) matrix...",tvars.NIBS)
    sprint("     - connectivity scale factor: %.3f"%inpvars._cfactor,tvars.NIBS)
    cmatrix,dmatrix,nbonds = intl.get_adjmatrix(xcc,symbols,scale=inpvars._cfactor,mode="bool")
    fragments = list(intl.get_fragments_from_adjmatrix(cmatrix))
    sprint("     - number of bonds: %i"%nbonds,tvars.NIBS)
    # print bonds
    nrows, ncols = dmatrix.shape
    bonds = []
    for idx1 in range(nrows):
        for idx2 in range(idx1+1,ncols):
            if not cmatrix[idx1][idx2]: continue
            bonds.append("%s%i-%s%i"%(symbols[idx1],idx1+1,symbols[idx2],idx2+1))
    ml = max([len(bond) for bond in bonds])
    for idx in range(0,len(bonds),5):
        line = "  ".join(["%%-%is"%ml%bond for bond in bonds[idx:idx+5]])
        sprint(line,tvars.NIBS+7)
    sprint()
    # fragments
    sprint("     - number of fragments: %i"%len(fragments),tvars.NIBS)
    for count,fragment in enumerate(fragments):
        frag_mformu = fncs.get_molformula([symbols[idx] for idx in fragment])
        sprint("       (%i) %s"%(count+1,frag_mformu),tvars.NIBS)
    sprint()

    # Get atoms of target torsions
    sprint("   * identifying target torsions in z-matrix...",tvars.NIBS)
    lzmat, zmatvals, zmatatoms = zmat
    for X,ic in zip(inpvars._ttorsions,inpvars._tic):
        if ic not in zmatatoms.keys():
            sprint("   ERROR: Unable to find '%s' in z-matrix..."%ic,tvars.NIBS)
            raise Exception
        inpvars._tatoms[X] = list(zmatatoms[ic])
        storsion = "-".join([symbols[atom]+"%i"%(atom+1) for atom in inpvars._tatoms[X]])
        sprint("  torsion%s (%s): %s"%(X,ic,storsion),tvars.NIBS+3)
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
#==================================================#








#==================================================#
def main():

    # Presentation and arguments
    argsbools = get_options_from_prompt()

    num_options = len([argsbools[option] for option in OPTIONS__.split(",") \
                                    if argsbools[option] is not False])

    # No options selected & no input- print info
    if num_options == 0 and not os.path.exists(tvars.IFILE):
       pp.print_notfound(tvars.IFILE)
       sprint("No options were given! Use --help for more information",tvars.NIBS,1)
       raise exc.END

    # Input & Templates creation
    if argsbools["inp"]: inpvars = deal_with_input(case="create"); raise exc.END
    else               : inpvars = deal_with_input(case="read")

    # No options selected & no zmatrix - print info
    if num_options == 0 and not os.path.exists(inpvars._zmatfile):
       pp.print_notfound(inpvars._zmatfile)
       sprint("No options were given! Use --help for more information",tvars.NIBS,1)
       raise exc.END

    # zmat file exists?
    if not os.path.exists(inpvars._zmatfile):
       pp.print_filenotfound(inpvars._zmatfile,tvars.NIBS,1)
       raise exc.END

    # Read zmat, check variables, generate templates and update current variables
    inpvars.prepare_variables()
    if len(inpvars._ttorsions) == 0:
       sprint("No torsions selected in input file!",tvars.NIBS)
       if num_options == 0:
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

    # Check opt mode
    if ndummy != 0 and inpvars._optmode != 0:
       sprint("WARNING! Dummy atoms detected! Keyword 'optmode' will be set to 0!!",tvars.NIBS,1)
       inpvars._optmode = 0

    #=========================#
    # Act according option(s) #
    #=========================#
    # No options selected - print info
    if num_options == 0:
       sprint("",tvars.NIBS)
       sprint("No options were given! Use --help for more information",tvars.NIBS,1)
       raise exc.END

    if argsbools["hlopt"] == "1":
       args = (inpvars,cmatrix,False)
       execute_code(highlevel_reopt,args,"Generating High-Level input files")
       raise exc.END

    # --mstor: Generate mstor files
    if argsbools["mstor"] in ["ll","all"]:
       args = (inpvars,"ll",dcorr,lCH3)
       execute_code(gen_mstor,args,"Generating MsTor input (LL)")
    if argsbools["mstor"] in ["hl","all"]:
       args = (inpvars,"hl",dcorr,lCH3)
       execute_code(gen_mstor,args,"Generating MsTor input (HL)")
    if argsbools["mstor"] is not False: raise exc.END

    # --regen: Regerate data 
    if argsbools["regen"]:
       args = (inpvars,cmatrix)
       execute_code(regen_from_tmp,args,"Regenerating domains")
       raise exc.END

    # --msho: Just checking
    if argsbools["msho"] in ["ll","all"]:
       argsLL = (inpvars,cmatrix,"LL",dcorr)
       execute_code(classify_files,argsLL,"Low-Level MSHO partition functions")
    if argsbools["msho"] in ["hl","all"]:
       argsHL = (inpvars,cmatrix,"HL",dcorr)
       execute_code(classify_files,argsHL,"High-Level MSHO partition functions")
    if argsbools["msho"] is not False: raise exc.END

    # LL search (+opt+freq)
    if argsbools["prec" ] is not False or argsbools["stoc"]:
        # Preconditioned search?
        inpvars._prec = argsbools["prec"]

        try: 
            args = (inpvars,zmat,symbols,cmatrix,dcorr,lNH2)
            execute_code(search_conformers,args,"Low-Level search")
        except exc.WrongDimension as exception:
            sprint("ERROR! Wrong dimension in file '%s'!"%inpvars._pcfile,tvars.NIBS2)
            sprint("First line of file contains %i torsion(s)!"%exception._ntor,tvars.NIBS2+7)
            sprint("--> %s "%exception._fline,tvars.NIBS2+7)
            sprint("Line which causes this error:",tvars.NIBS2+7)
            sprint("--> %s "%exception._line,tvars.NIBS2+7)
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

    # HL optimization+freq
    if argsbools["hlopt"] == "0":
        args = (inpvars,cmatrix,True)
        execute_code(highlevel_reopt,args,"High-Level calculations")
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
   except        :
       EXITMESS = "\n  ERROR while executing %s!\n"%tvars.PROGNAMEnopy
       print(EXITMESS)
       print_error()
#==================================================#
    

