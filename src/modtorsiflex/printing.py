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
| Sub-module :  printing           |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different
prints used in TorsiFlex
'''

#==================================================#
import os
import sys
import time
import getpass
import numpy as np
#--------------------------------------------------#
from   common.fncs          import time2human, fill_string, afreq2cm
from   common.physcons      import KCALMOL, KB
#--------------------------------------------------#
from   modtorsiflex.shelp   import HSTRING
from   modtorsiflex.version import PROGNAME, PROGNAMEnopy
from   modtorsiflex.version import VERSION, DATE
from   modtorsiflex.tfvars  import IFILE, DIRTEMPL, GAUTEMPL
from   modtorsiflex.tfvars  import NIBS, NIBS2, IBS, IBS2
#==================================================#

DESCRDOMAINS  = {\
"conf" :"Set of conformers",\
"enan" :"Set of enantiomers",\
"repe" :"Set of optimized points that match with a stored conformer",\
"excl" :"Set of optimized points that were excluded",\
"wimag":"Set of optimized points with a wrong number of imaginary frequencies",\
"fail" :"Set of guess points whose optimization failed",\
}



#==================================================#
def sprint(string="",nibs=0,nbl=0):
    print(" "*nibs+string+"\n"*nbl)
#==================================================#

#==================================================#
def print_option_error(option):
    info = []
    if   option == "--smiles":
         info.append("Option --smiles must be followed by two arguments;")
         info.append("* 1st argument: the SMILES code, between quotation marks;")
         info.append("* 2nd argument: the output zmat file.")
    elif option == "--cartesian":
         info.append("Option --cartesian must be followed by two arguments;")
         info.append("* 1st argument: name of the .xyz file (with Cartesian coordinates);")
         info.append("* 2nd argument: name of the output .zmat file.")
    elif option == "--input":
         info.append("Option --input admits only one argument (the name of the Z-matrix file).")
    elif option == "--prec":
         info.append("Arguments for --prec, if any, must be two integers")
    elif option == "--hlopt":
         info.append("Arguments for --hlopt, if any, must be one or two integers")
         info.append("and/or the string 'nocalc'")
    elif option == "--msho":
         info.append("Argument for --msho can only be 'll' or 'hl'")
    elif option == "--mstor":
         info.append("Argument for --mstor can only be 'll' or 'hl'")

    # Print info for this error
    print_welcome(PROGNAME)
    if len(info) == 0: return
    sprint("ERROR! %s"%info[0],NIBS)
    for line in info[1:]: sprint("       "+line,NIBS)
    sprint("")
#--------------------------------------------------#
def print_user_info():
    cdate = time.strftime("%Y-%m-%d")
    ctime = time.strftime("%H:%M:%S")
    user  = getpass.getuser()
    host  = os.uname()[1]
    pwd   = os.getcwd()
    vinfo = "%i.%i.%i"%(sys.version_info[0:3])
    sprint("--------------------------------",NIBS)
    sprint(" Program information",NIBS)
    sprint("  -name    : %s"%PROGNAMEnopy,NIBS)
    sprint("  -version : %s"%VERSION,NIBS)
    sprint("  -date    : %s"%DATE,NIBS)
    sprint("--------------------------------",NIBS)
    sprint(" Python    : %s"%vinfo,NIBS)
    sprint(" User      : %s"%user,NIBS)
    sprint(" Host      : %s"%host,NIBS)
    sprint(" Directory : %s"%pwd,NIBS)
    sprint(" Date      : %s (YY-MM-DD)"%cdate,NIBS)
    sprint(" Time      : %s (HH:MM:SS)"%ctime,NIBS)
    sprint("--------------------------------",NIBS)
    sprint()
#--------------------------------------------------#
def print_no_options():
    sprint("No options were given! Use --help for more information",NIBS,1)
#--------------------------------------------------#
def print_smiles_geomfile(smiles,geomfile):
    sprint("SMILES code   : %s"%smiles,NIBS)
    sprint("Z-matrix file : %s"%geomfile,NIBS)
    sprint("")
#--------------------------------------------------#
def print_geomfile_ext():
    sprint("ERROR! The geometry file MUST end with .zmat extension!",NIBS)
    sprint("")
#--------------------------------------------------#
def print_not_installed(pkg):
    sprint("Package not installed: %s"%pkg,NIBS)
#--------------------------------------------------#
def print_no_smiles():
    sprint("The use of SMILES codes is disabled!",NIBS,1)
#--------------------------------------------------#
def print_converting_smiles():
    sprint("Generating Z-matrix from SMILES code...",NIBS)
    sprint("")
#--------------------------------------------------#
def print_smiles_done():
    sprint("Done! SMILES code converted to Z-matrix!",NIBS)
    sprint("")
#--------------------------------------------------#
def print_no_smiles():
    sprint("SMILES code will not be converted to Z-matrix...",NIBS)
    sprint("")
#--------------------------------------------------#
def print_smileswarning(zmatfile):
    sprint("WARNING! %s file already exists!"%zmatfile,NIBS)
    sprint("         .xyz --> .zmat conversion will be omitted!",NIBS)
    sprint("")
#--------------------------------------------------#
def print_reading_file(filename):
    sprint("Reading file (%s) ..."%filename,NIBS)
    sprint("")
#--------------------------------------------------#
def print_file_exists(filename,message="ERROR"):
    sprint("%s! File '%s' already exists..."%(message,filename),NIBS)
    sprint("")
#--------------------------------------------------#
def print_question_defaults():
    sprint("System information (default):",NIBS)
    sprint("")
    sprint("System corresponds to a : 0 (i.e. minimum)",NIBS+4)
    sprint("Total charge            : 0",NIBS+4)
    sprint("Spin multiplicity       : 1",NIBS+4)
    sprint("Torsional enantiomerism : n (i.e. no)",NIBS+4)
    sprint("")
    sprint("If ok, press on Enter in the corresponding question!",NIBS)
    sprint("")
#--------------------------------------------------#
def print_question_answers(ts,ch,mtp,enantio):
    sprint("System information (by user):",NIBS)
    sprint("")
    if ts:
        sprint("System corresponds to a : 1 (i.e. transition state)",NIBS+4)
    else:
        sprint("System corresponds to a : 0 (i.e. minimum)",NIBS+4)
    sprint("Total charge            : %i"%ch,NIBS+4)
    sprint("Spin multiplicity       : %i"%mtp,NIBS+4)
    if enantio:
       sprint("Torsional enantiomerism : y (i.e. yes)",NIBS+4)
    else:
       sprint("Torsional enantiomerism : n (i.e. no)",NIBS+4)
    sprint("")
#--------------------------------------------------#
def print_question_ts(ts):
    if ts: current = "1"
    else : current = "0"
    sprint("Is your system a minimum (0) or a transition state (1)? Current = %s"%current,NIBS)
#--------------------------------------------------#
def print_question_ch(ch):
    sprint("Total charge? Current = %i"%ch,NIBS)
#--------------------------------------------------#
def print_question_mtp(mtp):
    sprint("Spin multiplicity? Current = %i"%mtp,NIBS)
#--------------------------------------------------#
def print_question_enantio(enantio):
    if enantio: current = "y"
    else      : current = "n"
    sprint("Does this molecule present torsional enantiomerism [y/N]? Current = %s"%current,NIBS)
#--------------------------------------------------#
def print_cartesian(xyzfile,zmatfile):
    sprint("Converting Cartesian coordinates to Z-matriz...",NIBS)
    sprint("")
    sprint("  - Cartesian coordinates: %s"%xyzfile,NIBS)
    sprint("  - Z-matrix  coordinates: %s"%zmatfile,NIBS)
    sprint("")
#--------------------------------------------------#
def print_torsions(torsions,cx3,zmatatoms,symbols):
    sprint("Identified torsions:",NIBS,1)
    nn = max([len(torsion) for torsion in torsions+cx3])
    for torsion in torsions+cx3:
        atoms = zmatatoms[torsion]
        chain = "-".join( ["%s%i"%(symbols[atom],atom+1) for atom in atoms] )
        line  = "  * %%-%is <--> %%s"%nn%(torsion,chain)
        if torsion in cx3: line += "  [CX3 rotation] <-- excluded"
        sprint(line,NIBS)
    sprint("")
#--------------------------------------------------#
def print_creating_input_torsi(fname):
    sprint("Creating input file: %s"%fname,NIBS)
    sprint("")
#--------------------------------------------------#
def print_domains(ddomains,DOMAINS):
    ntot = 0
    for domain in DOMAINS:
        npoints = len(ddomains[domain])
        ntot   += npoints
        sprint("- points in domain %-5s: %i"%(domain.upper(),npoints),NIBS2+2)
    sprint("- TOTAL NUMBER OF POINTS: %i"%(ntot),NIBS2+2)
    sprint()
    for domain in DOMAINS:
        description = DESCRDOMAINS[domain]
        sprint("- %-5s: %s"%(domain.upper(),description),NIBS2+2)
    sprint()
#--------------------------------------------------#
def print_filenotfound(filename,nibs=0,bl=0):
    string = "File NOT FOUND: %s"%filename
    sprint(string,nibs)
    for i in range(bl): sprint()
#--------------------------------------------------#
def print_dirnotfound(dirname,nibs=0,bl=0):
    string = "Folder NOT FOUND: %s"%dirname
    sprint(string,nibs)
    for i in range(bl): sprint()
#--------------------------------------------------#
def print_head(string,nibs=0):
    sprint()
    sprint("="*len(string),nibs)
    sprint(        string ,nibs)
    sprint("="*len(string),nibs)
    sprint()
#--------------------------------------------------#
def print_elapsed(delta_time,nibs):
    t,units = time2human(delta_time)
    sprint("--------------------------------",NIBS)
    sprint("Elapsed time: %.2f %s"%(t,units),nibs)
    sprint("--------------------------------",NIBS)
    sprint()
#--------------------------------------------------#
def print_tail():
    cdate = time.strftime("%Y-%m-%d")
    ctime = time.strftime("%H:%M:%S")
    sprint("--------------------------------",NIBS)
    sprint(" Date      : %s (YY-MM-DD)"%cdate,NIBS)
    sprint(" Time      : %s (HH:MM:SS)"%ctime,NIBS)
    sprint("--------------------------------",NIBS)
    sprint()
#--------------------------------------------------#
def print_createdir(folder,nibs):
    sprint("- the user should create folder: %s"%folder,nibs)
#--------------------------------------------------#
def print_template_notfound():
    sprint("ERROR! Template for Electronic Structure Calculations",NIBS)
    sprint("       NOT FOUND (%s)"%GAUTEMPL,NIBS)
    sprint()
#--------------------------------------------------#
def print_error_template(key):
    sprint("ERROR! Template associated to keyword '%s'"%key,NIBS)
    sprint("       was NOT FOUND inside '%s'"%GAUTEMPL,NIBS)
    sprint()
#--------------------------------------------------#
def print_template(template,key,scase,nibs):
    sprint("Template for %s (%s):"%(scase,key),nibs)
    sprint("---------------------------",nibs)
    for line in template: sprint(line[:-1],nibs)
    sprint("---------------------------",nibs)
    sprint()
#--------------------------------------------------#
def print_infoLLsearch(inpvars,keyoptLL,keyfrqLL,loptLL,lfrqLL):
    # (a) tolerance
    sprint("Similarity test tolerance: %5.1f degrees"%inpvars._dist1d,NIBS2)
    sprint("Redundancy test tolerance: %5.1f degrees"%inpvars._epsdeg,NIBS2)
    sprint()
    # (b) LL opt template
    print_template(loptLL,keyoptLL,"LL optimization",NIBS2)
    # (c) LL freq template
    print_template(loptLL,keyfrqLL,"LL frequency calculation",NIBS2)
#--------------------------------------------------#
def print_tests(tests,which=0):
    sprint("Tests:",NIBS2)

    if which == 0:
       sprint("  - Connectivity test on Guess structure: %s"%(tests[0][0] == 1),NIBS2)
       sprint("  - Similarity   test on Guess structure: %s"%(tests[0][1] == 1),NIBS2)
       sprint("  - H-Constraint test on Guess structure: %s"%(tests[0][2] == 1),NIBS2)
       sprint("  - S-Constraint test on Guess structure: %s"%(tests[0][3] == 1),NIBS2)
       sprint()

    sprint("  - Connectivity test on Opt structure  : %s"%(tests[1][0] == 1),NIBS2)
    sprint("  - Reduncancy   test on Opt structure  : %s"%(tests[1][1] == 1),NIBS2)
    sprint("  - H-Constraint test on Opt structure  : %s"%(tests[1][2] == 1),NIBS2)
    sprint("  - S-Constraint test on Opt structure  : %s"%(tests[1][3] == 1),NIBS2)

    sprint("")
#--------------------------------------------------#
def print_ifreqconstr(interval):
    if len(interval) == 0: return
    sprint("  - Imaginary frequency restriction DEFINED!",NIBS2)
    for if1,if2 in interval:
        sprint("    * From %7.2fi to %7.2fi"%(if1,if2),NIBS2)
    sprint("")
#--------------------------------------------------#
def print_excluded_ifreq(ifreq):
    TEXT_EXCL = "[excluded!] imaginary freq., %.2fi cm^-1, outside defined interval"
    sprint(TEXT_EXCL%ifreq,NIBS2+4)
#--------------------------------------------------#
def print_accepted_ifreq(ifreq):
    TEXT_ACCE = "[accepted!] imaginary freq., %.2fi cm^-1, inside defined interval"
    sprint(TEXT_ACCE%ifreq,NIBS2+4)
#--------------------------------------------------#
def print_numcalcs(count,nopt,nfrq,nsp):
    sprint("About number of calculations:",NIBS2)
    sprint("--> Number of guess points                 : %i"%count,NIBS2+4)
    sprint("--> Number of optimizations (nopt)         : %i"%nopt,NIBS2+4)
    sprint("--> Number of frequency calculations (nfrq): %i"%nfrq,NIBS2+4)
    sprint("--> Number of stationary points found (nsp): %i"%nsp,NIBS2+4)
    if nfrq != 0: sprint("--> nopt/nfrq = %.4f"%(float(nopt)/float(nfrq)),NIBS2+4)
    if nsp  != 0: sprint("--> nopt/nsp  = %.4f"%(float(nopt)/float(nsp )),NIBS2+4)
    if nsp  != 0: sprint("--> nfrq/nsp  = %.4f"%(float(nfrq)/float(nsp )),NIBS2+4)
    sprint()
#--------------------------------------------------#
def print_exploredvol(explored,ii,nn):
    sprint("Estimating explored volume by Monte-Carlo...",NIBS2)
    sprint("total number of points used for Monte-Carlo (nn): %i"%nn,NIBS2+4)
    sprint("number of points inside explored space      (ii): %i"%ii,NIBS2+4)
    sprint("estimation of the explored percentage    (ii/nn): %.1f%%"%explored,NIBS2+4)
    sprint()
#==================================================#




#==================================================#
# All the strings printed in this program          #
#==================================================#
def print_welcome(PROGNAME):
    line = " Welcome to %s "%(PROGNAME.split(".py")[0].upper())
    print_head(line)
#--------------------------------------------------#
def print_inpcreation(IFILE):
    sprint("Input file: %s"%IFILE,NIBS)
    sprint("- Creating file...",NIBS2)
    sprint("- Please, modify this file before executing again!",NIBS2)
    sprint()
#--------------------------------------------------#
def print_notfound(IFILE):
    sprint("File name: %s [NOT FOUND]"%IFILE,NIBS)
    sprint()
#--------------------------------------------------#
def print_found(filename):
    sprint("File name: %s [FOUND]"%filename,NIBS)
    sprint()
#--------------------------------------------------#
def print_sthwrong(filename):
    sprint("Something is wrong with file: %s"%filename,NIBS)
    sprint()
#--------------------------------------------------#
def print_molinfo(mformu,natoms,ttorsions,tatoms,symbols,num_dummy):
    sprint("Molecular formula     : %s"%mformu,NIBS+3)
    sprint("Number of atoms       : %i"%natoms,NIBS+3)
    sprint("Number of dummy atoms : %i"%num_dummy,NIBS+3)
    sprint("Number of torsions    : %i"%len(ttorsions),NIBS+3)
    for X in ttorsions:
        atoms = tatoms[X]
        torsion_ats   = "-".join([symbols[at]             for at in atoms])
        torsion_atswn = "-".join([symbols[at]+"%i"%(at+1) for at in atoms])
        sprint("- dtor%-2s --> %s (%s)"%(X,torsion_ats,torsion_atswn),NIBS+4)
    sprint()
#--------------------------------------------------#
def print_readerror(variable):
    sprint("ERROR reading input file!",NIBS)
    sprint("variable '%s' in not (or wrongly) defined!"%variable,NIBS+3)
    sprint()
#--------------------------------------------------#
def print_varchanged(variable,val1,val2):
    # deal with boolean
    if val1 is True : val1 = "yes"
    if val2 is True : val2 = "yes"
    if val1 is False: val1 = "no"
    if val2 is False: val2 = "no"
    # print info
    sprint("Variable '%s' changed with regard to previous execution!"%variable,NIBS)
    sprint("  previous: "+str(val2),NIBS)
    sprint("  current : "+str(val1),NIBS)
    sprint()
#--------------------------------------------------#
def print_vars(inpvars):
    lines = inpvars.string_vars().split("\n")
    for line in lines:
        line = line.strip()
        if line == "": continue
        sprint("    "+line,NIBS)
    sprint()
#--------------------------------------------------#
def print_unknown_option(arg,PROGNAME):
    sprint("Unknown option: %s"%arg)
    sprint("Execute %s with -h to know which are the valid options!"%PROGNAME)
#--------------------------------------------------#
def print_torsiangleclass():
    sprint("Torsional angle classification (in degrees):",NIBS2)
    sprint()
    sprint("C  --> angle =   0"       ,NIBS2+4)
    sprint("C+ --> angle in (  0, 30]",NIBS2+4)
    sprint("G+ --> angle in ( 30, 75]",NIBS2+4)
    sprint("g+ --> angle in ( 75, 90]",NIBS2+4)
    sprint("a+ --> angle in ( 90,105]",NIBS2+4)
    sprint("A+ --> angle in (105,150]",NIBS2+4)
    sprint("T+ --> angle in (150,180)",NIBS2+4)
    sprint("T  --> angle = 180"       ,NIBS2+4)
    sprint("T- --> angle in (180,210)",NIBS2+4)
    sprint("A- --> angle in [210,255)",NIBS2+4)
    sprint("a- --> angle in [255,270)",NIBS2+4)
    sprint("g- --> angle in [270,285)",NIBS2+4)
    sprint("G- --> angle in [285,330)",NIBS2+4)
    sprint("C- --> angle in [330,360)",NIBS2+4)
    sprint()
    sprint("Ref: E. Papajak, P. Seal, X. Xu, and D. G. Truhlar,",NIBS2+4)
    sprint("     J. Chem. Phys. 137, 104314 (2012)",NIBS2+4)
    sprint()
#--------------------------------------------------#
def print_mjtable(lMj,lsj,points,Ntot):
    sprint("number of trial points: %i"%Ntot,NIBS2+4)
    max_sigma = max(lsj)
    sprint("max{sigma(Mj)} = %.4f"%max_sigma,NIBS2+4)
    sprint()
    ml = max(len(str(points[0])),len("Conformer"))
    tlineformat = " %%%is | %%6s | %%9s "%ml
    head = ("Conformer","Mj","sigma(Mj)")
    thead = tlineformat%head
    tdivi = "-"*len(thead)
    sprint(tdivi,NIBS2+4)
    sprint(thead,NIBS2+4)
    sprint(tdivi,NIBS2+4)
    lines = []
    for idx,geom in enumerate(points):
        data = (geom,"%6.4f"%lMj[idx],"%6.4f"%lsj[idx])
        lines.append( (lMj[idx],tlineformat%data) )
    lines.sort()
    for mj,line in lines: sprint(line,NIBS2+4)
    sprint(tdivi,NIBS2+4)
    sprint()
#--------------------------------------------------#
def print_notvisited(notvisited):
    sprint("Conformers not visited during Monte-Carlo:",NIBS2+4)
    for geom in notvisited: sprint(geom,NIBS2+6)
    sprint("These conformers will be ignored for MsTor",NIBS2+4)
    sprint()
#--------------------------------------------------#
def print_msgGauTempl(fname):
    sprint("Creating file with Gaussian templates: %s"%fname,NIBS)
    sprint("")
    sprint("When modifying this file, take into",NIBS2)
    sprint("account that:",NIBS2)
    sprint("  * Gaussian command line must start with #p",NIBS2)
    sprint("  * The iop(99/9=1,99/14=3) is needed",NIBS2)
    sprint("  * scf=(incore) is recommended when using Hartree-Fock",NIBS2)
    sprint()
#--------------------------------------------------#
def print_msgGauTempl_old(nfg,DIRTEMPL):
    if nfg == 0: return
    sprint("Creating Gaussian templates...",NIBS)
    sprint("Number of generated templates in '%s' is: %i"%(DIRTEMPL,nfg),NIBS2)
    sprint("The user should modify these files and should",NIBS2)
    sprint("take into account that:",NIBS2)
    sprint("  * Gaussian command line must start with #p",NIBS2)
    sprint("  * The iop(99/9=1,99/14=3) is needed",NIBS2)
    sprint("  * scf=(incore) is recommended when using Hartree-Fock",NIBS2)
    sprint()
#--------------------------------------------------#
def print_nogauexe(keyGau):
    sprint('* Gaussian path must be exported with the "%s"'%keyGau,NIBS+6)
    sprint('  environment variable in your .bashrc file!',NIBS+6)
    sprint('  e.g. export GauExe="/home/programs/Gaussian16/g16"',NIBS+6)
    sprint()
#--------------------------------------------------#
def get_nomenclature(vec,nomenclatures):
    nomenclature = vec.nomenclature()
    idx = 1
    while nomenclature in nomenclatures:
        idx += 1
        nomenclature = vec.nomenclature()+"_%i"%idx
    return nomenclature
#--------------------------------------------------#
def print_table_conformers(data,inpvars,fscal,EXCLUDED,extra,npre,nsto,SPACE="   "):
    '''
    Print table with information of each conformer
    '''

    # Columns in table
    COLS = ["","angles","conf","weight","E","E+ZPE"]

    # Conformer of minimum V0
    minV0, idx0 = sorted([(data_i[1],idx) for idx,data_i in enumerate(data)])[0]
    vec0 = data[idx0][0]
    # Conformer of minimum V1
    minV1, idx1 = sorted([(data_i[2],idx) for idx,data_i in enumerate(data)])[0]
    vec1 = data[idx1][0]

    # Some integers for table formatting
    iii   = len(str(len(data)))
    jjj   = max(len(str(vec0)),len(COLS[1]))
    kkk   = max(inpvars._ntorsions*2+2,len(COLS[2]))

    #--------------------#
    # Reference Energies #
    #--------------------#
    string  = "Conformer of minimum total energy\n"
    string += SPACE+"point: %s\n"%str(vec0)
    string += SPACE+"value: %.7f hartree\n"%minV0
    string += "\n"
    string += "Conformer of minimum total energy + zero-point energy\n"
    string += SPACE+"point: %s\n"%str(vec1)
    string += SPACE+"value: %.7f hartree\n"%minV1
    string += "\n"
    string += ""
    string += "Table with relative energies\n\n"
    string += SPACE+"weight : weight of the conformer (1 or 2)\n"
    string += SPACE+"E      : relative total energy (kcal/mol)\n"
    string += SPACE+"E+ZPE  : relative total energy plus zero-point energy (kcal/mol);\n"
    string += SPACE+"         [scale factor for frequencies is %.4f]\n"%fscal
    if inpvars._ts:
       string += SPACE+"ifreq : transition state imaginary frequency in cm^-1\n"
    string += "\n"

    #---------------#
    # HEAD OF TABLE #
    #---------------#
    LENS = [iii,jjj,kkk,6,7,7]
    if inpvars._ts: LENS += [8]; COLS += ["ifreq"]
    COLS = tuple(COLS)
    LENS = tuple(LENS)
    COLS = " "+" | ".join(fill_string(col,length) for col,length in zip(COLS,LENS))+" "
    DIVI = "-"*len(COLS)
    string += SPACE+DIVI+"\n"
    string += SPACE+COLS+"\n"
    string += SPACE+DIVI+"\n"

    #--------------------#
    # Add lines to TABLE #
    #--------------------#
    nomenclatures = []
    for idx,data_i in enumerate(data):
        vec,V0,V1,weight,rotcons,ifreq,Qrv,G,lzmat,zmatvals,fccards,gts = data_i

        # Relative energies
        relV0 = (V0-minV0)*KCALMOL
        relV1 = (V1-minV1)*KCALMOL
        # nomenclature
        nomenclature = get_nomenclature(vec,nomenclatures)
        nomenclatures.append(nomenclature)
        # Data for line
        args  = [fill_string(str(idx+1),iii)]
        args += [fill_string(str(vec),jjj)]
        args += [fill_string(nomenclature,kkk)]
        args += ["%2i"%weight]
        args += ["%7.3f"%relV0,"%7.3f"%relV1]
        if inpvars._ts:
           try   : sifreq = "%6.1fi"%abs(afreq2cm(ifreq))
           except: sifreq = "   -   "
           args += [sifreq]
        line = " "+" | ".join(fill_string(col,length) \
                         for col,length in zip(args,LENS))
        # preconditioned or stochastic?
        if gts.startswith("stoc"): line += "   "
        else                     : line += " * "
        # Info about checking
        if   gts in EXCLUDED["wrongconn"  ]: line += " <-- wrong connectivity"
        elif gts in EXCLUDED["outofdomain"]: line += " <-- invalid domain"
        elif gts in EXCLUDED["norestr"    ]: line += " <-- constraint not fulfilled"
        elif gts in EXCLUDED["wimag"      ]: line += " <-- wrong imag. freq."
        elif gts in EXCLUDED["repeated"   ]:
           (same,enantio,VECj) = extra[idx]
           if same   : line += " <-- same as %s"%str(VECj)
           if enantio: line += " <-- enantio. of %s"%str(VECj)

        # add line to string
        string += SPACE+line+"\n"
    string += SPACE+DIVI+"\n"
    string += SPACE+"* From preconditioned search\n"
    string += "\n"

    #---------------------------------------#
    # Counting Preconditioned vs Stochastic #
    #---------------------------------------#
    nconfs = npre+nsto
    string += "Number of conformers\n"
    string += "\n"
    string += SPACE+"total         : %3i\n"%nconfs
    string += SPACE+"preconditioned: %3i (%6.2f%%)\n"%(npre,100.0*npre/nconfs)
    string += SPACE+"stochastic    : %3i (%6.2f%%)\n"%(nsto,100.0*nsto/nconfs)

    #---------------------------------#
    # Print table with initial spaces #
    #---------------------------------#
    for line in string.split("\n"): print(IBS2+line)
#--------------------------------------------------#
def print_table_rotcons(data,SPACE="   "):
    # Columns in table
    COLS = ("","angles","B1","B2","B3")

    # Some integers for nice format
    iii   = len(str(len(data)))
    jjj   = max(len(str(data[0][0])),len(COLS[1]))

    # Prepare column names
    LENS = (iii,jjj,7,7,7,7)
    COLS = " "+" | ".join(fill_string(col,length) for col,length in zip(COLS,LENS))+" "
    HEAD = "-"*len(COLS)

    string  = "Table with rotational constants (Bi, in GHz)\n"
    string += "\n"
    string += HEAD+"\n"
    string += COLS+"\n"
    string += HEAD+"\n"
    for idx,data_i in enumerate(data):
        args = [str(idx+1),str(data_i[0])]
        #if dipole is None: dipole = "   -   "
        #else             : dipole = "%7.3f"%dipole
        rotcons = sorted(data_i[4])
        if rotcons is None: rotcons = ["   -   "  for  i in range(3)]
        else              : rotcons = ["%7.3f"%bi for bi in rotcons]
        args += rotcons
        line = " "+" | ".join([fill_string(col,length) for col,length in zip(args,LENS)])
        string += line+"\n"
    string += HEAD+"\n"
    for line in string.split("\n"):print(IBS2+SPACE+line)
#--------------------------------------------------#
def print_partition_function(data,temps,fscal,enantio,SPACE="   "):

    minV0, idx0 = sorted([(data_i[1],idx) for idx,data_i in enumerate(data)])[0]
    minV1, idx1 = sorted([(data_i[2],idx) for idx,data_i in enumerate(data)])[0]
    vec0 = data[idx0][0]
    vec1 = data[idx1][0]

    tmp  = [sorted([(data_i[7][idxT],i) for i,data_i in enumerate(data)])[0] \
            for idxT,T in enumerate(temps)]
    idxG = [idx for G,idx in tmp]
    minG = [ G  for G,idx in tmp]
    del tmp


    #-----------------------------------------#
    # Total ro-vibrational partition function #
    #-----------------------------------------#
    QRV  = np.zeros(len(temps))
    beta = 1.0/np.array(temps)/KB
    for idx,data_i in enumerate(data):
        point,V0,V1,weight,rotcons,ifreq,Qrv,gibbs,lzmat,zmatvals,sfccards,gts = data_i
        QRV += weight * Qrv * np.exp(-(V1-minV1)*beta)

    #-----------------------------#
    # Contribution to p. function #
    #-----------------------------#
    xs = []
    sum_prec = np.zeros(len(temps))
    sum_stoc = np.zeros(len(temps))
    for idx,data_i in enumerate(data):
        point,V0,V1,weight,rotcons,ifreq,Qrv,gibbs,lzmat,zmatvals,sfccards,gts = data_i
        x_i = Qrv * np.exp(-(V1-minV1)*beta) / QRV
        if   gts.startswith("prec."): sum_prec += weight * x_i
        elif gts.startswith("stoc."): sum_stoc += weight * x_i
        xs.append(x_i)

    #-------------------#
    # Total p. function #
    #-------------------#
    print_mshotable(temps,QRV,sum_prec,sum_stoc,fscal)

    #----------------#
    # Now, the table #
    #----------------#

    # Columns in table
    COLS = ("","angles","Qrv","X","G")

    # Some integers for nice format
    iii   = len(str(len(data)))
    jjj   = max(len(str(data[0][0])),len(COLS[1]))

    # Prepare column names
    LENS = (iii,jjj,11,7,7)
    COLS = " "+" | ".join(fill_string(col,length) for col,length in zip(COLS,LENS))+" "
    HEAD = "-"*len(COLS)

    string  = "INDIVIDUAL RRHO rovibrational partition functions\n"
    string += "\n"
    string += "  * Qrv : rovib pf of the conformer with regards to  its E + ZPE\n"
    string += "  * X   : contribution of the conformer to the total partition function\n"
    string += "  * G   : relative Gibbs free energy (in kcal/mol)\n"
    string += "\n"
    if enantio:
       string += "  --> Qrv and G do NOT include the conformational weight due to enantiomerism\n"
       string += "\n"
    for line in string.split("\n"):print(IBS2+line)


    for idxT,T in enumerate(temps):

        string  = HEAD+"\n"
        string += fill_string("T = %.3f K"%T,len(HEAD))+"\n"
        string += HEAD+"\n"
        string += COLS+"\n"
        string += HEAD+"\n"
        #print("#------------#")
        #print("# %8.2f K #"%T)
        #print("#------------#")
        for ii,data_i in enumerate(data):
            point  = data_i[0]
            V0     = data_i[1]
            V1     = data_i[2]
            weight = data_i[3]
            Qrv    = data_i[6][idxT]
            Gibbs  = (data_i[7][idxT] - minG[idxT]) * KCALMOL
            gts    = data_i[11]
            xi     = weight*xs[ii][idxT]

            args = (str(ii+1),str(point),"%11.3E"%Qrv,"%.5f"%xi,"%7.3f"%Gibbs)
            line = " "+" | ".join([fill_string(col,length) for col,length in zip(args,LENS)])
            string += line+"\n"
        string += HEAD+"\n"
        string += "\n"
        for line in string.split("\n"):print(IBS2+SPACE+line)

    return
#--------------------------------------------------#
def print_mshotable(temps,QMSHO,x_prec,x_stoc,freqscal):
    sprint("TOTAL MSHO rovibrational partition function:",NIBS2)
    sprint()
    sprint("- Calculated from vibrational zero-point energy",NIBS2+4)
    sprint("- Scale parameter for harmonic frequencies: %.4f"%freqscal,NIBS2+4)
    sprint()
    sprint("-----------------------------------------------------",NIBS2+6)
    sprint("  T(K)  |   Qrv(MSHO)  | % precondit. | % stochastic ",NIBS2+6)
    sprint("-----------------------------------------------------",NIBS2+6)
    for idx,T in enumerate(temps):
        datainline = (T,QMSHO[idx],x_prec[idx]*100,x_stoc[idx]*100)
        sprint("%7.2f | %12.4E |   %7.3f    |   %7.3f    "%datainline,NIBS2+6)
    sprint("-----------------------------------------------------",NIBS2+6)
    sprint()
#==================================================#


