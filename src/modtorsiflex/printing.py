'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: TorsiFlex
Version     : 2021.3
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
| Sub-module :  printing           |
| Last Update:  2021/11/22 (Y/M/D) |
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
#--------------------------------------------------#
from   common.fncs          import time2human
from   modtorsiflex.tfvars import PROGVER,PROGNAME, PROGNAMEnopy, VERSION, DATE
from   modtorsiflex.tfvars import FDOMAINS, IFILE, DIRTEMPL
from   modtorsiflex.tfvars import NIBS, NIBS2, IBS, IBS2
#==================================================#

DESCRDOMAINS  = {\
"conf" :"Set of conformers",\
"enan" :"Set of enantiomers",\
"repe" :"Set of optimized points that match with a stored conformer",\
"excl" :"Set of optimized points that were excluded",\
"wimag":"Set of optimized points with a wrong number of imaginary frequencies",\
"fail" :"Set of guess points whose optimization failed",\
}

HSTRING = '''
Current version: %s %s

Description:

  A program to seek the conformers of flexible molecules
  by means of a combined preconditioned-stochastic algorithm.
  The located conformers, calculated with a low-level (LL)
  electronic structure method, can be re-optimized using
  a high-level (HL) method.

Main authors:

  Dr.   David Ferro-Costas
  Prof. Antonio Fernandez-Ramos

  Centro Singular de Investigacion en Quimica Bioloxica
  e Materiais Moleculares (CIQUS), Universidade de
  Santiago de Compostela, Galicia, Spain

Execution:

  %s [--help /-h]  [--version/-v]  [--inp/--input]
  %s [--prec    ]  [--stoc      ]  [--hlopt      ]
  %s [--msho    ]  [--mstor     ]  [--regen      ]

-------------------
  Program options  
-------------------

   --inp / --input

     Generates the input file and the templates for Gaussian.


   --prec [M m]

     Uses the preconditioned algorithm for the conformer location.

     In order to divide the preconditioned guesses into M groups
     and deal with the m-th group, use this option as follows:
          --hlopt M m
     For example:
          torsiflex --hlopt 10 2
     divides the guesses into 10 groups and only carries out
     the calculations associated to the 2nd group.

     --> Calculations with Gaussian are carried out <--


   --stoc

     Uses the stochastic algorithm for the conformer location.

     --> Calculations with Gaussian are carried out <--


   --hlopt [nocalc]

     Re-optimizes LL conformers at HL.

     When followed by 'nocalc':
          --hlopt nocalc
     it generates the gjf files (Gaussian inputs) without
     carrying out calculations.
     Useful to send the calculations on your own.

     --> Calculations with Gaussian are carried out <--


   --msho [ll/hl]

     Checks the located conformers and calculates the
     multi-structure harmonic-oscillator (MS-HO) partition
     functions.

     This option can be carried out exclusively for the
     low-level (ll) or the high-level (hl) conformers if
     followed by the corresponding abbreviation:
          --msho ll
          --msho hl


   --mstor [ll/hl]

     Generates the MsTor input files.

     This option can be carried out exclusively for the
     low-level (ll) or the high-level (hl) conformers if
     followed by the corresponding abbreviation:
          --mstor ll
          --mstor hl


   --regen

     Regenerates the %s file using the temporal files.


   --help (also -h)

     Prints this help message.


   --version (also -v)

     Prints the program version.

---------------------
  Extra Information  
---------------------

(a) Assert the path to the Gaussian executable
    is defined in your .bashrc file under the name
    'GauExe' and export it. Example:

    export GauExe='/home/programs/Gaussian/g09'

(b) Modify the Gaussian templates (in %s),
    taking into account that:

    * [nproc], [mem], [level], [optmode], [charge],
      [multipl], [zmat], [modred], and [fccards]
      are %s indications.
      They should not be removed.

    * Gaussian command line must start with '#p'

    * The command line includes 'iop(99/9=1,99/14=3)'.
      This is mandatory and should not be deleted:
         99/9=1  --> rotates to Z-matrix orientation first
         99/14=3 --> expresses final optimized structure
                     in terms of the input Z-matrix

    * 'scf=(incore)' is recommended with Hartree-Fock
      calculations. With it, Gaussian stores the full
      integral list in memory, speeding up the calculation.

(c) We highly recommend to firstly carry out the
    preconditioned search, which should found the
    chemically-intuitive conformers.
    After it, the stochastic algorithm should be used.
'''%(PROGNAMEnopy,PROGVER,PROGNAME," "*len(PROGNAME)," "*len(PROGNAME),FDOMAINS,DIRTEMPL,PROGNAMEnopy)



#==================================================#
def sprint(string="",nibs=0,nbl=0):
    print(" "*nibs+string+"\n"*nbl)
#==================================================#

#==================================================#
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
def print_template(template,filename,case,nibs):
    sprint("Template for %s (%s):"%(case,filename),nibs)
    sprint("---------------------------",nibs)
    for line in template: sprint(line[:-1],nibs)
    sprint("---------------------------",nibs)
    sprint()
#--------------------------------------------------#
def print_infoLLsearch(inpvars,TEMPLOPTLL,TEMPLFRQLL,loptLL,lfrqLL):
    # (a) tolerance
    sprint("Similarity test tolerance: %5.1f degrees"%inpvars._dist1d,NIBS2)
    sprint("Redundancy test tolerance: %5.1f degrees"%inpvars._epsdeg,NIBS2)
    sprint()
    # (b) LL opt template
    print_template(loptLL,TEMPLOPTLL,"LL optimization",NIBS2)
    # (c) LL freq template
    print_template(loptLL,TEMPLFRQLL,"LL frequency calculation",NIBS2)
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
def print_msgGauTempl(nfg,DIRTEMPL):
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
def print_mshotable(temps,QMSHO,QMSHO_pre,QMSHO_sto,freqscal):
    sprint("MSHO rovibrational partition function:",NIBS2)
    sprint()
    sprint("- Calculated from vibrational zero-point energy",NIBS2+4)
    sprint("- Scale parameter for harmonic frequencies: %.4f"%freqscal,NIBS2+4)
    sprint()
    sprint("-----------------------------------------------------",NIBS2+6)
    sprint("  T(K)  |   Qrv(MSHO)  | % precondit. | % stochastic ",NIBS2+6)
    sprint("-----------------------------------------------------",NIBS2+6)
    for idx,T in enumerate(temps):
        perc_pre = 100.0*QMSHO_pre[idx]/QMSHO[idx]
        perc_sto = 100.0*QMSHO_sto[idx]/QMSHO[idx]
        datainline = (T,QMSHO[idx],perc_pre,perc_sto)
        sprint("%7.2f | %12.4E |   %7.3f    |   %7.3f    "%datainline,NIBS2+6)
    sprint("-----------------------------------------------------",NIBS2+6)
    sprint()
#==================================================#


