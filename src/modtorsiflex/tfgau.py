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
| Sub-module :  tfgau              |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains things
related to Gaussian
'''

#==================================================#
import os
import numpy as np
import sys
#--------------------------------------------------#
import common.fncs       as fncs
import common.physcons   as pc
import common.gaussian   as gau
from   common.internal   import zmat2xcc
from   common.Molecule   import Molecule
from   common.files      import mkdir_recursive
from   common.Exceptions import ExeNotDef, ExeNotFound, LevelNotFound
from   common.pgs        import get_pgs
#--------------------------------------------------#
import modtorsiflex.tfvars      as tvars
import modtorsiflex.tfrw        as rw
import modtorsiflex.tfhelper    as tfh
from   modtorsiflex.tpespoint   import TorPESpoint
#==================================================#


KEYGAU = "GauExe"

#============================#
#   TEMPLATES FOR GAUSSIAN   #
#============================#
gt_opt_min1 = '''%nproc=[nproc]
%mem=[mem]
#p [level]
scf=(tight)
iop(99/9=1,99/14=3)
opt=([optmode],tight,maxcycle=200)

--Optimization of minimum--

[charge],[multipl]
[zmat]
[modred]
[fccards]

'''
#----------------------------#
gt_opt_min2 = '''%nproc=[nproc]
%mem=[mem]
#p [level]
scf=(verytight)
iop(99/9=1,99/14=3)
opt=([optmode],verytight,maxcycle=200)

--Optimization of minimum--

[charge],[multipl]
[zmat]
[modred]
[fccards]

'''
#----------------------------#
gt_opt_ts1 = '''%nproc=[nproc]
%mem=[mem]
#p [level]
scf=(tight)
iop(99/9=1,99/14=3)
opt=([optmode],tight,calcfc,ts,noeigentest,maxcycle=200)

--Optimization of transition state--

[charge],[multipl]
[zmat]
[modred]
[fccards]

'''
#----------------------------#
gt_opt_ts2 = '''%nproc=[nproc]
%mem=[mem]
#p [level]
scf=(verytight)
iop(99/9=1,99/14=3)
opt=([optmode],verytight,calcfc,ts,noeigentest,maxcycle=200)

--Optimization of transition state--

[charge],[multipl]
[zmat]
[modred]
[fccards]

'''
#----------------------------#
gt_freq1 = '''%nproc=[nproc]
%mem=[mem]
#p [level]
scf=(tight)
iop(99/9=1,99/14=3)
freq=noraman

--Frequency calculation--

[charge],[multipl]
[zmat]

'''
#----------------------------#
gt_freq2 = '''%nproc=[nproc]
%mem=[mem]
#p [level]
scf=(verytight)
iop(99/9=1,99/14=3)
freq=noraman

--Frequency calculation--

[charge],[multipl]
[zmat]

'''
#============================#


#==================================================#
#          FUNCTIONS RELATED TO TEMPLATES          #
#==================================================#
def generate_templates(mode=1):
    '''
    '''
    #--------------------#
    def aux_templates(key,string_i):
        string  = ""
        string += "#--------------#\n"
        string += "start_%s\n"%key
        string += string_i
        string += "end_%s\n"%key
        string += "#~~~~~~~~~~~~~~#\n\n"
        return string
    #--------------------#
    ngen = 0
    #--------------------#
    if mode == 0:
       # Create folder for templates
       mkdir_recursive(tvars.DIRTEMPL)
       for file_i, string_i in [(tvars.TEMPLMINOPTLL,gt_opt_min1),\
                                (tvars.TEMPLMINOPTHL,gt_opt_min2),\
                                (tvars.TEMPLMINFRQLL,gt_freq1   ),\
                                (tvars.TEMPLMINFRQHL,gt_freq2   ),\
                                (tvars.TEMPLTSOPTLL , gt_opt_ts1),\
                                (tvars.TEMPLTSOPTHL , gt_opt_ts2),\
                                (tvars.TEMPLTSFRQLL , gt_freq1  ),\
                                (tvars.TEMPLTSFRQHL , gt_freq2  )]:
           # check existency of files
           if os.path.exists(file_i): continue
           # create file
           with open(file_i,'w') as asdf: asdf.write(string_i)
           ngen += 1
    #--------------------#
    if mode == 1 and not os.path.exists(tvars.GAUTEMPL):
       string  = aux_templates("MINOPTLL",gt_opt_min1)
       string += aux_templates("MINOPTHL",gt_opt_min2)
       string += aux_templates("MINFRQLL",gt_freq1   )
       string += aux_templates("MINFRQHL",gt_freq2   )
       string += aux_templates("TSOPTLL" ,gt_opt_ts1 )
       string += aux_templates("TSOPTHL" ,gt_opt_ts2 )
       string += aux_templates("TSFRQLL" ,gt_freq1   )
       string += aux_templates("TSFRQHL" ,gt_freq2   )
       # create file
       with open(tvars.GAUTEMPL,'w') as asdf: asdf.write(string)
       ngen += 1
    #--------------------#
    return ngen
#--------------------------------------------------#
def read_templates():
    template = {}
    if not os.path.exists(tvars.GAUTEMPL): return template
    # Read lines in file
    with open(tvars.GAUTEMPL,'r') as asdf: lines = asdf.readlines()
    # Extract data
    current = None
    for line in lines:
        if line.startswith("start_"):
           current = line.split("_")[1].strip().lower()
           template[current] = []
        elif line.startswith("end_"):
           current = None
        elif current is not None:
           template[current].append( line )
    return template
#==================================================#

#==================================================#
#--------------------------------------------------#
def check_executable():
    print(tvars.IBS+"   Checking Gaussian executable file...")
    try:
       gau.set_EXE()
       print(tvars.IBS+"      * %s: %s"%(gau.KEYGAU,gau.EXE))
       gau.check_EXE()
       end_program = False
    except ExeNotDef:
       print(tvars.IBS+"      * %s: Environment variable NOT DEFINED!!!"%(gau.KEYGAU))
       end_program = True
    except ExeNotFound:
       print(tvars.IBS+"        - file not found!")
       end_program = True
    print("")
    return end_program
#--------------------------------------------------#
def generate_gjf(template,i_mname,inpvars,lzmat,zmatvals,case,string_fccards=None):
    # Name of files/folders according to case
    if   case == "LL":
         dirtmp  = inpvars._tmpll
         dirsave = inpvars._dirll
         level   = inpvars._lowlevel
         nproc   = "%i"%inpvars._nprocll
         mem     = "%s"%inpvars._memll
    elif case == "HL":
         dirtmp  = inpvars._tmphl
         dirsave = inpvars._dirhl
         level   = inpvars._highlevel
         nproc   = "%i"%inpvars._nprochl
         mem     = "%s"%inpvars._memhl
    else             : raise Exception
    # variables of interest from inpvars
    charge  = "%i"%inpvars._charge
    multipl = "%i"%inpvars._multipl
    optmode = int(inpvars._optmode)
    tatoms  = inpvars._tatoms
    # Name of files
    wname, ifile, ofile, chk, fchk, err = gau.iofiles(i_mname,dirtmp)
    # (a) Log file already exists and ends normally (calculation is NOT needed)
    if os.path.exists(ofile) and gau.normal_termination(ofile):
       return False, ifile, ofile, err

    # (b) Input file must be created (calculation is needed)

    # (b.1) LL Hessian matrix 
    fccards = None
    if inpvars._fccards and (string_fccards is not None):
       try   : fccards = string_fccards
       except: fccards = None
    # (b.2) Write Gaussian input file
    args   = (ifile,template,level,charge,multipl,nproc,mem,lzmat,zmatvals,optmode,tatoms,fccards)
    string = write_gjf(*args)
    return True, ifile, ofile, err
#--------------------------------------------------#
def write_gjf(ifile,template,level,charge,multipl,nproc,mem,lzmat,zmatvals,optmode,tatoms,fccards=None):
    # Generate input string
    string = ""
    for line in template:
        # replace [optmode]
        if "[optmode]" in line:
            if optmode == 0: optstring = "z-matrix"
            else           : optstring = "modredundant"
            if fccards is not None: optstring += ",fccards"
            line = line.replace("[optmode]",optstring)
        # replace [nproc]
        if "[nproc]" in line: line = line.replace("[nproc]",nproc)
        # replace [mem]
        if "[mem]" in line: line = line.replace("[mem]",mem)
        # replace [level]
        if "[level]" in line: line = line.replace("[level]",level)
        # replace [charge]
        if "[charge]" in line: line = line.replace("[charge]",charge)
        # replace [multipl]
        if "[multipl]" in line: line = line.replace("[multipl]",multipl)
        # replace [zmat]
        if "[zmat]" in line:
            szmat = gau.gen_zmatrix_string(lzmat,zmatvals)
            line = line.replace("[zmat]",szmat)
        # replace [fccards]
        if "[fccards]" in line:
            if fccards is None: line = ""
            else: line = line.replace("[fccards]",fccards)
        # replace [modred]
        if "[modred]" in line:
            if optmode != 1: line = ""
            else:
               modredlines = ""
               for X,atoms in tatoms.items():
                   modredlines += "D %s B\n"%" ".join(["%i"%(at+1) for at in atoms])
               line = line.replace("[modred]",modredlines)
        # increase string
        string += line
    # Write file
    with open(ifile,'w') as asdf: asdf.write(string)
#--------------------------------------------------#
def read_normal_log(ofile,inpvars,case):
    # keylvl from inpvars according to case
    if   case is None: keylvl  = None
    elif case == "LL": keylvl  = inpvars._keyll
    elif case == "HL": keylvl  = inpvars._keyhl
    else             : raise Exception
    # reading data
    exception_line = "ERROR! Level '%s' not found in archive section!\n"%keylvl
    try: logdata = gau.read_gaussian_log(ofile,target_level=keylvl)
    except LevelNotFound: print(" "*11+exception_line); sys.exit()
    commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energy,statusFRQ,ozmat,level = logdata
    # rotational sigma
    symbols,atonums = fncs.symbols_and_atonums(symbols)
    masses          = fncs.symbols2masses(symbols)
    pgroup,rotsigma = get_pgs(atonums,masses,xcc)
    # it may happen that in FRQ calculation zmatrix is not found!
    try:
       lzmat,zmatvals,zmatatoms = gau.convert_zmat(ozmat)[0]
       # put angles in (0,360) (sometimes Gaussian returns an angle > 360)
       for ic in inpvars._tic: zmatvals[ic] %= 360
       # Save info
       vec  = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)
    except:
       zmatvals  = None
       zmatatoms = None
       vec       = None
    print("           energy (%s) = %.6f hartree"%(level,energy))
    print("           final vector: %s"%str(vec))
    if statusFRQ >= 0: print("           num imag frq: %i"%statusFRQ)
    # return data
    basic = (xcc,atonums,ch,mtp,energy,pgroup,rotsigma,gcc,Fcc)
    return statusFRQ,vec,lzmat,zmatvals,zmatatoms,basic
#--------------------------------------------------#
def execute(template,i_mname,inpvars,lzmat,zmatvals,case,log_ll=None):
    '''
    About statusOPT, it returns:
       * -1 Normal termination not found
       *  0 Normal termination found
    About statusFRQ, it returns:
       * -1       if the frequencies were not found in output (or fails)
       * num_imag if the frequencies were found
    '''
    # Setup Gaussian
    gau.set_EXE()
    gau.check_EXE()
    # Generate input file
    generated, ifile, ofile, err = generate_gjf(template,i_mname,inpvars,lzmat,zmatvals,case,log_ll)
    print("           input  file: %s"%ifile)
    print("           output file: %s"%ofile)
    #------------------#
    # Execute Gaussian #
    #------------------#
    if not generated:
       print("           output file already exists and ends with Normal termination!")
       statusOPT = 0
    else:
       try   : gau.execute(ifile,ofile,err)
       except KeyboardInterrupt: sys.exit()
       except: pass
       # Normal termination?
       if gau.normal_termination(ofile):
          print("           output file ends with Normal termination!")
          statusOPT = 0
       else:
          print("           output file does not end with Normal termination!")
          statusOPT = -1
    #-------------------#
    # Read Gaussian log #
    #-------------------#
    vec, statusFRQ = None, -1
    # Read output file
    if statusOPT == 0: statusFRQ,vec,lzmat,zmatvals,zmatatoms,basic = read_normal_log(ofile,inpvars,case)
    else             : zmatatoms,basic = None, None
    # return data
    return (ofile,statusOPT,statusFRQ,vec,zmatvals,zmatatoms,basic)
#--------------------------------------------------#
def status_from_log(logfile,inpvars,dcorr):
    '''
    About statusOPT, it returns:
       * -1 if the calculation failed
       *  0 if the calculation succeeded
    About statusFRQ, it returns:
       * -1       if the frequencies were not found in output (or fails)
       * num_imag if the frequencies were found
    '''
    assert os.path.exists(logfile)


    vec1,evec1 = None, None
    vec2,evec2 = None, None
    statusOPT  = -1
    statusFRQ  = -1

    # Input vector
    try:
       zmat = gau.zmat_from_loginp(logfile)
       (lzmat,zmatvals,zmatatoms), symbols = gau.convert_zmat(zmat)
       vec1 = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)
       if inpvars._enantio:
          evec1 = tfh.enantiovec(lzmat,zmatvals,dcorr,inpvars)
    except:
       vec1,evec1 = None, None

    # Log file does not end with Normal termination?
    if gau.normal_termination(logfile):
       logdata   = gau.read_gaussian_log(logfile)
       statusOPT = 0
       statusFRQ = logdata[ 9]
       zmatlines = logdata[10]
       if zmatlines is not None:
          (lzmat,zmatvals,zmatatoms), symbols = gau.convert_zmat(zmatlines)
          # put angles in (0,360) (sometimes Gaussian returns an angle > 360)
          for ic in inpvars._tic: zmatvals[ic] %= 360
          # final vector
          vec2 = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)
          if inpvars._enantio:
             evec2 = tfh.enantiovec(lzmat,zmatvals,dcorr,inpvars)
    # Return data
    return statusOPT,statusFRQ,(vec1,evec1),(vec2,evec2)
#==================================================#



#==================================================#
def prepare_log_data(logdata):
    ch      = logdata[2]
    mtp     = logdata[3]
    symbols = logdata[4]
    xcc     = logdata[5]
    gcc     = logdata[6]
    Fcc     = logdata[7]
    V0      = logdata[8]
    zmat    = logdata[10]
    # MM freq files --> z-matrix not shown in output; so get it from input
    if zmat is None: zmat = gau.zmat_from_loginp(folder+log)
    # z-matrix?
    if zmat is None:
        print("WARNING: file %s does not contains Z-matrix!"%(folder+log))
        raise Exception
    # data from zmat lines
    (lzmat,zmatvals,zmatatoms), symbols = gau.convert_zmat(zmat)
    # correct symbols (just in case)
    symbols,atonums = fncs.symbols_and_atonums(symbols)
    # Data without dummies
    symbols_wo,xcc_wo = fncs.clean_dummies(symbols,xcc=list(xcc))
    symbols_wo,gcc_wo = fncs.clean_dummies(symbols,gcc=list(gcc))
    symbols_wo,Fcc_wo = fncs.clean_dummies(symbols,Fcc=list(Fcc))
    # return data
    return ch,mtp,V0,symbols,symbols_wo,xcc_wo,gcc_wo,Fcc_wo,lzmat,zmatvals
#--------------------------------------------------#
def get_imag_freq(log,fscal):
    try:
        # Read log file
        logdata = gau.read_gaussian_log(log)
        logdata = prepare_log_data(logdata)
        ch,mtp,V0,symbols,symbols_wo,xcc_wo,gcc_wo,Fcc_wo,lzmat,zmatvals = logdata

        # diagonalize Fcc
        molecule = Molecule()
        molecule.setvar(xcc=xcc_wo,Fcc=Fcc_wo,V0=V0)
        molecule.setvar(fscal=fscal)
        molecule.setvar(symbols=symbols_wo)
        molecule.setvar(ch=ch,mtp=mtp)
        molecule.prepare()
        molecule.setup()
        molecule.ana_freqs()
        ifreq = molecule._ccimag[0][0]
    except: ifreq = None
    # return
    return ifreq
#==================================================#


