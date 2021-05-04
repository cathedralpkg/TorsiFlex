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
| Sub-module :  tfgau              |
| Last Update:  2020/12/21 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains things
related to Gincore,aussian
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
#--------------------------------------------------#
import modtorsiflex.tfvars        as tvars
import modtorsiflex.tfrw     as rw
import modtorsiflex.tfhelper      as tfh
from   modtorsiflex.tpespoint   import TorPESpoint
#==================================================#


KEYGAU = "GauExe"

#============================#
#   TEMPLATES FOR GAUSSIAN   #
#============================#
gt_opt_min = '''%nproc=[nproc]
%mem=[mem]
#p [level]
scf=(verytight)
iop(99/9=1,99/14=3)
opt=([optmode],MaxCycles=200)

--Optimization of minimum--

[charge],[multipl]
[zmat]
[modred]
[fccards]

'''
#----------------------------#
gt_opt_ts = '''%nproc=[nproc]
%mem=[mem]
#p [level]
scf=(verytight)
iop(99/9=1,99/14=3)
opt=([optmode],calcfc,ts,noeigentest,MaxCycles=200)

--Optimization of transition state--

[charge],[multipl]
[zmat]
[modred]
[fccards]

'''
#----------------------------#
gt_freq = '''%nproc=[nproc]
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
def generate_templates():
    # Create folder for templates
    mkdir_recursive(tvars.DIRTEMPL)
    # write templates
    ngen = 0
    for file_i, string_i in [(tvars.TEMPLMINOPTLL,gt_opt_min),\
                             (tvars.TEMPLMINOPTHL,gt_opt_min),\
                             (tvars.TEMPLMINFRQLL,gt_freq   ),\
                             (tvars.TEMPLMINFRQHL,gt_freq   ),\
                             (tvars.TEMPLTSOPTLL , gt_opt_ts),\
                             (tvars.TEMPLTSOPTHL , gt_opt_ts),\
                             (tvars.TEMPLTSFRQLL , gt_freq  ),\
                             (tvars.TEMPLTSFRQHL , gt_freq  )]:
        # check existency of files
        if os.path.exists(file_i): continue
        # create file
        with open(file_i,'w') as asdf: asdf.write(string_i)
        ngen += 1
    #--------------------#
    return ngen
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
def fccards_from_log(log):
    gcc,Fcc = gau.read_gaussian_log(log)[6:8]
    return gau.get_fccards_string(gcc,Fcc)
#--------------------------------------------------#
def generate_gjf(template,i_mname,inpvars,lzmat,zmatvals,case,log_ll=None):
    # Name of files/folders according to case
    if   case == "LL":
         dirtmp  = inpvars._tmpll
         dirsave = inpvars._dirll
         level   = inpvars._lowlevel
         log_ll  = None
    elif case == "HL":
         dirtmp  = inpvars._tmphl
         dirsave = inpvars._dirhl
         level   = inpvars._highlevel
    else             : raise Exception
    # variables of interest from inpvars
    nproc   = "%i"%inpvars._nproc
    mem     = "%s"%inpvars._mem
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
    if inpvars._fccards and (log_ll is not None):
       try   : fccards = fccards_from_log(log_ll)
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
    if   case == "LL": keylvl  = inpvars._keyll
    elif case == "HL": keylvl  = inpvars._keyhl
    else             : raise Exception
    # reading data
    exception_line = "ERROR! Level '%s' not found in archive section!\n"%keylvl
    try: logdata = gau.read_gaussian_log(ofile,target_level=keylvl)
    except LevelNotFound: print(" "*11+exception_line); sys.exit()
    commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energy,statusFRQ,ozmat,level = logdata
    # it may happen that in FRQ calculation zmatrix is not found!
    try:
       zmatvals = gau.convert_zmat(ozmat)[0][1]
       # put angles in (0,360) (sometimes Gaussian returns an angle > 360)
       for ic in inpvars._tic: zmatvals[ic] %= 360
       # Save info
       vec  = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)
    except:
       zmatvals = None
       vec      = None
    print("           energy (%s) = %.6f hartree"%(level,energy))
    if statusFRQ >= 0: print("           num imag frq: %i"%statusFRQ)
    print("           final vector: %s"%str(vec))
    return statusFRQ,vec,zmatvals
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
    if statusOPT == 0: statusFRQ,vec,zmatvals = read_normal_log(ofile,inpvars,case)
    # Update tracking file
    return (ofile,statusOPT,statusFRQ,vec,zmatvals)
#--------------------------------------------------#
def status_from_log(logfile,cmatrix,inpvars):
    '''
    About statusOPT, it returns:
       * -1 if the calculation failed
       *  0 if the calculation succeeded
    About statusFRQ, it returns:
       * -1       if the frequencies were not found in output (or fails)
       * num_imag if the frequencies were found
    '''
    assert os.path.exists(logfile)

    # Log file does not end with Normal termination?
    if not gau.normal_termination(logfile):
       statusOPT = -1
       statusFRQ = -1
       vec       = None
    # Ends with normal termination
    else:
       logdata   = gau.read_gaussian_log(logfile)
       statusOPT = 0
       statusFRQ = logdata[ 9]
       zmatlines = logdata[10]
       if zmatlines is None:
          statusOPT = -1
          statusFRQ = -1
          vec       = None
       else:
          zmat,zmatvals,symbols = rw.data_from_zmat_lines(zmatlines)
          # put angles in (0,360) (sometimes Gaussian returns an angle > 360)
          for ic in inpvars._tic: zmatvals[ic] %= 360
          # final vector
          vec = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)
    # Return data
    return statusOPT,statusFRQ,vec
#==================================================#



#==================================================#
def prepare_log_data(logdata):
    ch      = logdata[2]
    mtp     = logdata[3]
    symbols = logdata[4]
    xcc     = logdata[5]
    Fcc     = logdata[7]
    V0      = logdata[8]
    zmat    = logdata[10]
    # MM freq files --> z-matrix not shown in output; so get it from input
    if zmat is None: zmat = gau.zmat_from_loginp(folder+log)
    # z-matrix?
    if zmat is None:
        print("WARNING: file %s does not contains z-matrix!"%(folder+log))
        raise Exception
    # data from zmat lines
    (lzmat,zmatvals,zmatatoms), symbols = gau.convert_zmat(zmat)
    # correct symbols (just in case)
    symbols,atonums = fncs.symbols_and_atonums(symbols)
    # Data without dummies
    symbols_wo,xcc_wo = fncs.clean_dummies(symbols,xcc=list(xcc))
    symbols_wo,Fcc_wo = fncs.clean_dummies(symbols,Fcc=list(Fcc))
    # return data
    return ch,mtp,V0,symbols,symbols_wo,xcc_wo,Fcc_wo,lzmat,zmatvals
#--------------------------------------------------#
def folder_data(folder,inpvars,fscal):
    if inpvars._ts: fmode = -1
    else          : fmode = 0
    # the log files
    logs = [fname for fname in os.listdir(folder) \
            if fname.endswith(".log")]
    # First preconditioned, then stochastic!
    logs_prec = [log for log in logs if "prec." in log]
    logs_stoc = [log for log in logs if "stoc." in log]
    logs = sorted(logs_prec)+sorted(logs_stoc)
    # Calculate partition functions at different temperatures
    for log in logs:
        name    = log.split(".")[1]
        point   = TorPESpoint(name,inpvars._tlimit)
        # Read log file
        logdata = gau.read_gaussian_log(folder+log)
        # prepare data
        try   : logdata = prepare_log_data(logdata)
        except: continue
        ch,mtp,V0,symbols,symbols_wo,xcc_wo,Fcc_wo,lzmat,zmatvals = logdata

        # generate Molecule instance
        molecule = Molecule()
        molecule.setvar(xcc=xcc_wo,Fcc=Fcc_wo,V0=V0)
        molecule.setvar(fscal=fscal)
        molecule.setvar(symbols=symbols_wo)
        molecule.setvar(ch=ch,mtp=mtp)
        molecule.prepare()
        molecule.setup()
        molecule.ana_freqs()

        # Calculate partition functions for the temperatures
        qtot, V1, qis = molecule.calc_pfns(inpvars._temps,fmode=fmode)
        qtr,qrot,qvib,qele = qis
        Qrv = np.array([xx*yy for xx,yy in zip(qrot,qvib)])
        # Calculate partition functions at target temperature
        qtot, V1, qis = molecule.calc_pfns([inpvars._temp],fmode=fmode)
        # remove rotsigma for Gibbs
        gibbs = V1-pc.KB*inpvars._temp*np.log(qtot[0]*molecule._rotsigma)

        # imaginary frequency
        if   len(molecule._ccimag) == 0: ifreq = None
        elif len(molecule._ccimag) == 1: ifreq = [ifreq for ifreq,ivec in list(molecule._ccimag)][0]
        else                           : ifreq = None

        # weight for this conformer
        if inpvars._enantio and molecule._pgroup.lower() == "c1": weight = 2
        else                                                    : weight = 1

        # append and yield data
        conf_tuple = (point,V0,V1,gibbs,weight,Qrv,ifreq,lzmat,zmatvals,log)
        yield (conf_tuple,symbols)
#==================================================#


