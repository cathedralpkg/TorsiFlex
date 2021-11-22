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
| Sub-module :  optHLOPT           |
| Last Update:  2021/11/22 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''

#--------------------------------------------------#
#      Importing basic and personal libraries      #
#--------------------------------------------------#
import os
from   shutil import copyfile, move
#--------------------------------------------------#
import common.gaussian   as gau
import common.fncs       as fncs
import common.Exceptions as exc
from   common.physcons   import KCALMOL
from   common.files      import mkdir_recursive
from   common.files      import write_zmat
#--------------------------------------------------#
import modtorsiflex.printing    as pp
import modtorsiflex.tfvars        as tvars
import modtorsiflex.tfgau    as itf
import modtorsiflex.tfhelper      as tfh
import modtorsiflex.tfrw     as rw
from   modtorsiflex.printing    import sprint
from   modtorsiflex.tpespoint   import TorPESpoint
#--------------------------------------------------#




#==================================================#
# Refinement of geoms at high-level                #
#==================================================#
def read_minG():
    if not os.path.exists(tvars.MINGTXT): return None
    # read file
    with open(tvars.MINGTXT,'r') as asdf: lines = asdf.readlines()
    try:
      for line in lines:
          if "minG" not in line: continue
          dummy, minG, temp = line.split()
          return (float(minG),float(temp))
    except: pass
    return None
#--------------------------------------------------#
#def deal_with_ll_conformer(conftuple):
#--------------------------------------------------#
def highlevel_reopt(inpvars,cmatrix,\
                    lowerconformer=0,\
                    upperconformer=float("inf"),\
                    calculate=False,\
                    lNH2=[]):

    dirll, dirhl = inpvars._dirll , inpvars._dirhl
    if not os.path.exists(dirll):
       print("       Folder NOT FOUND: %s\n"%dirll)
       return
    # Create folders
    if not os.path.exists(dirhl):
       try   : mkdir_recursive(dirhl)
       except: pp.print_createdir(dirhl,tvars.NIBS2); exit()
    # Check Gaussian software
    if calculate:
       end_program = itf.check_executable()
       if end_program: return

    # print status of tests
    pp.print_tests(inpvars._tests,1)
    pp.print_ifreqconstr(inpvars._ifqrangeHL)

    #------------------------#
    # Templates for Gaussian #
    #------------------------#
    # HL opt & HL freq templates for MIN
    if not inpvars._ts: file1,file2 = tvars.TEMPLMINOPTHL,tvars.TEMPLMINFRQHL
    # HL opt & HL freq templates for TS
    else              : file1,file2 = tvars.TEMPLTSOPTHL,tvars.TEMPLTSFRQHL
    # (a) do templates exist??
    if not os.path.exists(file1): pp.print_filenotfound(file1,tvars.NIBS2,1); raise exc.END
    if not os.path.exists(file2): pp.print_filenotfound(file2,tvars.NIBS2,1); raise exc.END
    # (b) read templates
    with open(file1,'r') as asdf: loptHL = asdf.readlines()
    pp.print_template(loptHL,file1,"HL optimization",tvars.NIBS2)
    with open(file2,'r') as asdf: lfrqHL = asdf.readlines()
    pp.print_template(lfrqHL,file2,"HL frequency calculation",tvars.NIBS2)
    # (c) temporal folder
    if not os.path.exists(inpvars._tmphl):
        try   : mkdir_recursive(inpvars._tmphl)
        except: pp.print_createdir(inpvars._tmphl,tvars.NIBS2); exit()


    # read LL energies
    if not os.path.exists(tvars.ENERGYSUMLL):
        sprint("- file NOT found: %s"%tvars.ENERGYSUMLL,tvars.NIBS2)
        sprint("  execute %s using --msho ll"%tvars.PROGNAMEnopy,tvars.NIBS2)
        sprint("")
        return
    sprint("- reading file: %s"%tvars.ENERGYSUMLL,tvars.NIBS2)
    sprint("")
    llenergies,  lltemp = rw.read_llenergies()
    # check temperature
    if abs(lltemp-inpvars._temp) > 1e-3:
        sprint("  something went wrong!",tvars.NIBS2)
        sprint("")
        sprint("  temperature in file         : %.3f K"%(lltemp),tvars.NIBS2)
        sprint("  temperature in %s: %.3f K"%(tvars.IFILE,inpvars._temp),tvars.NIBS2)
        sprint("")
        return
    sprint("  temperature for Gibbs free energy: %.3f K"%lltemp,tvars.NIBS2)
    # check number of conformers
    logs = [fname for fname in os.listdir(inpvars._dirll) if fname.endswith(".log")]
    if len(logs) != len(llenergies):
        sprint("  something went wrong!",tvars.NIBS2)
        sprint("  number of conformer differs!",tvars.NIBS2)
        sprint("")
        sprint("  # of conformers in file     : %i"%len(llenergies),tvars.NIBS2)
        sprint("  # of conformers in LL folder: %i"%len(logs),tvars.NIBS2)
        sprint("")
        return
    sprint("  number of low-level conformers   : %i"%len(llenergies),tvars.NIBS2)
    sprint("")
    # check conformers
    svecs_llenergies = [svec for count,V0,G,svec in llenergies]
    for log in logs:
        svec_i = log.split(".")[1]
        if svec_i not in svecs_llenergies:
           sprint("  something went wrong!",tvars.NIBS2)
           sprint("")
           sprint("  conformer NOT included: %s"%log,tvars.NIBS2)
           sprint("")
           return

    #---------------------#
    # Create HL gjf files #
    #---------------------#
    nopt  = 0
    nfrq  = 0
    nsp   = 0
    minV0 = llenergies[0][0]
    for count,relV0_kcalmol,relG_kcalmol,svec in llenergies:
        if not (lowerconformer <= count <= upperconformer): continue
        # read LL/HL-correlation file
        correlations = rw.read_llhlcorr()
        # LL log file
        log_prec = "prec.%s.log"%svec
        log_stoc = "stoc.%s.log"%svec
        if   os.path.exists(inpvars._dirll+log_prec):
             log     = log_prec
             prefix1 = "optprec."
             prefix2 = "frqprec."
             prefix3 = "prec."
        elif os.path.exists(inpvars._dirll+log_stoc):
             log     = log_stoc
             prefix1 = "optstoc."
             prefix2 = "frqstoc."
             prefix3 = "stoc."
        else: raise Exception
        # read log file
        conftuple,symbols,string_fccards = itf.log_data(log,inpvars,inpvars._freqscalLL,dirll)
        # unpack conftuple
        vecA,V0,V1,G,weight,Qrv,ifreq,lzmat,zmatvals,log = conftuple
        sprint("(%i) %s [%s]"%(count,str(vecA),prefix3[:-1]),tvars.NIBS2)
        sprint("")


        # ONLY CONFORMERS IN ENERGY WINDOWS
        sprint("LL electronic energy: %.4f kcal/mol"%relV0_kcalmol,tvars.NIBS2+4)

        # ONLY CONFORMERS BELOW HLCUTOFF IN GIBBS
        if inpvars._hlcutoff is not None:
           relG = (G-minG)*KCALMOL
           sprint("Gibbs energy: %.3f kcal/mol"%relG,tvars.NIBS2+4)
           if relG > inpvars._hlcutoff:
              sprint("hlcutoff    : %.3f kcal/mol"%inpvars._hlcutoff,tvars.NIBS2+4)
              sprint("HL optimization will not be carried out!",tvars.NIBS2+4)
              sprint()
              continue
           sprint()

        # SKIP IF ALREADY CALCULATED
        if str(vecA) in correlations.keys():
           vecB = correlations[str(vecA)]
           # possible names for log files
           saved_log1 = inpvars._dirhl+"prec.%s.log"%str(vecB)
           saved_log2 = inpvars._dirhl+"stoc.%s.log"%str(vecB)
           # check existence
           bool1 = os.path.exists(saved_log1)
           bool2 = os.path.exists(saved_log2)
           # exists? then skip
           if bool1 or bool2:
              sprint("already in %s"%inpvars._dirhl,tvars.NIBS2+4)
              sprint("optimized conformer: %s"%str(vecB),tvars.NIBS2+4)
              sprint()
              continue
        sprint()

        #-----------------#
        # OPT CALCULATION #
        #-----------------#
        # Perform HL optimization
        data = (loptHL,prefix1+str(vecA),inpvars,lzmat,zmatvals,"HL",string_fccards)
        if calculate:
           sprint("optimization...",tvars.NIBS2+4)
           # execute Gaussian
           ofileOPT,statusOPT,dummy,vecB,zmatvals,zmatatoms = itf.execute(*data)
           sprint()
           # opt ends ok?
           if statusOPT == -1: continue
        else:
           sprint("optimization file:",tvars.NIBS2+4)
           generated, ifile, ofile, err = itf.generate_gjf(*data)
           sprint("- %s"%ifile,tvars.NIBS2+8)
           if generated:
              sprint("- generated!",tvars.NIBS2+8,1)
              nopt += 1
              continue
           sprint("- not generated!",tvars.NIBS2+8,1)
           logdata = gau.read_gaussian_log(ofile)
           commands,comment,ch,mtp,symbols,xcc,gcc,Fcc,energy,statusFRQ,ozmat,level = logdata
           lzmat,zmatvals,zmatatoms = gau.convert_zmat(ozmat)[0]
           for ic in inpvars._tic: zmatvals[ic] %= 360
           vecB = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)

        # The optimization succedeed!
        sprint("analysing geometry...",tvars.NIBS2+4)
        sprint("- optimized point: %s"%str(vecB),tvars.NIBS+8,1)

        # NH2 inversion
        zmatvals,inversion = tfh.correct_NH2_inversion(lzmat,zmatvals,zmatatoms,lNH2)
        if inversion:
           svecB1 = str(vecB)
           vecB   = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)
           svecB2 = str(vecB)
           if svecB1 != svecB2:
              sprint("- NH2 inversion detected! Exchanging H atoms...",tvars.NIBS+8)
              sprint("- final vector: %s"%svecB2,tvars.NIBS+8,1)

        bools = tfh.precheck_geom(lzmat,zmatvals,cmatrix,inpvars)
        if (inpvars._tests[1][0] == 1) and not bools[0]:
           sprint("- connectivity test is negative!",tvars.NIBS+8,1)
           continue
        if (inpvars._tests[1][1] == 1) and not bools[1]:
           sprint("- torsion angle out of domain!",tvars.NIBS+8,1)
           continue
        if (inpvars._tests[1][2] == 1) and not bools[2]:
           sprint("- hard constraint test is negative!",tvars.NIBS+8,1)
           continue
        if (inpvars._tests[1][3] == 1) and not bools[3]:
           sprint("- soft constraint test is negative!",tvars.NIBS+8,1)
           continue

        # Already saved?
        pointsHL = tfh.folder_points(inpvars._dirhl)
        inlist, equalto = vecB.is_in_list(pointsHL,inpvars._epsdeg)
        if inlist:
           # exists?
           sprint("already in %s (%s)"%(inpvars._dirhl,str(equalto)),tvars.NIBS2+4)
           sprint()
           # rename stoc --> prec if required!
           for ext in "log,zmat".split(","):
               # possible names for log files
               saved_log1 = inpvars._dirhl+"prec.%s.%s"%(str(equalto),ext)
               saved_log2 = inpvars._dirhl+"stoc.%s.%s"%(str(equalto),ext)
               # check existence
               bool1 = os.path.exists(saved_log1)
               bool2 = os.path.exists(saved_log2)
               # if this is prec and the one obtained is stoc, replace name!!
               if bool2 and "prec" in prefix3: move(saved_log2,saved_log1)
           # save LL-->HL correlation
           correlations[str(vecA)] = str(vecB)
           rw.add_llhlcorr(str(vecA),str(vecB))
           continue
        sprint()

        #-----------------#
        # FRQ CALCULATION #
        #-----------------#
        data = (lfrqHL,prefix2+str(vecA),inpvars,lzmat,zmatvals,"HL")
        if not calculate:
           sprint("frequency file:",tvars.NIBS2+4)
           generated, ifile, ofile, err = itf.generate_gjf(*data)
           sprint("- %s"%ifile,tvars.NIBS2+8)
           if generated:
              sprint("- generated!\n",tvars.NIBS2+8,1)
              nfrq += 1
              continue
           # log with normal termination!
           else:
              sprint("- not generated!\n",tvars.NIBS2+8,1)
              statusFRQ,vecB_frq,dummy1,dummy2 = itf.read_normal_log(ofile,inpvars,"HL")
              ofileFRQ = ofile

        # Perform HL frequency
        else:
           sprint("frequency calculation...",tvars.NIBS2+4)
           ofileFRQ,dummy,statusFRQ,vecB_frq,dummy,dummy = itf.execute(*data)

        #============#
        # Save files #
        #============#
        bool1 = (statusFRQ == 0 and not inpvars._ts)
        bool2 = (statusFRQ == 1 and     inpvars._ts)

        # check imag freq
        if bool2 and inpvars._ifqrangeHL != []:
           ifreq = abs( fncs.afreq2cm(itf.get_imag_freq(ofileFRQ,inpvars._freqscalHL)) )
           isok  = fncs.float_in_domain(ifreq,inpvars._ifqrangeHL)
           if not isok:
              sprint("imaginary frequency is: %.2fi cm^-1"%ifreq,tvars.NIBS2+4)
              bool2 = False

        if bool1 or bool2:
           nsp += 1
           # create z-matrix file
           dst = inpvars._dirhl+prefix3+"%s.zmat"%str(vecB)
           write_zmat(dst,lzmat, zmatvals)
           # copy frq file
           src = ofileFRQ
           dst = inpvars._dirhl+prefix3+"%s.log"%str(vecB)
           if not os.path.exists(dst): copyfile(src, dst)
           # save LL-->HL correlation
           correlations[str(vecA)] = str(vecB)
           rw.add_llhlcorr(str(vecA),str(vecB))

        sprint()


    # Print number of calculations
    pp.print_numcalcs(count,nopt,nfrq,nsp)

    if not calculate:
       sprint("Number of opt gjf files generated: %i"%nopt,tvars.NIBS2)
       sprint("Number of frq gjf files generated: %i"%nfrq,tvars.NIBS2)
       sprint()
    else:
       # rewrite correlations to remove repetitions
       correlations = rw.read_llhlcorr()
       rw.write_llhlcorr(correlations)
#==================================================#

