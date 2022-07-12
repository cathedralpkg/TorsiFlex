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
| Sub-module :  optHLOPT           |
| Last Update:  2022/07/12 (Y/M/D) |
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
from   common.files      import write_zmat,write_gtsfile, read_zmat
#--------------------------------------------------#
import modtorsiflex.printing    as pp
import modtorsiflex.tfvars      as tvars
import modtorsiflex.tfgau       as itf
import modtorsiflex.tfhelper    as tfh
import modtorsiflex.tfrw        as rw
import modtorsiflex.version     as version
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

    #============================================#
    # Print information + Templates single point #
    #============================================#
    if not os.path.exists(tvars.GAUTEMPL):
       pp.print_template_notfound()
       raise exc.END
    templates = itf.read_templates()
    if inpvars._ts: keyoptHL,keyfrqHL = "tsopthl" ,"tsfrqhl"
    else          : keyoptHL,keyfrqHL = "minopthl","minfrqhl"
    try   : loptHL = templates[keyoptHL]
    except: pp.print_error_template(keyoptHL); raise exc.END
    try   : lfrqHL = templates[keyfrqHL]
    except: pp.print_error_template(keyfrqHL); raise exc.END
    pp.print_template(loptHL,keyoptHL,"HL optimization"         ,tvars.NIBS2)
    pp.print_template(lfrqHL,keyfrqHL,"HL frequency calculation",tvars.NIBS2)

    # temporal folder
    if not os.path.exists(inpvars._tmphl):
        try   : mkdir_recursive(inpvars._tmphl)
        except: pp.print_createdir(inpvars._tmphl,tvars.NIBS2); exit()


    # Read data from LL folder
    data, symbols = rw.folder_data(dirll,inpvars,verbose=False)
    # Sort according energy
    data.sort(key=lambda x:x[1])
    sprint("  number of low-level conformers   : %i"%len(data),tvars.NIBS2,1)

    #---------------------#
    # Create HL gjf files #
    #---------------------#
    nopt  = 0
    nfrq  = 0
    nsp   = 0
    minV0 = data[0][1]

    for idx,data_i in enumerate(data):
        
        count    = idx+1
        vecA     = data_i[0]
        V0       = data_i[1]
        lzmat    = data_i[8]
        zmatvals = data_i[9]
        sfccards = data_i[10]
        gts      = data_i[11]

        if not (lowerconformer <= count <= upperconformer): continue

        # read LL/HL-correlation file
        correlations = rw.read_llhlcorr()

        if   gts.split(".")[0] == tvars.PREF3p[:-1] :
            prefix1, prefix2, prefix3 = tvars.PREF1p, tvars.PREF2p, tvars.PREF3p
        elif gts.split(".")[0] == tvars.PREF3s[:-1] :
            prefix1, prefix2, prefix3 = tvars.PREF1s, tvars.PREF2s, tvars.PREF3s
        else: raise Exception

     #  # LL log file
     #  log_prec = "prec.%s.log"%svec
     #  log_stoc = "stoc.%s.log"%svec
     #  if   os.path.exists(inpvars._dirll+log_prec):
     #       log     = log_prec
     #       prefix2 = "frqprec."
     #       prefix3 = "prec."
     #  elif os.path.exists(inpvars._dirll+log_stoc):
     #       log     = log_stoc
     #       prefix1 = "optstoc."
     #       prefix2 = "frqstoc."
     #       prefix3 = "stoc."
     #  else: raise Exception
     #  # read log file
     #  conftuple,symbols,sfccards = itf.log_data(log,inpvars,inpvars._freqscalLL,dirll)
     #  # unpack conftuple
     #  vecA,V0,V1,G,weight,Qrv,ifreq,lzmat,zmatvals,log = conftuple

        sprint("(%i) %s [%s]"%(count,str(vecA),prefix3[:-1]),tvars.NIBS2)
        sprint("")

        # ONLY CONFORMERS IN ENERGY WINDOWS
        relV0_kcalmol = (V0-minV0)*KCALMOL
        sprint("LL electronic energy: %.4f kcal/mol"%relV0_kcalmol,tvars.NIBS2+4)

        # ONLY CONFORMERS BELOW HLCUTOFF IN GIBBS
        if inpvars._hlcutoff is not None:
           if relV0_kcalmol > inpvars._hlcutoff:
              sprint("hlcutoff    : %.3f kcal/mol"%inpvars._hlcutoff,tvars.NIBS2+4)
              sprint("HL optimization will not be carried out!",tvars.NIBS2+4)
              sprint()
              continue
           sprint()

        # SKIP IF ALREADY CALCULATED
        if str(vecA) in correlations.keys():
           vecB = correlations[str(vecA)]
           # possible names for gts files
           saved_gts1 = inpvars._dirhl+"%s.%s.gts"%(tvars.PREF3p,str(vecB))
           saved_gts2 = inpvars._dirhl+"%s.%s.gts"%(tvars.PREF3s,str(vecB))
           # check existence
           bool1 = os.path.exists(saved_gts1)
           bool2 = os.path.exists(saved_gts2)
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
        data = (loptHL,prefix1+str(vecA),inpvars,lzmat,zmatvals,"HL",sfccards)
        if calculate:
           sprint("optimization...",tvars.NIBS2+4)
           # execute Gaussian
           ofileOPT,statusOPT,dummy,vecB,zmatvals,zmatatoms,basic_opt = itf.execute(*data)
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

        bools,extra = tfh.precheck_geom(lzmat,zmatvals,cmatrix,inpvars)
        if (inpvars._tests[1][0] == 1) and not bools[0]:
           sprint("- connectivity test is negative!",tvars.NIBS+8,1)
           continue
        if (inpvars._tests[1][1] == 1) and not bools[1]:
           sprint("- torsion angle out of domain (torsion%i)!"%(extra[1]+1),tvars.NIBS+8,1)
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
           for ext in "gts,zmat,molden".split(","):
               # possible names for log files
               saved_file1 = inpvars._dirhl+"%s.%s.%s"%(tvars.PREF3p,str(equalto),ext)
               saved_file2 = inpvars._dirhl+"%s.%s.%s"%(tvars.PREF3s,str(equalto),ext)
               # check existence
               bool1 = os.path.exists(saved_file1)
               bool2 = os.path.exists(saved_file2)
               # if this is prec and the one obtained is stoc, replace name!!
               if bool2 and tvars.PREF3p in gts: move(saved_file2,saved_file1)
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
              statusFRQ,vecB_frq,xx,xx,xx,basic_frq = itf.read_normal_log(ofile,inpvars,"HL")
              ofileFRQ = ofile

        # Perform HL frequency
        else:
           sprint("frequency calculation...",tvars.NIBS2+4)
           ofileFRQ,dummy,statusFRQ,vecB_frq,dummy,dummy,basic_frq = itf.execute(*data)


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

        #============#
        # Save files #
        #============#
        if bool1 or bool2:
           # create z-matrix file
           file1 = inpvars._dirhl+prefix3+"%s.zmat"%str(vecB)
           if not os.path.exists(file1): write_zmat(file1,lzmat, zmatvals)
           # Create file with main info
           file2 = inpvars._dirhl+prefix3+"%s.gts"%str(vecB)
           if not os.path.exists(file2): write_gtsfile(*tuple(list(basic_frq)+[file2]))
           # create molden file
           file3 = inpvars._dirhl+prefix3+"%s.molden"%str(vecB)
           if not os.path.exists(file3): rw.gts2molden(file2,file3)
           # save LL-->HL correlation
           correlations[str(vecA)] = str(vecB)
           rw.add_llhlcorr(str(vecA),str(vecB))
           # Update counter
           nsp += 1

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

