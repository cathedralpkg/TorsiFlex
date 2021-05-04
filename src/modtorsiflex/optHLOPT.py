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
| Sub-module :  optHLOPT           |
| Last Update:  2020/12/21 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''

#--------------------------------------------------#
#      Importing basic and personal libraries      #
#--------------------------------------------------#
import os
from   shutil import copyfile
#--------------------------------------------------#
import common.gaussian   as gau
from   common.physcons   import KCALMOL
from   common.files      import mkdir_recursive
import common.Exceptions as exc
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
def highlevel_reopt(inpvars,cmatrix,calculate=False):

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


    # read LL/HL-correlation file
    correlations = rw.read_llhlcorr(inpvars._dirll,inpvars._dirhl)

    # Read LL-log files
    readall = False
    if inpvars._hlcutoff is not None:
       tminG = read_minG()
       if tminG is None: readall = True
       else:
          minG = tminG[0]
          if abs(tminG[1]-inpvars._temp) >= 0.01: readall = True

    sprint("Reading log files inside '%s'..."%dirll,tvars.NIBS2)
    if readall:
       data      = [(ctuple,symbols) for ctuple,symbols in \
                              itf.folder_data(dirll,inpvars,inpvars._freqscalLL)]
       symbols   = data[0][1]
       dataconfs = [ctuple           for ctuple,symbols in data]
#      dataconfs.sort(key=lambda x:x[3])
       del data
       sprint("Number of LL-log files: %i"%len(dataconfs),tvars.NIBS2+4)
       if len(dataconfs) == 0: return
       minG  = min([data[3] for data in dataconfs])
       with open(tvars.MINGTXT,'w') as asdf: asdf.write("minG %.8f %.2f\n"%(minG,inpvars._temp))
    else:
       dataconfs = itf.folder_data(dirll,inpvars,inpvars._freqscalLL)
    sprint()

    #---------------------#
    # Create HL gjf files #
    #---------------------#
    nopt  = 0
    nfrq  = 0
    nsp   = 0
    count = 0
    for conftuple in dataconfs:
        count += 1
        if len(conftuple) == 2: conftuple,symbols = conftuple
        # unpack conftuple
        vecA,V0,V1,G,weight,Qrv,ifreq,lzmat,zmatvals,log = conftuple
        sprint("(%i) %s"%(count,str(vecA)),tvars.NIBS2)
        if str(vecA) in correlations.keys():
           sprint("already in %s --> %s"%(inpvars._dirhl,correlations[str(vecA)]),tvars.NIBS2+4)
           sprint()
           continue
        sprint()

        # preconditioned or stochastic??
        if log.startswith("stoc."): prefix1,prefix2,prefix3 = "optstoc.","frqstoc.","stoc."
        else                      : prefix1,prefix2,prefix3 = "optprec.","frqprec.","prec."

        if inpvars._hlcutoff is not None:
           relG = (G-minG)*KCALMOL
           sprint("Gibbs energy: %.3f kcal/mol"%relG,tvars.NIBS2+4)
           if relG > inpvars._hlcutoff:
              sprint("hlcutoff    : %.3f kcal/mol"%inpvars._hlcutoff,tvars.NIBS2+4)
              sprint("HL optimization will not be carried out!",tvars.NIBS2+4)
              sprint()
              continue
           sprint()
        sprint("guess vector: %s"%str(vecA),tvars.NIBS2+4)
        sprint()

        #-----------------#
        # OPT CALCULATION #
        #-----------------#
        # Perform HL optimization
        data = (loptHL,prefix1+str(vecA),inpvars,lzmat,zmatvals,"HL",dirll+log)
        if calculate:
           sprint("optimization...",tvars.NIBS2+4)
           # execute Gaussian
           ofileOPT,statusOPT,dummy,vecB,zmatvals = itf.execute(*data)
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
           zmatvals = gau.convert_zmat(ozmat)[0][1]
           for ic in inpvars._tic: zmatvals[ic] %= 360
           vecB = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)


        # The optimization succedeed!
        sprint("analysing geometry...",tvars.NIBS2+4)
        sprint("- optimized point: %s"%str(vecB),tvars.NIBS2+8,1)
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
        inlist, equalto  = vecB.is_in_list(pointsHL,inpvars._epsdeg)
        if inlist:
           sprint("already in %s (%s)"%(inpvars._dirhl,str(equalto)),tvars.NIBS2+4)
           sprint()
           # save LL-->HL correlation
           correlations[str(vecA)] = str(vecB)
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
              statusFRQ,vecB_frq,dummy = itf.read_normal_log(ofile,inpvars,"HL")
              ofileFRQ = ofile

        # Perform HL frequency
        else:
           sprint("frequency calculation...",tvars.NIBS2+4)
           ofileFRQ,dummy,statusFRQ,vecB_frq,dummy = itf.execute(*data)

        #============#
        # Save files #
        #============#
        bool1 = (statusFRQ == 0 and not inpvars._ts)
        bool2 = (statusFRQ == 1 and     inpvars._ts)
        if bool1 or bool2:
           nsp += 1
           # copy frq file
           src = ofileFRQ
           dst = inpvars._dirhl+prefix3+"%s.log"%str(vecB)
           if not os.path.exists(dst): copyfile(src, dst)
           # save LL-->HL correlation
           correlations[str(vecA)] = str(vecB)
        sprint()


    # Print number of calculations
    pp.print_numcalcs(count,nopt,nfrq,nsp)

    if not calculate:
       sprint("Number of opt gjf files generated: %i"%nopt,tvars.NIBS2)
       sprint("Number of frq gjf files generated: %i"%nfrq,tvars.NIBS2)
       sprint()
    else:
       rw.write_llhlcorr(correlations)
#==================================================#

