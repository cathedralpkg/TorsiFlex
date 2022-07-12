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
| Sub-module :  optSEARCH          |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''


#----------------------------------------------------------#
#          Importing basic and personal libraries          #
#----------------------------------------------------------#
import numpy as np
#----------------------------------------------------------#
import os
from   shutil import copyfile
#----------------------------------------------------------#
import common.gaussian          as gau
import common.internal          as intl
import common.fncs              as fncs
import common.Exceptions        as exc
from   common.physcons          import ANGSTROM
from   common.files             import write_zmat,write_gtsfile, read_zmat
#----------------------------------------------------------#
import modtorsiflex.printing    as pp
import modtorsiflex.tfvars      as tvars
import modtorsiflex.tfgau       as itf
import modtorsiflex.tfhelper    as tfh
import modtorsiflex.tfrw        as rw
from   modtorsiflex.printing    import sprint
from   modtorsiflex.tpespoint   import TorPESpoint
#----------------------------------------------------------#


#==================================================#
def gen_guess_zmatrix(vec,i_zmatvals,inpvars,zmat,distance):
    # read z-matrix from closest domain
    use_default = True
    try:
       lzmat,zmatvals,zmatatoms = read_zmat(inpvars._dirll+zmat)[0]
       sprint("using Z-matrix from stored file (%s)"%zmat,tvars.NIBS2+4)
       sprint("distance to stored Z-matrix: %.0f degrees"%distance,tvars.NIBS2+4)
       use_default = False
    except: use_default = True
    # get closest z-matrix
    if use_default:
       sprint("using default Z-matrix (%s)"%inpvars._zmatfile,tvars.NIBS2+4)
       zmatvals = {k:v for (k,v) in i_zmatvals.items()}
    # Apply vector to current z-matrix
    for angle,ic in zip(vec._fvec,inpvars._tic): zmatvals[ic] = angle
    # return zmatvals
    return zmatvals
#--------------------------------------------------#
def search_conformers(inpvars,zmat,symbols,cmatrix,dcorr,lNH2):
    
    ntorsions = len(inpvars._ttorsions)
    lzmat,zmatvals,zmatatoms = zmat
    # reference torsional vector
    vecR = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)
    inpvars._vecR = vecR

    # Create folders
    try   : tfh.create_dir(inpvars._tmpll)
    except: raise exc.END
    try   : tfh.create_dir(inpvars._dirll)
    except: raise exc.END

    #============================================#
    # Print information + Templates single point #
    #============================================#
    if not os.path.exists(tvars.GAUTEMPL):
       pp.print_template_notfound()
       raise exc.END
    templates = itf.read_templates()
    if inpvars._ts: keyoptLL,keyfrqLL = "tsoptll" ,"tsfrqll"
    else          : keyoptLL,keyfrqLL = "minoptll","minfrqll"
    try   : loptLL = templates[keyoptLL]
    except: pp.print_error_template(keyoptLL); raise exc.END
    try   : lfrqLL = templates[keyfrqLL]
    except: pp.print_error_template(keyfrqLL); raise exc.END
    pp.print_infoLLsearch(inpvars,keyoptLL,keyfrqLL,loptLL,lfrqLL)
    

    # Check gaussian software
    end_program = itf.check_executable()
    if end_program: return

    pp.print_tests(inpvars._tests)
    pp.print_ifreqconstr(inpvars._ifqrangeLL)

    #=======================#
    # Deal with domain file #
    #=======================#
    ddomains = rw.read_domains()
    if os.path.exists(tvars.FDOMAINS):
        sprint("Domain file exists (%s)!"%tvars.FDOMAINS,tvars.NIBS2,1)
        pp.print_domains(ddomains,tvars.DOMAINS)
    else: rw.initialize_domains(ntorsions)

    #=========================#
    # Search in torsional PES #
    #=========================#
    if inpvars._prec  is not False:
       sprint("Preconditioned algorithm will be used!",tvars.NIBS2)
       prefix1,prefix2,prefix3 = tvars.PREF1p, tvars.PREF2p, tvars.PREF3p
       if inpvars._prec != (1,1):
          sprint(" - Number of blocks: %i"%inpvars._prec[0],tvars.NIBS2)
          sprint(" - Selected block  : %i"%inpvars._prec[1],tvars.NIBS2)
       sprint()
    else:
       sprint("Stochastic algorithm will be used!",tvars.NIBS2,1)
       prefix1,prefix2,prefix3 = tvars.PREF1s, tvars.PREF2s, tvars.PREF3s

    # save initial zmat
    i_lzmat = list(lzmat)
    i_zmatvals = {k:v for (k,v) in zmatvals.items()}
    # Calculate for each random angle
    count = 0
    nopt  = 0
    nfrq  = 0
    nsp   = 0
    svec2zmat = tfh.get_zmat_vecs(inpvars,dcorr,svec2zmat={})
    for vecA in tfh.yield_angles(inpvars):

        count += 1
        # read domains
        ddomains = rw.read_domains()
        # print target point
        sprint("(%i) %s"%(count,str(vecA)),tvars.NIBS2)
        sprint()

        #------------------#
        # analyze geometry #
        #------------------#
        sprint("preparing Z-matrix...",tvars.NIBS2+4)
        args = (vecA,inpvars,ddomains,svec2zmat,dcorr)
        bool_guess, closest, svec2zmat, distance = tfh.test_similarity_redundancy(*args)
        if bool_guess == -1:
           sprint("already listed in %s"%tvars.FDOMAINS,tvars.NIBS2+4)
           sprint()
           continue
        if bool_guess == 0 and (inpvars._tests[0][1] == 1):
           sprint("similarity test is negative: in the domain of %s"%str(closest),tvars.NIBS2+4)
           sprint()
           continue
        closest_zmat = svec2zmat.get(str(closest),None)
        args = (vecA,i_zmatvals,inpvars,closest_zmat,distance)
        zmatvals = gen_guess_zmatrix(*args)
        if inpvars._enantio:
           vecA_enantio = tfh.enantiovec(lzmat,zmatvals,dcorr,inpvars)
           sprint("enantiomer: %s"%str(vecA_enantio),tvars.NIBS2+4)

        #- - - - - - - - - #
        # 1st set of TESTS #
        #- - - - - - - - - #
        bools,extra = tfh.precheck_geom(lzmat,zmatvals,cmatrix,inpvars)
        if (inpvars._tests[0][0] == 1) and (bools[0] is False):
           sprint("connectivity test is negative: different connectivity!",tvars.NIBS+7,1)
           continue
        if (bools[1] is False):
           sprint("torsion angle out of domain (torsion%i)..."%(extra[1]+1),tvars.NIBS+7,1)
           continue
        if (inpvars._tests[0][2] == 1) and (bools[2] is False):
           sprint("hard constraint test is negative!",tvars.NIBS+7,1)
           continue
        if (inpvars._tests[0][3] == 1) and (bools[3] is False):
           sprint("soft constraint test is negative!",tvars.NIBS+7,1)
           continue
        sprint()

        #------------------#
        # optimization     #
        #------------------#
        sprint("optimization...",tvars.NIBS2+4)
        data = (loptLL,prefix1+str(vecA),inpvars,lzmat,zmatvals,"LL")
        ofileOPT,statusOPT,dummy,vecB,zmatvals,zmatatoms,basic_opt = itf.execute(*data)
        nopt += 1

        # optimization was carried out --> save guess
        if statusOPT != 0:
           rw.add_domain("fail",vecA,None,-1)
          #if inpvars._enantio: rw.add_domain("fail",vecA_enantio,None,-1)
           sprint()
           continue

        #===========================#
        #    check NH2 inversion    #
        #===========================#
        zmatvals,inversion = tfh.correct_NH2_inversion(lzmat,zmatvals,zmatatoms,lNH2)
        if inversion:
           # recalculate vector
           svecB1 = str(vecB)
           vecB   = TorPESpoint(tfh.zmat2vec(zmatvals,inpvars._tic),inpvars._tlimit)
           svecB2 = str(vecB)
           if svecB1 != svecB2:
              sprint("NH2 inversion detected! Exchanging H atoms...",tvars.NIBS2+4)
              sprint("final vector: %s"%svecB2,tvars.NIBS2+4)
        #===========================#

        if inpvars._enantio:
           vecB_enantio = tfh.enantiovec(lzmat,zmatvals,dcorr,inpvars)
           sprint("enantiomer  : %s"%str(vecB_enantio),tvars.NIBS2+4)

        change = int(round(vecA.distance(vecB)))
        sprint("distance from guess point: %i degrees"%change,tvars.NIBS2+4)
        sprint()

        #------------------#
        # analyze geometry #
        #------------------#
        sprint("analysing geometry...",tvars.NIBS2+4)
        #- - - - - - - - - #
        # 2nd set of TESTS #
        #- - - - - - - - - #
        bools,extra = tfh.precheck_geom(lzmat,zmatvals,cmatrix,inpvars)
        if (inpvars._tests[1][0] == 1) and not bools[0]:
           sprint("connectivity test is negative: different connectivity!",tvars.NIBS+7,1)
           rw.add_domain("excl",vecA,vecB,-1)
           if inpvars._enantio: rw.add_domain("excl",vecA_enantio,vecB_enantio,-1)
           continue
        if not bools[1]:
           sprint("torsion angle out of domain (torsion%i)..."%(extra[1]+1),tvars.NIBS+7,1)
           rw.add_domain("excl",vecA,vecB,-1)
           if inpvars._enantio: rw.add_domain("excl",vecA_enantio,vecB_enantio,-1)
           continue
        if (inpvars._tests[1][2] == 1) and not bools[2]:
           sprint("hard constraint test is negative!",tvars.NIBS+7,1)
           rw.add_domain("excl",vecA,vecB,-1)
           if inpvars._enantio: rw.add_domain("excl",vecA_enantio,vecB_enantio,-1)
           continue
        if (inpvars._tests[1][3] == 1) and not bools[3]:
           sprint("soft constraint test is negative!",tvars.NIBS+7,1)
           rw.add_domain("excl",vecA,vecB,-1)
           if inpvars._enantio: rw.add_domain("excl",vecA_enantio,vecB_enantio,-1)
           continue
        # current conformations
       #confs = tfh.folder_points(inpvars._dirll)
        confs  = [ fvec for ivec,fvec,frqstatus in ddomains["conf"]]
        confs += [ fvec for ivec,fvec,frqstatus in ddomains["enan"]]
        # already stored?
        inlist, equalto  = vecB.is_in_list(confs,inpvars._epsdeg)
        if inlist and (inpvars._tests[1][1] == 1):
           sprint("redundancy test is negative: same as stored point (%s)"%equalto,tvars.NIBS2+4,1)
           rw.add_domain("repe",vecA,vecB,-1)
           if inpvars._enantio: rw.add_domain("repe",vecA_enantio,vecB_enantio,-1)
           continue
        sprint("new geometry found!",tvars.NIBS2+4)
        sprint()

        #------------------#
        # freq calculation #
        #------------------#
        sprint("frequency calculation...",tvars.NIBS2+4)
        data = (lfrqLL,prefix2+str(vecA),inpvars,lzmat,zmatvals,"LL")
        ofileFRQ,dummy,statusFRQ,vecB_frq,dummy,dummy,basic_frq = itf.execute(*data)
        nfrq += 1

        #------------------------------#
        # Save points into domain file #
        #------------------------------#
        bool1 = (statusFRQ == 0 and not inpvars._ts)
        bool2 = (statusFRQ == 1 and     inpvars._ts)
        bool_conf = bool1 or bool2
        # Save optimized point
        if bool_conf: seldomain = "conf"
        else        : seldomain = "wimag"

        # check imag freq
        if seldomain == "conf" and inpvars._ifqrangeLL != []:
           ifreq = abs( fncs.afreq2cm(itf.get_imag_freq(ofileFRQ,inpvars._freqscalLL)) )
           isok  = fncs.float_in_domain(ifreq,inpvars._ifqrangeLL)
           if not isok:
              seldomain,bool_conf = "wimag", False
              pp.print_excluded_ifreq(ifreq)
           else: pp.print_accepted_ifreq(ifreq)

        rw.add_domain(seldomain,vecA,vecB,statusFRQ)
        if   inpvars._enantio and seldomain == "conf":
           rw.add_domain("enan",vecA_enantio,vecB_enantio,statusFRQ)
        elif inpvars._enantio and seldomain == "wimag":
           rw.add_domain("wimag",vecA_enantio,vecB_enantio,statusFRQ)
        #============#
        # Save files #
        #============#
        if bool_conf:
           # create z-matrix file
           file1 = inpvars._dirll+prefix3+"%s.zmat"%str(vecB)
           if not os.path.exists(file1): write_zmat(file1,lzmat, zmatvals)
           # Create file with main info
           file2 = inpvars._dirll+prefix3+"%s.gts"%str(vecB)
           if not os.path.exists(file2): write_gtsfile(*tuple(list(basic_frq)+[file2]))
           # create molden file
           file3 = inpvars._dirll+prefix3+"%s.molden"%str(vecB)
           if not os.path.exists(file3): rw.gts2molden(file2,file3)
           # update counter
           nsp += 1
        sprint()

    # Print number of calculations
    pp.print_numcalcs(count,nopt,nfrq,nsp)
#==================================================#

