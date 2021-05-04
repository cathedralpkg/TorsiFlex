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
| Sub-module :  optSEARCH          |
| Last Update:  2021/02/21 (Y/M/D) |
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
def gen_guess_zmatrix(vec,i_zmatvals,inpvars,closest):
    # read z-matrix from closest domain
    use_default = True
    if closest is not None:
       logfileS = inpvars._dirll+"prec.%s.log"%str(closest) # preconditioned
       logfileR = inpvars._dirll+"stoc.%s.log"%str(closest) # random
       if   os.path.exists(logfileS): logfile = logfileS
       elif os.path.exists(logfileR): logfile = logfileR
       else                         : logfile = None
       if logfile is None: use_default = True
       else:
          try:
             logdata   = gau.read_gaussian_log(logfile)
             zmatlines = logdata[10]
             zmat, zmatvals, symbols = rw.data_from_zmat_lines(zmatlines)
             sprint("using z-matrix from stored file (%s)"%logfile,tvars.NIBS2+4)
             use_default = False
          except: use_default = True
    # get closest z-matrix
    if use_default:
       sprint("using default z-matrix (%s)"%inpvars._zmatfile,tvars.NIBS2+4)
       zmatvals = {k:v for (k,v) in i_zmatvals.items()}
    # Apply vector to current z-matrix
    for angle,ic in zip(vec._fvec,inpvars._tic):
        zmatvals[ic] = angle
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

    #====================================#
    # Print information before searching #
    #====================================#
    if inpvars._ts: file1,file2 = tvars.TEMPLTSOPTLL,tvars.TEMPLTSFRQLL
    else          : file1,file2 = tvars.TEMPLMINOPTLL,tvars.TEMPLMINFRQLL
    # do templates exist??
    if not os.path.exists(file1): pp.print_filenotfound(file1,tvars.NIBS2,1); raise exc.END
    if not os.path.exists(file2): pp.print_filenotfound(file2,tvars.NIBS2,1); raise exc.END
    # read templates
    with open(file1,'r') as asdf: loptLL = asdf.readlines()
    with open(file2,'r') as asdf: lfrqLL = asdf.readlines()
    pp.print_infoLLsearch(inpvars,file1,file2,loptLL,lfrqLL)

    # Check gaussian software
    end_program = itf.check_executable()
    if end_program: return

    pp.print_tests(inpvars._tests)
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
       prefix1,prefix2,prefix3 = "optprec.","frqprec.","prec."
       if inpvars._prec != (1,1):
          sprint(" - Number of blocks: %i"%inpvars._prec[0],tvars.NIBS2)
          sprint(" - Selected block  : %i"%inpvars._prec[1],tvars.NIBS2)
       sprint()
    else:
       sprint("Stochastic algorithm will be used!",tvars.NIBS2,1)
       prefix1,prefix2,prefix3 = "optstoc.","frqstoc.","stoc."

    # save initial zmat
    i_lzmat = list(lzmat)
    i_zmatvals = {k:v for (k,v) in zmatvals.items()}
    # Calculate for each random angle
    count = 0
    nopt  = 0
    nfrq  = 0
    nsp   = 0
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
        sprint("preparing z-matrix...",tvars.NIBS2+4)
        bool_guess, closest = tfh.test_similarity_redundancy(vecA,inpvars,ddomains)
        if bool_guess == -1:
           sprint("already listed in %s"%tvars.FDOMAINS,tvars.NIBS2+4)
           sprint()
           continue
        if bool_guess == 0 and (inpvars._tests[0][1] == 1):
           sprint("similarity test is negative: in the domain of %s"%str(closest),tvars.NIBS2+4)
           sprint()
           continue
        args = (vecA,i_zmatvals,inpvars,closest)
        zmatvals = gen_guess_zmatrix(*args)
        if inpvars._enantio:
           vecA_enantio = tfh.enantiovec(lzmat,zmatvals,dcorr,inpvars)
           sprint("enantiomer: %s"%str(vecA_enantio),tvars.NIBS2+4)

        #- - - - - - - - - #
        # 1st set of TESTS #
        #- - - - - - - - - #
        bools = tfh.precheck_geom(lzmat,zmatvals,cmatrix,inpvars)
        if (inpvars._tests[0][0] == 1) and (bools[0] is False):
           sprint("connectivity test is negative: different connectivity!",tvars.NIBS+7,1)
           continue
        if (bools[1] is False):
           sprint("torsion angle out of domain...",tvars.NIBS+7,1)
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
        ofileOPT,statusOPT,dummy,vecB,zmatvals = itf.execute(*data)
        nopt += 1

        #===========================#
        # if NH2, check inversion!! #
        #===========================#
        if len(lNH2) != 0:
           xcc0 = intl.zmat2xcc(lzmat,zmatvals)
           inversion = False
           for HNRH_atoms,HNRH_value,HNRH_bool in lNH2:
               # check inversion
               cvalue,cbool = tfh.deal_with_HNRH(HNRH_atoms,xcc0)
               if HNRH_bool == cbool: continue
               # EXCHANGE HIDROGEN ATOMS!!
               inversion = True
               idxH1     = HNRH_atoms[0]
               idxH2     = HNRH_atoms[3]
               x1        = xcc0[3*idxH1:3*idxH1+3]
               x2        = xcc0[3*idxH2:3*idxH2+3]
               xcc0[3*idxH1:3*idxH1+3] = x2
               xcc0[3*idxH2:3*idxH2+3] = x1
               # recalculate zmatrix values (only those involving H1 or H2)
               for ic,atoms in zmatatoms.items():
                   if not (idxH1 in atoms or idxH2 in atoms): continue
                   xs = (xcc0[3*at:3*at+3] for at in atoms)
                   if len(atoms) == 2: zmatvals[ic] = ANGSTROM * fncs.distance(*xs)
                   if len(atoms) == 3: zmatvals[ic] = np.rad2deg(fncs.angle(*xs)   )
                   if len(atoms) == 4: zmatvals[ic] = np.rad2deg(fncs.dihedral(*xs))
           # recalculate vector
           if inversion:
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

        # optimization was carried out --> save guess
        if statusOPT != 0:
           rw.add_domain("fail",vecA,None,-1)
          #if inpvars._enantio: rw.add_domain("fail",vecA_enantio,None,-1)
           sprint()
           continue
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
        bools = tfh.precheck_geom(lzmat,zmatvals,cmatrix,inpvars)
        if (inpvars._tests[1][0] == 1) and not bools[0]:
           sprint("connectivity test is negative: different connectivity!",tvars.NIBS+7,1)
           rw.add_domain("excl",vecA,vecB,-1)
           if inpvars._enantio: rw.add_domain("excl",vecA_enantio,vecB_enantio,-1)
           continue
        if not bools[1]:
           sprint("torsion angle out of domain...",tvars.NIBS+7,1)
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
        ofileFRQ,dummy,statusFRQ,vecB_frq,dummy = itf.execute(*data)
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
        rw.add_domain(seldomain,vecA,vecB,statusFRQ)
        if   inpvars._enantio and seldomain == "conf":
           rw.add_domain("enan",vecA_enantio,vecB_enantio,statusFRQ)
        elif inpvars._enantio and seldomain == "wimag":
           rw.add_domain("wimag",vecA_enantio,vecB_enantio,statusFRQ)
        #============#
        # Save files #
        #============#
        if bool_conf:
           # copy frq file
           src = ofileFRQ
           dst = inpvars._dirll+prefix3+"%s.log"%str(vecB)
           if not os.path.exists(dst): copyfile(src, dst)
           # update counter
           nsp += 1
        sprint()

    # Print number of calculations
    pp.print_numcalcs(count,nopt,nfrq,nsp)
#==================================================#

