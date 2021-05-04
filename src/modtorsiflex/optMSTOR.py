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
| Sub-module :  optMSTOR           |
| Last Update:  2020/12/21 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''

#==================================================#
import os
import numpy as np
#--------------------------------------------------#
import common.internal   as intl
import common.gaussian   as gau
import common.fncs       as fncs
from   common.pgs        import get_pgs
#--------------------------------------------------#
import modtorsiflex.tfrw     as rw
import modtorsiflex.printing    as pp
import modtorsiflex.tfhelper      as tfh
from   modtorsiflex.tfvars        import NIBS2, DIRMSTOR
from   modtorsiflex.printing    import sprint
from   modtorsiflex.tpespoint   import TorPESpoint
#==================================================#


#==================================================#
# Functions for writing mstor files                #
#==================================================#
def clean_torsions_from_ics(ics,torsions):
    # central atoms of torsions
    central_bonds = [ tuple(sorted([b,c])) for (a,b,c,d) in torsions]
    # keep torsions without those central atoms
    ics_clean = []
    for ictype,ic in ics:
        if len(ic) == 4:
           at1,at2,at3,at4 = ic
           if tuple(sorted([at2,at3])) in central_bonds: continue
        ics_clean.append( (ictype,ic) )
    # Add torsions
    ics_clean += [("pt",(at1,at2,at3,at4)) for (at1,at2,at3,at4) in torsions]
    # return
    return ics_clean
#--------------------------------------------------#
def valid_redundant_set(list_xcc,list_gcc,list_Fcc,list_ccfrqs,symbols,masses,torsions):
    for trial in range(2):
        for idx1 in range(len(list_xcc)):
            xcc1    = list_xcc[idx1]
            gcc1    = list_gcc[idx1]
            Fcc1    = list_Fcc[idx1]
            ccfrqs1 = list_ccfrqs[idx1]
            ics     = intl.ics_from_geom(xcc1,symbols)
            if trial == 0: ics = intl.ics_depure(ics)
            # add target torsions and depure the rest
            ics     = clean_torsions_from_ics(ics,list(torsions))
            # check it is valid
            valid = True
            for idx2 in range(len(list_xcc)):
                if idx1 == idx2: continue
                xcc2    = list_xcc[idx2]
                gcc2    = list_gcc[idx2]
                Fcc2    = list_Fcc[idx2]
                ccfrqs2 = list_ccfrqs[idx2]
                # calc ic-freqs
                icfrqs2 = intl.calc_icfreqs(Fcc2,masses,xcc2,gcc2,ics)[0]
                # compare
                valid = fncs.same_freqs(ccfrqs2,icfrqs2)
                if not valid: break
            if valid: break
        if valid: break
    if valid: return ics
    else    : return None
#--------------------------------------------------#
def valid_nonredundant_set(ics,list_xcc,list_gcc,list_Fcc,list_ccfrqs,symbols,masses,torsions):
    if ics is None: return None
    ics_torsions = [("pt",tuple(torsion)) for torsion in  torsions]
    valid        = False
    for idx1 in range(len(list_xcc)):
        xcc1    = list_xcc[idx1]
        gcc1    = list_gcc[idx1]
        Fcc1    = list_Fcc[idx1]
        ccfrqs1 = list_ccfrqs[idx1]
        # generate non-redundant from ics
        nrics, same = intl.nonredundant(xcc1,masses,gcc1,Fcc1,ics,ccfrqs1,unremov=list(ics_torsions))
        if not same: continue
        # clean torsions
        nrics = clean_torsions_from_ics(ics,torsions)
        # check it is valid
        valid = True
        for idx2 in range(len(list_xcc)):
            if idx1 == idx2: continue
            xcc2    = list_xcc[idx2]
            gcc2    = list_gcc[idx2]
            Fcc2    = list_Fcc[idx2]
            ccfrqs2 = list_ccfrqs[idx2]
            # calc ic-freqs
            icfrqs2 = intl.calc_icfreqs(Fcc2,masses,xcc2,gcc2,nrics)[0]
            # compare
            valid = fncs.same_freqs(ccfrqs2,icfrqs2)
            if not valid: break
        if valid: break
    if valid: return nrics
    else    : return None
#--------------------------------------------------#
def calculate_mj(inpvars,strconfs,weights,allconfs,denantio):
    ttorsions = list(inpvars._ttorsions)
    ntor      = len(ttorsions)
    nconfs    = len(strconfs)
    # convert each point in points to string
    # initialize number of points at each point
    lNj     = np.array([        0     for conf in strconfs])
    lwj     = np.array([weights[conf] for conf in strconfs])
    
    # update frequency for (Mj,sigma_j) & max number of iterations
    #UPDATEFREQ, MAXITER = ntor  *10**3, ntor  *(2*10**5)
    UPDATEFREQ, MAXITER = nconfs*10**2, nconfs*10**5
    # Calculate Mj's
    Ntot = 0
    while True:
        Ntot += 1
        # generate random point
        rvec = tfh.random_point(inpvars._tdomain,inpvars._tlimit)
        # check which is closest structure
        closest = str(rvec.closest(allconfs))
        # enantiomer??
        closest = denantio.get(closest,closest)
        # locate point
        idx = strconfs.index(closest)
        # update frequency of the structure
        lNj[idx] += 1
        # Calculate Mj and sigma_j
        if Ntot % UPDATEFREQ == 0:
           # correct Nj with weight
           llNj  = [1.0*Nj/wj for Nj,wj in zip(lNj,lwj)]
           # calculate Mj and sigma(Mj)
           llMj  = [(1.0*Ntot/Nj)**(1.0/ntor) if Nj != 0 else np.inf for Nj in llNj]
           llsj  = [ Mj/ntor/np.sqrt(Nj) if Nj != 0 else 0 for Nj,Mj in zip(llNj,llMj)]
           # ignore infinity for max(sigma)
           max_sigma = max(llsj)
           sprint("iteration %6i --> max(sigma_j)= %5.3f"%(Ntot,max_sigma),NIBS2+4)
           if max_sigma <= inpvars._sigmamj: break

        # if too many, break anyway
        if Ntot == MAXITER: sprint("too many iterations...",NIBS2+4); break
    # print values
    sprint()
    pp.print_mjtable(llMj,llsj,strconfs,Ntot)
    # Any conformer not visited?
    notvisited = [strconfs[idx] for idx,Nj in enumerate(lNj) if Nj == 0]
    if len(notvisited) != 0: pp.print_notvisited(notvisited)
    # return data
    return {strconf:Mj for strconf,Mj in zip(strconfs,llMj)}, notvisited
#--------------------------------------------------#
def gen_mstor(inpvars,case,dcorr,lCH3=[]):
    if case == "ll":
        folder   = inpvars._dirll
        freqscal = inpvars._freqscalLL
    elif case == "hl":
        folder   = inpvars._dirhl
        freqscal = inpvars._freqscalHL

    # Folder with data does not exist
    if not os.path.exists(folder):
       pp.print_dirnotfound(folder,NIBS2,1)
       return
    # Create folder for mstor files
    tfh.create_dir(folder+DIRMSTOR)

    # list of points
    logs   = [log for log in os.listdir(folder)  if log.endswith(".log")]
    points = [TorPESpoint(log.split(".")[1]) for log in logs]

    # no files?
    if len(logs) == 0:
       sprint("There are no log files inside folder '%s'"%folder,NIBS2,1)
       return

    #----------------#
    # Read log files #
    #----------------#
    sprint("Reading log files...",NIBS2,1)
    geoms    = []
    weights  = {}
    strconfs = []
    allconfs = []
    denantio = {}
    for log,point in zip(logs,points):
        logdata   = gau.read_gaussian_log(folder+log,target_level=None)

        ch        = logdata[2]
        mtp       = logdata[3]
        symbols   = logdata[4]
        xcc       = logdata[5]
        gcc       = logdata[6]
        Fcc       = logdata[7]
        V0        = logdata[8]
        zmat      = logdata[10]

        if zmat is None:
           zmat = gau.zmat_from_loginp(folder+log)

        (lzmat,zmatvals,zmatatoms), symbols = gau.convert_zmat(zmat)

        # same current point
        strconfs.append(str(point))
        allconfs.append(point)
        # correct symbols (just in case)
        symbols,atonums = fncs.symbols_and_atonums(symbols)
        # Data without dummies
        symbols_wo,xcc_wo = fncs.clean_dummies(symbols,xcc=list(xcc))
        symbols_wo,gcc_wo = fncs.clean_dummies(symbols,gcc=list(gcc))
        symbols_wo,Fcc_wo = fncs.clean_dummies(symbols,Fcc=list(Fcc))

        # calculate weight
        atonums_wo = fncs.symbols2atonums(symbols_wo)
        masses_wo  = fncs.symbols2masses(symbols_wo)
        pgroup, rotsigma = get_pgs(atonums_wo,masses_wo,xcc_wo)
        if inpvars._enantio and pgroup.lower() == "c1":
           weight  = 2
           enantio = tfh.enantiovec(lzmat,zmatvals,dcorr,inpvars)
           allconfs.append( enantio)
           denantio[str(enantio)] = str(point)
        else:
           weight = 1
        # save data
        geoms.append( (V0,xcc_wo,gcc_wo,Fcc_wo,str(point)) )
        weights[str(point)] = weight
    geoms.sort()


    #----------------------------#
    # calculate M_{j,tau} values #
    #----------------------------#
    sprint("Calculating M_{j,tau} values using Monte-Carlo...",NIBS2)
    dMj,notvisited = calculate_mj(inpvars,strconfs,weights,allconfs,denantio)
    sprint()
    # clean geoms
    geoms = [geom_tuple for geom_tuple in geoms if geom_tuple[-1] not in notvisited]

    #----------------------#
    # Internal coordinates #
    #----------------------#
    nrics = None
    symbols = symbols_wo
    if True:
       # Generate non-redundant internal coordinates
       torsions    = [inpvars._tatoms[X] for X in inpvars._ttorsions]
       masses      = fncs.symbols2masses(symbols)
       list_xcc    = [fncs.shift2com(xcc,masses) for (V0,xcc,gcc,Fcc,point) in geoms]
       list_gcc    = [gcc                        for (V0,xcc,gcc,Fcc,point) in geoms]
       list_Fcc    = [Fcc                        for (V0,xcc,gcc,Fcc,point) in geoms]
       list_ccfrqs = [fncs.calc_ccfreqs(Fcc,masses,xcc)[0] for xcc,Fcc in zip(list_xcc,list_Fcc)]

       # Update torsions with CH3s
       torsions = torsions+lCH3
       sprint()
       # (a) generate redundant internal coordinates
       sprint("Trying to generate valid redundant set of internal coordinates...",NIBS2)
       nvdof = len(list_ccfrqs[0])
       sprint("number of vib. degrees of freedom: %i"%nvdof,NIBS2+4)
       idata_rics = (list_xcc,list_gcc,list_Fcc,list_ccfrqs,symbols,masses,torsions)
       rics       = valid_redundant_set(*idata_rics)
       nrics      = None
       num_ics    = None
       if rics is not None:
          num_ics = intl.count_ics(rics)
          sprint("number of internal coordinates   : %i"%num_ics,NIBS2+4)
          sprint()
       
       # (b) generate non-redundant internal coordinates
       if num_ics == nvdof: nrics = rics
       elif rics is not None:
          sprint("Trying to generate valid non-redundant set of internal coordinates...",NIBS2)
          sprint("number of vib. degrees of freedom: %i"%nvdof,NIBS2+4)
          idata_nrics = (rics,list_xcc,list_gcc,list_Fcc,list_ccfrqs,symbols,masses,torsions)
          nrics       = valid_nonredundant_set(*idata_nrics)
          if nrics is not None:
             num_ics = intl.count_ics(nrics)
             sprint("number of internal coordinates   : %i"%num_ics,NIBS2+4)
          sprint()

    #-------------------------#
    # Write ms-tor input file #
    #-------------------------#
    sprint("Writing mstor files...",NIBS2)
    rw.write_mstorinp(geoms,symbols,mtp,freqscal,inpvars,dMj,nrics,folder,lCH3)
    sprint()
#==================================================#

