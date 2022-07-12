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
| Sub-module :  tfrw               |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different
functions for reading and writing
files in TorsiFlex
'''

#==================================================#
import os
import numpy as np
#--------------------------------------------------#
import common.files              as ff
import common.fncs               as fncs
import common.physcons           as pc
import common.internal           as intl
from   common.gaussian           import get_fccards_string
#--------------------------------------------------#
from   common.Molecule   import Molecule
from   common.pgs        import get_pgs
from   common.Exceptions import WrongDimension
#--------------------------------------------------#
import modtorsiflex.tfvars    as tvars
#--------------------------------------------------#
from   modtorsiflex.tpespoint import TorPESpoint
#==================================================#


#==================================================#
# Convert imoms to rotcons (in GHz)
IMOM2GHZ = pc.H_SI/(8*np.pi**2 * pc.KG*pc.METER**2)/1E9
#==================================================#

#==================================================#
# Functions for reading/writing torsiflex files    #
#==================================================#
#--------------------------------------------------#
def initialize_domains(ntorsions):
    if os.path.exists(tvars.FDOMAINS): return
    # info about statusFRQ
    info   = "# if statusFRQ <  0 --> frequency calculation was not required or failed\n"
    info  += "# if statusFRQ >= 0 --> number of imaginary frequencies\n"
    # table columns
    length = max(4*ntorsions-1,14)
    cols   = ["domain","%%-%is"%length%"initial_vector","%%-%is"%length%"final_vector","statusFRQ"]
    cols   = "#"+"  ".join([" %s "%col for col in cols])+"\n"
    with open(tvars.FDOMAINS,'w') as asdf: asdf.write(info+"\n"+cols)
#--------------------------------------------------#
def get_domainline(domain,ipoint,fpoint,statusFRQ):
    ll = max(max(len(str(ipoint)),len(str(fpoint))),14)
    line = "  %%-6s    %%-%is    %%-%is        %%2i\n"%(ll,ll)
    if ipoint is None: ipoint = "-"
    if fpoint is None: fpoint = "-"
    return line%(domain,str(ipoint),str(fpoint),statusFRQ)
#--------------------------------------------------#
def add_domain(domain,ipoint,fpoint,statusFRQ):
    # Get line to append to file
    line = get_domainline(domain,ipoint,fpoint,statusFRQ)
    # Write at the end of file
    with open(tvars.FDOMAINS,'a') as asdf: asdf.write(line)
#--------------------------------------------------#
def write_domains(ntorsions,ddomains):
    if os.path.exists(tvars.FDOMAINS): return
    initialize_domains(ntorsions)
    for domain in ddomains.keys():
        for ipoint,fpoint,statusFRQ in ddomains[domain]:
            add_domain(domain,ipoint,fpoint,statusFRQ)
#--------------------------------------------------#
def read_llhlcorr():
    correlations = {}
    if not os.path.exists(tvars.LLHLCORRFILE): return correlations
    with open(tvars.LLHLCORRFILE,'r') as asdf: lines = asdf.readlines()
    for line in lines:
        line = line.split("#")[0].strip()
        if line == "": continue
        # remove arrow
        line = line.replace("-->"," ")
        # get vectors
        vecLL,vecHL = line.split()
        # update dict
        correlations[vecLL] = vecHL
    return correlations
#--------------------------------------------------#
def add_llhlcorr(vecLL,vecHL):
    line = "%s  -->  %s\n"%(vecLL,vecHL)
    with open(tvars.LLHLCORRFILE,'a') as asdf: asdf.write(line)
#--------------------------------------------------#
def write_llhlcorr(correlations):
    string = "".join("%s  -->  %s\n"%(vecLL,vecHL) for (vecLL,vecHL) in correlations.items())
    with open(tvars.LLHLCORRFILE,'w') as asdf: asdf.write(string)
#--------------------------------------------------#
def create_xyz(xcc,symbols,inpvars):
    xyzfile = inpvars._zmatfile+".xyz"
    print(tvars.IBS+"   - xyz file: %s"%xyzfile)
    if not os.path.exists(xyzfile):
       string = "%i\nxyz file from %s\n"%(len(symbols),inpvars._zmatfile)
       for at,symbol in enumerate(symbols):
           x,y,z = fncs.xyz(xcc,at)
           x *= pc.ANGSTROM
           y *= pc.ANGSTROM
           z *= pc.ANGSTROM
           string += "%2s  %+12.8f  %+12.8f  %+12.8f\n"%(symbol,x,y,z)
       with open(xyzfile,'w') as asdf: asdf.write(string)
#--------------------------------------------------#
def gts2molden(gtsfile,moldenfile):
    molecule = Molecule()
    molecule.set_from_gts(gtsfile)
    molecule.setup()
    molecule.genfile_molden(moldenfile)
#--------------------------------------------------#
def write_molden_allconfs(data,folder):
    string = ""
    # reference energy
    minV0 = data[0][1]
    # loop
    for idx,data_i in enumerate(data):
        vec,V0,V1,weight,rotcons,ifreq,Qrv,G,lzmat,zmatvals,fccards,gts = data_i
        # convert to xyz
        V0_rel  = (V0-minV0)*pc.KCALMOL
        info    = " * E = %+7.3f kcal/mol ; (%i) %s"%(V0_rel,idx+1,str(vec))
        # Symbols and xcc from z-matrix
        allsymbols = [symbol for symbol,conns,keys in lzmat]
        allxcc = intl.zmat2xcc(lzmat,zmatvals)
        # Remove dummies
        symbols = [symbol for symbol in allsymbols if symbol.upper() not in ("X","XX")]
        xcc     = [fncs.xyz(allxcc,idx) for idx,symbol in enumerate(allsymbols) \
                                   if symbol.upper() not in ("X","XX")]
        xcc     = fncs.flatten_llist(xcc)
        natoms  = len(symbols)
        # add to string
        string += "%i\n"%natoms
        string += info
        if not info.endswith("\n"): string += "\n"
        for at,symbol in enumerate(symbols):
            x,y,z = fncs.xyz(xcc,at)
            x *= pc.ANGSTROM
            y *= pc.ANGSTROM
            z *= pc.ANGSTROM
            string += "%2s  %+12.8f  %+12.8f  %+12.8f\n"%(symbol,x,y,z)
    # Write file with all the conformers
    with open(folder+tvars.ALLCONFS,'w') as asdf: asdf.write(string)
#--------------------------------------------------#
def write_mstorinp(geoms,symbols,mtp,freqscal,inpvars,dMj,nrics,folder,lCH3=[]):
    str_inp  = ""
    str_inp += "#--------------------------------------------#\n"
    str_inp += "# Ms-Tor input file generated with TorsiFlex #\n"
    str_inp += "#--------------------------------------------#\n"
    str_inp += "\n"

    torsions  = inpvars._tatoms
    ntorsions = inpvars._ntorsions
    Lenantio  = inpvars._enantio
    natoms  = len(symbols)
    atonums = fncs.symbols2atonums(symbols)
    masses  = fncs.symbols2masses(symbols)
    Eref    = geoms[0][0]
    #-----------------#
    # general section #
    #-----------------#
    str_inp += "$general\n"
    str_inp += "  natoms    %i\n"%natoms
    str_inp += "  nstr      %i\n"%len(geoms)
    str_inp += "  ntor      %i\n"%(ntorsions+len(lCH3))
    str_inp += "  freqscale %.4f\n"%freqscal
    str_inp += "  deltat    0.5\n"
    str_inp += "  elec\n"
    str_inp += "    %i  0.0\n"%mtp
    str_inp += "  end\n"
    str_inp += "end\n"
    str_inp += "\n"

    #---------------------#
    # temperature section #
    #---------------------#
    str_inp += "$temp\n"
    for idx in range(0,len(inpvars._temps),6):
        str_inp += "  ".join("%7.2f"%T for T in inpvars._temps[idx:idx+6])+"\n"
    str_inp += "end\n"
    str_inp += "\n"

    #----------------------#
    # internal coordinates #
    #----------------------#
    str_inp += "$intdef\n"
    if nrics is not None:
        for ic in nrics:
            ictype,icatoms = ic
            # which torsion?
            dtorX = ""
            if ictype == "pt":
               for X,tatoms in inpvars._tatoms.items():
                   if icatoms in (tuple(tatoms),tuple(tatoms[::-1])):
                      dtorX = "# torsion%s"%X
                      break
               for CH3 in lCH3:
                   if icatoms in (tuple(CH3),tuple(CH3[::-1])):
                      dtorX = "# CH3"
                      break
            # add to string
            str_inp += "  %s %s\n"%(intl.ic2string(ic).replace("_","-"),dtorX)
    str_inp += "end\n"
    str_inp += "\n"

    #---------------------#
    # coordinates&hessian #
    #---------------------#
    count = 0
    str_hess = ""
    for E,xcc,gcc,Fcc,point in geoms:
        count += 1
        # relative energy
        erel = (E-Eref)*pc.KCALMOL
        # rotational symmetry number
        pgroup, rotsigma = get_pgs(atonums,masses,xcc)
        if Lenantio and pgroup.lower() == "c1": weight = 2
        else                                  : weight = 1
        # The string
        str_inp += "$structure  %i\n"%count
        str_inp += "geom\n"
        for idx,at in enumerate(atonums):
            x,y,z = [coord*pc.ANGSTROM for coord in xcc[3*idx:3*idx+3]]
            str_inp += "%-2s  %15.8E  %15.8E  %15.8E\n"%(at,x,y,z)
        str_inp += "end\n"
        # M_j,tau values
        str_inp += "mtor\n"
        part_ttorsions = " ".join(["%.4f"%dMj[str(point)] for torsion in inpvars._ttorsions])
        part_ch3       = " ".join("3.0000" for CH3 in lCH3)
        str_inp += "  %s %s \n"%(part_ttorsions,part_ch3)
        str_inp += "end\n"
        # weight, rel energy and rotsigma
        str_inp += "weight    %i\n"%weight
        str_inp += "energy    %.5f\n"%erel
        str_inp += "rotsigma  %i # point group: %s\n"%(rotsigma,pgroup)
        str_inp += "end\n"
        str_inp += "\n"
        # hessian
        str_hess += "$hess  %2i\n"%count
        for idx in range(0,len(Fcc),6):
            values = "".join([" %+11.8f"%value for value in Fcc[idx:idx+6]])
            str_hess += values + "\n"
        str_hess += "end\n\n"

    # write files
    print(tvars.IBS2+"    - %s"%(folder+tvars.MSTORIF1))
    with open(folder+tvars.MSTORIF1,"w") as asdf: asdf.write(str_inp)
    print(tvars.IBS2+"    - %s"%(folder+tvars.MSTORIF2))
    with open(folder+tvars.MSTORIF2,"w") as asdf: asdf.write(str_hess)
#==================================================#


#==================================================#
def folder_data(folder,inpvars,fscal=1.0,verbose=True):
    gtss = [fi for fi in os.listdir(folder) if fi.endswith(".gts")]
    # First preconditioned, then stochastic!
    gtss.sort()
    # Mode for calculation of partition functions
    if inpvars._ts: fmode = -1
    else          : fmode =  0
    # Calculate partition functions at different temperatures
    data    = []
    symbols = None
    for gts in gtss:
        name       = gts.split(".")[-2]
        zmatfile   = gts.replace(".gts",".zmat")
        moldenfile = gts.replace(".gts",".molden")
        point      = TorPESpoint(name,inpvars._tlimit)
        # Read zmat file
        (lzmat,zmatvals,zmatatoms), symbols, masses = ff.read_zmat(folder+zmatfile)
        # Molecule instance from gts
        molecule = Molecule()
        molecule.set_from_gts(folder+gts)
        molecule.setvar(fscal=fscal)
        molecule.prepare()
        V0 = molecule._V0

        # String fccards for gaussian
        gcc = list(molecule._gcc)
        Fcc = fncs.matrix2lowt(molecule._Fcc)
        sfccards = get_fccards_string(gcc,Fcc)

        # Initialize variables
        weight   = None
        rotcons  = None
        ifreq    = None
        Qrv      = None
        gibbs    = None
        V1       = None

        # Calculate previous variables if asked for
        if verbose:
           molecule.setup()
           if not os.path.exists(folder+moldenfile): molecule.genfile_molden(folder+moldenfile)
           molecule.ana_freqs()
           # Rotational constants (GHz)
           try   : rotcons = IMOM2GHZ / np.array(molecule._imoms)
           except: rotcons = None
           # Calculate partition functions for the temperatures
           qtot, V1, qis = molecule.calc_pfns(inpvars._temps,fmode=fmode)
           qtr,qrot,qvib,qele = qis
           Qrv = np.array(qrot) * np.array(qvib)
           # remove rotsigma for Gibbs
           gibbs = V1-pc.KB*np.array(inpvars._temps)*np.log(qtot*molecule._rotsigma)
           # Weight of conformer
           if inpvars._enantio and molecule._pgroup.lower() == "c1": weight = 2
           else                                                    : weight = 1
           # imaginary frequency
           ccimag = molecule._ccimag
           if   len(ccimag) == 0: ifreq = None
           elif len(ccimag) == 1: ifreq = [ifreq for ifreq,ivec in list(molecule._ccimag)][0]
           else                 : ifreq = None
           V1 = molecule._ccV1
        # data to save
        data_i = (point,V0,V1,weight,rotcons,ifreq,Qrv,gibbs,lzmat,zmatvals,sfccards,gts)
        data.append(data_i)
    # Return data
    return data,symbols
#==================================================#



def read_pcfile(pcfile):
    # read file
    with open(pcfile,'r') as asdf: lines = asdf.readlines()
    # Clean lines
    lines = [line.split("#")[0].strip() for line in lines]
    lines = [line for line in lines if line != ""]
    # select torsions
    sel_torsions = lines[0].replace("_"," ").lower().split()
    # add name of torsion, just in case
    sel_torsions = [t if t.startswith("torsion") else "torsion"+t for t in sel_torsions]
    # for each line, get vector
    for line in lines[1:]:
        # replace "_" by " "
        values = line.replace("_"," ").split()
        # check dimension
        if len(values) != len(sel_torsions):
           exception = WrongDimension(Exception)
           exception._ntor = len(sel_torsions)
           exception._fline = lines[0]
           exception._line = line
           raise exception
        # get values
        values = {tor:int(float(val)%360) for tor,val in zip(sel_torsions,values) \
                  if val !="-"}
        yield values
#--------------------------------------------------#
def read_domains():
    # Initialize domains
    ddomains  = {domain:set([]) for domain in tvars.DOMAINS}
    # does file exist? read it!
    if os.path.exists(tvars.FDOMAINS):
       with open(tvars.FDOMAINS,'r') as asdf:
          lines = asdf.readlines()
    else: lines = []
    # Analyze lines
    for line in lines:
        # ignore comments and empty lines
        line = line.split("#")[0].strip()
        if line == "": continue
        # Get data in line
        domain, ipoint, fpoint, statusFRQ = line.split()
        domain = domain.lower()
        # save to dictionary
        if ipoint == "-": ipoint = None
        else            : ipoint = TorPESpoint(ipoint)
        if fpoint == "-": fpoint = None
        else            : fpoint = TorPESpoint(fpoint)
        ddomains[domain].add( (ipoint,fpoint,int(statusFRQ)) )
    # Return domains
    return ddomains




#==================================================#
# 01. Adjacency matrix                             #
#==================================================#
def read_adjacency(fugdict):
    ''' reads file with connectivities '''
    ugdict = {}
    # read lines in file
    with open(fugdict,'r') as asdf: lines =asdf.readlines()
    # extract data
    for line in lines:
        line = line.split("#")[0].strip()
        if line == "": continue
        node, neighbors = line.split(":")
        neighbors = neighbors.split()
        ugdict[int(node)-1] = set([int(ni)-1 for ni in neighbors])
    # return dictionary
    return ugdict
#--------------------------------------------------#
def write_adjacency(fugdict,ugdict):
    ''' writes file with connectivities '''
    string = ""
    # some operations for nice formatting
    nn = fncs.num_of_digits(len(ugdict.keys()))
    # add data to string
    for node,neighbors in ugdict.items():
        string += "%%-%ii : "%nn%(node+1)
        string += " ".join(["%%%ii"%nn%(ni+1) for ni in neighbors])
        string += "\n"
    # write file
    with open(fugdict,'w') as asdf: asdf.write(string)
#==================================================#


def read_llenergies():
    with open(tvars.ENERGYSUMLL,'r') as asdf: lines = asdf.readlines()
    data = []
    for line in lines:
        line = line.split("#")[0].strip()
        if line == "": continue
        # temperature
        if "tempGibbs" in line:
           T = float(line.split()[1])
           continue
        #data
        line = line.replace("|"," ")
        idx, V0, G,svec = line.split()
        idx = int(idx)
        V0  = float(V0)
        G   = float(G)
        data.append( (idx,V0,G,svec) )
    data.sort()
    return data,T
#--------------------------------------------------#
def write_llenergies(dataconfs,temp):
    ni = max(len(str(dataconfs[0][0])),len("point"))
    string  = "# Gibbs free energy temperature\n"
    string += "tempGibbs %9.3f K\n"%temp
    string += "\n"
    string += "              #--------------------------#\n"
    string += "              #        kcal/mol          #\n"
    string += "#----------------------------------------#\n"
    string += "#  conformer  |     V0     |    Gibbs    #\n"
    string += "#----------------------------------------#\n"
    minV0 = min([conftuple_i[1] for conftuple_i in dataconfs])
    minG  = min([conftuple_i[3]  for conftuple_i in dataconfs])
    for idx_i,conftuple_i in enumerate(dataconfs):
        vec_i,V0_i,V1_i,G_i = conftuple_i[0:4]
        relV0_i = (V0_i-minV0)*pc.KCALMOL
        relG_i  = (G_i -minG )*pc.KCALMOL
        string += "   %-9i  |  %8.3f  |  %9.3f  | %s \n"%(idx_i+1,relV0_i,relG_i,vec_i)
    string += "#----------------------------------------#\n"
    with open(tvars.ENERGYSUMLL,'w') as asdf: asdf.write(string)

