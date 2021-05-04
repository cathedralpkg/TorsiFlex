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
| Sub-module :  tfrw               |
| Last Update:  2020/12/21 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different
functions for reading and writing
files in TorsiFlex
'''

#==================================================#
import os
import common.files              as ff
import common.fncs               as fncs
import common.physcons           as pc
import common.internal           as intl
#--------------------------------------------------#
from   common.Molecule   import Molecule
from   common.pgs        import get_pgs
from   common.Exceptions import WrongDimension
#--------------------------------------------------#
import modtorsiflex.tfhelper  as tfh
import modtorsiflex.tfvars    as tvars
#--------------------------------------------------#
from   modtorsiflex.tpespoint import TorPESpoint
#==================================================#


#==================================================#
def readfile_xyz(fname):
    string = ""
    zmat, symbols, masses = ff.read_xyz_zmat(fname)
    xcc = intl.zmat2xcc(zmat[0],zmat[1])
    # Print more information
    ndummy = symbols.count("XX")
    natoms = len(symbols)-ndummy
    # without dummies
    symbols_wo , xcc_wo    = fncs.clean_dummies(symbols,xcc=xcc)
    symbols_wo , masses_wo = fncs.clean_dummies(symbols,masses=masses)
    molecule = Molecule()
    molecule.setvar(xcc=xcc_wo,symbols=symbols_wo,masses=masses_wo,ch=0,mtp=1)
    molecule.prepare()
    molecule.setup()
    string += "Molecular formula     : %s\n"%molecule._mform
    string += "Number of atoms       : %i\n"%natoms
    string += "Number of dummy atoms : %i\n"%ndummy
    string += "Vibrational d.o.f.    : %i\n"%molecule._nvdof
    string += "\n"
    # Cartesian Coordinates
    string += "Cartesian coordinates (in Angstrom):\n"
    sidx = "%%%si"%len(str(len(symbols)))
    at_wo = -1
    dwithout = {}
    for at,symbol in enumerate(symbols):
        xi,yi,zi = [value*pc.ANGSTROM for value in xcc[3*at : 3*at+3]]
        mass     = masses[at]*pc.AMU
        datainline = (sidx%(at+1),symbol,xi,yi,zi,mass)
        if symbol != "XX": at_wo += 1; dwithout[at] = at_wo
        else: dwithout[at] = None
        string += "[%s]  %-2s  %+10.5f  %+10.5f  %+10.5f  (%7.3f amu)\n"%datainline
    string += "\n"
    return xcc,zmat,symbols,masses,(dwithout,molecule,natoms,ndummy),string
#==================================================#


#==================================================#
# Functions for reading/writing torsiflex files    #
#==================================================#
def data_from_zmat_lines(zmatlines,torsions={}):
    # Strip lines and remove commas and equal signs: "," --> " "; "=" --> " "
    zmatlines = [line.strip().replace(","," ").replace("="," ") for line in zmatlines]
    # Store correct lines
    zmatlines = [line for line in zmatlines if line != ""]
    # Convert selected torsions to dtorX
    for idx,line in enumerate(zmatlines):
        words = line.split()
        for X,ic in torsions.items():
            dtorX = tfh.X2tname(X)
            if ic in words:
                words[words.index(ic)] = dtorX
                zmatlines[idx] = " ".join(words)
    # Get symbols from lines
    allsymbols = [line.split()[0] for line in zmatlines if len(line.split()) != 2]
    # Put symbols in correct format
    allsymbols = fncs.correct_symbols(allsymbols)
    # Get zmat list and values of each coordinate
    zmat     = [line for line in zmatlines if len(line.split()) != 2]
    zmatvals = {line.split()[0]:float(line.split()[1]) for line in zmatlines if len(line.split()) == 2}
#   # Just in case, check dtorX
#   for line in zmat:
#       if "dtor" in line:
#           idx = line.find("dtor")
#           ttorsion = line[idx:].split()[0]
#           if ttorsion not in zmatvals.keys(): zmatvals[ttorsion] = 180.0
    # Return data
    return zmat, zmatvals, allsymbols
#--------------------------------------------------#
def read_zmat(filename,torsions={}):
    with open(filename,'r') as asdf: lines = asdf.readlines()
    # Read data
    zmat, zmatvals, symbols = data_from_zmat_lines(lines,torsions)
    # Return data
    return zmat, zmatvals, symbols
#--------------------------------------------------#
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
def read_llhlcorr(dirll,dirhl):
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
def write_llhlcorr(correlations):
    string = "".join("%s  -->  %s\n"%(vecLL,vecHL) for (vecLL,vecHL) in correlations.items())
    with open(tvars.LLHLCORRFILE,'w') as asdf: asdf.write(string)
#--------------------------------------------------#
def string4xyzfile(xcc,symbols,info=""):
    string  = "%i\n"%len(symbols)
    string += info
    if not info.endswith("\n"): string += "\n"
    for at,symbol in enumerate(symbols):
        x,y,z = fncs.xyz(xcc,at)
        x *= pc.ANGSTROM
        y *= pc.ANGSTROM
        z *= pc.ANGSTROM
        string += "%2s  %+9.5f  %+9.5f  %+9.5f\n"%(symbol,x,y,z)
    return string
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
           string += "%2s  %+9.5f  %+9.5f  %+9.5f\n"%(symbol,x,y,z)
       with open(xyzfile,'w') as asdf: asdf.write(string)
#--------------------------------------------------#
def write_molden_allconfs(dataconfs,Eref,allsymbols,folder):
       print(tvars.IBS2+"Generating xyz file with all geometries...")
       print(tvars.IBS2+"    filename: %s"%(folder+tvars.ALLCONFS))
       string = ""
       for idx,conftuple in enumerate(dataconfs):
           vec,V0,V1,G,weight,Qrv,ifreq,zmat,zmatvals,log = conftuple
           # convert to xyz
           V0_rel = (V0-Eref)*pc.KCALMOL
           description = " * E = %+7.3f kcal/mol ; (%i) %s"%(V0_rel,idx+1,str(vec))
           symbols   = [symbol for symbol in allsymbols if symbol.upper() not in ("X","XX")]
           natoms    = len(symbols)
           xcc = intl.zmat2xcc(zmat,zmatvals)
           # remove dummy atoms (if there is any)
           xcc = fncs.flatten_llist( [fncs.xyz(xcc,idx) for idx,symbol in enumerate(allsymbols) \
                                      if symbol.upper() not in ("X","XX")])
           # add to string
           string += string4xyzfile(xcc,symbols,info=description)
       with open(folder+tvars.ALLCONFS,'w') as asdf: asdf.write(string)
       print("")
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







