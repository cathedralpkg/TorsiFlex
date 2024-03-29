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
| Sub-module :  inpvars            |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the InpVars
class for TorsiFlex
'''

#==================================================#
import os
import getpass
import itertools
#--------------------------------------------------#
from  common.fncs       import float_in_domain
from  common.Exceptions import ErrorTorsion1
from  common.Exceptions import ErrorTorsionN
from  common.Exceptions import ErrorTorsionRepe
from  common.Exceptions import UnknownKey
#--------------------------------------------------#
import modtorsiflex.tfhelper  as tfh
import modtorsiflex.tfrw as rw
#==================================================#


#==================================================#
# A class to keep track of all the variables       #
#==================================================#
class InpVars():

      def __init__(self):
          # booleans #
          self._enantio      = False
          self._ts           = False
          self._fccards      = False
          # strings #
          self._zmatfile     = None
          self._keyll        = None
          self._keyhl        = None
          self._tmpll        = None
          self._tmphl        = None
          self._dirll        = "files_LL/"
          self._dirhl        = "files_HL/"
          self._pcfile       = None
          self._lowlevel     = "HF    3-21G              "
          self._highlevel    = "B3LYP 6-31G int=ultrafine"
          self._memll        = "1GB"
          self._memhl        = "1GB"
          # int/floats #
          self._nprocll      = 1
          self._nprochl      = 1
          self._ncycles      = 200
          self._optmode      = 1
          self._dist1d       = 7
          self._cfactor      = 1.3
          self._epsdeg       = 2
          self._charge       = 0
          self._multipl      = 1
          self._hlcutoff     = None
          self._freqscalLL   = 1.000
          self._freqscalHL   = 1.000
          self._sigmamj      = 0.02 # max value for sigma(mj)
          # lists
          self._hconstr      = []
          self._sconstr      = []
          self._ifqrangeLL   = []
          self._ifqrangeHL   = []
          self._temps        = []
          self._skipconn     = []
          self._stoc         = []
          # dictionaries
          self._tic          = {} # now dict, later list
          self._tsigma       = {} # now dict, later list
          self._precond      = {} # now dict, later list
          self._tdomain      = {} # now dict, later list

          # not in input file
          self._ttorsions    = None
          self._ntorsions    = 0
          self._tlimit       = None
          self._vecR         = None

          # to add externally
          self._tatoms       = {}
          self._prec         = (1,1)

          # TESTS
          self._tests        = {}
          self._tests[0]     = [1,1,1,1] # tests for guess     geom
          self._tests[1]     = [1,1,1,1] # tests for optimized geom (LL and HL)

      def __str__(self):
          try   : return self.string_vars()
          except: return ""

      def write_default(self,filename,PROGNAME,torsions=None,rewrite=False):
          '''
          Generates input if it does not exist
            * returns 0 if input exists
            * return  1 if input was created
          '''
          # Does input file exist?
          if os.path.exists(filename) and not rewrite: return 0
          if torsions is None: torsions = ["--"]
          # any zmat file?
          files = [fname for fname in os.listdir(".") if \
                   fname.lower().endswith(".zmat")]
          if len(files) == 1:
             zmatfile = files[0]
          else:
              zmatfile = None
              for fname in files:
                  if fname.endswith(".zmat"): zmatfile = fname
              if zmatfile is None: zmatfile = "zmfile.zmat"

          # Name for tmp folders
          mname    = ".".join(zmatfile.split(".")[:-1])
          try   : username = getpass.getuser()
          except: username = "user"
          tmpll    = "/scratch/%s/LL_%s/"%(username,mname)
          tmphl    = "/scratch/%s/HL_%s/"%(username,mname)
          while True:
              if not os.path.exists(tmpll): break
              if not os.path.exists(tmphl): break
              randname = tfh.random_name()
              tmpll    = "/scratch/%s/LL_%s/"%(username,randname)
              tmphl    = "/scratch/%s/HL_%s/"%(username,randname)

          if self._enantio: senantio = "yes"
          else            : senantio = "no "
          if self._ts     : sts      = "yes"
          else            : sts      = "no "
          if self._fccards: sffcards = "yes"
          else            : sffcards = "no "
          stestsG = " ".join("%s"%val for val in self._tests[0])
          stestsO = " ".join("%s"%val for val in self._tests[1])

          # Create input file
          string  = "# This is a %s input file\n"%(PROGNAME.split(".py")[0])
          string += "\n"
          string += "#-----------------------#\n"
          string += "#        System         #\n"
          string += "#-----------------------#\n"
          string += "zmatfile    %-11s # Z-matrix file\n"%zmatfile
          string += "charge      %-3i         # charge of the system\n"%self._charge
          string += "multipl     %-3i         # multiplicity of the system\n"%self._multipl
          string += "enantio     %3s         # yes if torsional enantioners, no otherwise\n"%senantio
          string += "ts          %3s         # yes if transition state, no otherwise\n"%sts
          string += "cfactor     %.1f         # controls the connectivity criterium\n"%self._cfactor
          string += "#skipcon    (1,2)       # skips connectivity between pairs of atoms\n"
          string += "\n"
          string += "\n"
          string += "#-----------------------#\n"
          string += "#        Storage        #\n"
          string += "#-----------------------#\n"
          string += "dirll       %s   # folder to store LL conformers\n"%self._dirll
          string += "dirhl       %s   # folder to store HL conformers\n"%self._dirhl
          string += "tmpll       %s # folder for LL temporal files\n"%tmpll
          string += "tmphl       %s # folder for HL temporal files\n"%tmphl
          string += "\n"
          string += "\n"
          string += "#-----------------------#\n"
          string += "#    Target torsions    #\n"
          string += "#-----------------------#\n"
          for idx,torsion in enumerate(torsions):
              ni = idx+1
              string += "torsion%i   %-12s # name of target torsion number %i in the Z-matrix file\n"%(ni,torsion,ni)
          string += "\n"
          for ni in range(1,len(torsions)+1):
              string += "precond%i   60 180 300   # precond angles for torsion%i\n"%(ni,ni)
          string += "\n"
          for ni in range(1,len(torsions)+1):
              string += "tdomain%i   (0,360)      # domain for torsion%i\n"%(ni,ni)
          string += "\n"
          for ni in range(1,len(torsions)+1):
              string += "tsigma%i    1            # symmetry number for torsion%i\n"%(ni,ni)
          string += "\n"
          string += "#pcfile     precond.txt # file with precond. angles\n"
          string += "\n"
          string += "\n"
          string += "#-----------------------#\n"
          string += "#   Search  Procedure   #\n"
          string += "#-----------------------#\n"
          string += "ncycles     %3i         # number of steps of stochastic algorithm\n"%self._ncycles
          string += "\n"
          string += "\n"
          string += "#-----------------------#\n"
          string += "#   HL reoptimization   #\n"
          string += "#-----------------------#\n"
          string += "#hlcutoff    5.0        # Gibbs energy cutoff (kcal/mol)\n"
          string += "\n"
          string += "\n"
          string += "#-----------------------#\n"
          string += "#   Validation  tests   #\n"
          string += "#-----------------------#\n"
          string += "testsG      %s     # for Guess geom (Conn, Simil, Hard, Soft)\n"%stestsG
          string += "testsO      %s     # for Opt   geom (Conn, Redun, Hard, Soft)\n"%stestsO
          string += "dist1D      %i           # domain size about each point (degrees)\n"%self._dist1d
          string += "epsdeg      %i           # max diff between two identical angles (degrees)\n"%self._epsdeg
          string += "#hconstr    ic domain   # hard constraint (see manual)\n"
          string += "#sconstr    ic domain   # soft constraint (see manual)\n"
          string += "#ifqrangeLL domain      # restricts LL imaginary-frequency interval\n"
          string += "#ifqrangeHL domain      # restricts HL imaginary-frequency interval\n"
          string += "\n"
          string += "\n"
          string += "#-----------------------#\n"
          string += "# Gaussian calculations #\n"
          string += "#-----------------------#\n"
          string += "optmode     %i           # 0:opt(z-matrix) , 1:opt(modredundant)\n"%self._optmode
          string += "fccards     %3s         # Use LL Hessian in HL opt (yes/no)\n"%sffcards
          string += "lowlevel    %s # low-level of calculation\n"%self._lowlevel
          string += "highlevel   %s # high-level of calculation\n"%self._highlevel
          string += "nprocll     %i           # Number of threads (low-level)\n"%self._nprocll
          string += "memll       %s         # dynamic memory (low-level)\n"%self._memll
          string += "nprochl     %i           # Number of threads (high-level)\n"%self._nprochl
          string += "memhl       %s         # dynamic memory (high-level)\n"%self._memhl
      #   string += "# Energy keyword for Gaussian ouput files\n"
      #   string += "#keyll      hf\n"
      #   string += "#keyhl      hf\n"
          string += "\n"
          string += "\n"
          string += "#-----------------------#\n"
          string += "#  Partition functions  #\n"
          string += "#-----------------------#\n"
          string += "temps     100.00 200.00 # temperatures (K) for part. functions\n"
          string += "temps     298.15 300.00 # temperatures (K) for part. functions\n"
          string += "temps     500.00 750.00 # temperatures (K) for part. functions\n"
          string += "temps           1000.00 # temperatures (K) for part. functions\n"
          string += "temps           2000.00 # temperatures (K) for part. functions\n"
          string += "temps           2500.00 # temperatures (K) for part. functions\n"
          string += "freqscalLL  %.3f       # freq. scaling factor (LL)\n"%self._freqscalLL
          string += "freqscalHL  %.3f       # freq. scaling factor (HL)\n"%self._freqscalHL
          string += "sigmamj     %.2f        # max value for sigma(Mj); >= 0.01\n"%self._sigmamj
          string += "\n"
          string += "\n"
          # write string to file
          with open(filename,'w') as asdf: asdf.write(string)
          return 1

      def read_input(self,filename):
          # Read lines from file
          with open(filename,'r') as asdf: lines = asdf.readlines()
          # Look for variables in lines
          for line in lines:
              line = line.split("#")[0].strip()
              if line == "": continue
              # split line into variable and value
              # and update variable
              try:
                variable = line.split()[0].lower()
                value    = " ".join(line.split()[1:]).strip()
                self.setvar(variable,value)
              except ErrorTorsionRepe as error: raise error
              except UnknownKey       as error: raise error
              except: continue
          # check torsions
          Xs = [int(X) for X in self._tic.keys()]
          if len(Xs) != 0:
             minX = min(Xs); error1 = ErrorTorsion1; error1._var = minX
             maxX = max(Xs); errorN = ErrorTorsionN; errorN._var = (len(Xs),maxX)
             if minX != 1      : raise error1
             if maxX != len(Xs): raise errorN

      def setvar(self,variable,value):
          # lower case, to avoid problems
          variable = variable.lower()

          #-----------------#
          # System-specific #
          #-----------------#
          if   variable == "zmatfile":
             self._zmatfile = value.strip()
          elif variable == "enantio":
             if value.lower() == "yes": self._enantio = True
          elif variable == "ts":
             if value.lower() == "yes": self._ts = True
          elif variable == "cfactor":
             self._cfactor   = float(value)
          elif variable == "skipcon":
              pairs = [pair.replace("(","").strip() for pair in value.split(")")]
              pairs = [pair for pair in pairs if "," in pair]
              pairs = [tuple(int(atom)-1 for atom in pair.split(",")) for pair in pairs]
              self._skipconn += pairs

          #-----------------#
          # Target torsions #
          #-----------------#
          elif variable.startswith("torsion"):
             X = variable.split("torsion")[1].strip()
             if X in self._tic.keys():
                errorRepe = ErrorTorsionRepe
                errorRepe._var = X
                raise errorRepe
             if value != "--": self._tic[X] = value
          elif variable.startswith("tdomain"):
             X = variable.split("tdomain")[1].strip()
             self._tdomain[X] = value
          elif variable.startswith("tsigma"):
             X = variable.split("tsigma")[1].strip()
             self._tsigma[X] = int(value)
          elif variable.startswith("precond"):
             X = variable.split("precond")[1].strip()
             value = value.replace(","," ")
             self._precond[X] = [int(round(float(val)))%360 for val in value.split()]
          elif variable == "pcfile":
             self._pcfile = value.strip()

          #------------------#
          # Search procedure #
          #------------------#
          elif variable == "ncycles":
             self._ncycles = int(value)

          #-------------------#
          # HL reoptimization #
          #-------------------#
          elif variable == "hlcutoff":
             self._hlcutoff = float(value)

          #------------------#
          # Validation tests #
          #------------------#
          elif variable == "testsg":
             value = [int(val_i) for val_i in value.replace(","," ").split()]
             if len(value) == 4 and value.count(0)+value.count(1) == 4:
                 self._tests[0] = value
          elif variable == "testso":
             value = [int(val_i) for val_i in value.replace(","," ").split()]
             if len(value) == 4 and value.count(0)+value.count(1) == 4:
                 self._tests[1] = value
          elif variable == "dist1d":
             self._dist1d = max(float(value),2.0)
          elif variable == "epsdeg":
             self._epsdeg  = min(float(value),self._epsdeg)
          elif variable == "hconstr":
             ic,icdomain = value.split()
             self._hconstr.append( (ic,icdomain) )
          elif variable == "sconstr":
             ic,icdomain = value.split()
             self._sconstr.append( (ic,icdomain) )
          elif variable == "ifqrangell":
             self._ifqrangeLL = tfh.str2ifreqinterval(value.strip())
          elif variable == "ifqrangehl":
             self._ifqrangeHL = tfh.str2ifreqinterval(value.strip())

          #----------------#
          # Gaussian calcs #
          #----------------#
          elif variable == "optmode":
             self._optmode = int(value)
          elif variable == "fccards":
             if value.lower() == "yes": self._fccards = True
          elif variable == "charge":
             self._charge = int(value)
          elif variable == "multipl":
             self._multipl = int(value)
          elif variable == "lowlevel":
             self._lowlevel = value.strip()
          elif variable == "highlevel":
             self._highlevel = value.strip()
          elif variable == "nproc":
             self._nprocll = int(value)
             self._nprochl = int(value)
          elif variable == "nprocll":
             self._nprocll = int(value)
          elif variable == "nprochl":
             self._nprochl = int(value)
          elif variable == "mem":
             self._memll   = value.strip()
             self._memhl   = value.strip()
          elif variable == "memll":
             self._memll   = value.strip()
          elif variable == "memhl":
             self._memhl   = value.strip()
        # elif variable == "keyll":
        #    self._keyll = value.strip()
        # elif variable == "keyhl":
        #    self._keyhl = value.strip()

          #-----------------#
          # Part. functions #
          #-----------------#
          elif variable == "temps":
             self._temps += [float(val_i) for val_i in value.split()]
          elif variable == "freqscalll":
             self._freqscalLL = float(value)
          elif variable == "freqscalhl":
             self._freqscalHL = float(value)
          elif variable == "sigmamj":
             self._sigmamj = max(float(value),0.01)

          #---------#
          # Storage #
          #---------#
          elif variable == "tmpll":
             self._tmpll = value.strip()
             if not self._tmpll.endswith("/"): self._tmpll += "/"
          elif variable == "tmphl":
             self._tmphl = value.strip()
             if not self._tmphl.endswith("/"): self._tmphl += "/"
          elif variable == "dirll":
             self._dirll = value.strip()
             if not self._dirll.endswith("/"): self._dirll += "/"
          elif variable == "dirhl":
             self._dirhl = value.strip()
             if not self._dirhl.endswith("/"): self._dirhl += "/"

          else:
              print("    UNKNOWN KEYWORD: %s"%variable)
              print("")
              raise UnknownKey

      def prepare_variables(self):
          self._ttorsions = list(self._tic.keys())
          # Sort them considering int(X)
          self._ttorsions.sort(key=int)
          # prepare temperatures for partition functions
          if self._temps == []: self._temps = list(range(100,2501,100))
          self._temps.sort()
          # number of torsions
          self._ntorsions = len(self._ttorsions)
          for X in self._ttorsions:
              # (a) torsions symmetry number
              if X not in self._tsigma.keys() : self._tsigma[X] = 1
              # (b) domain for each torsion
              if X not in self._tdomain.keys(): self._tdomain[X] = "(0,%i)"%(360/self._tsigma[X])
              self._tdomain[X] = tfh.str2interval(self._tdomain[X],bint=True)
              # (c) default angles
              default = self._precond.get(X,[60,180,300])
              self._precond[X] = [val for val in default if float_in_domain(val,self._tdomain[X])]
          # Convert to lists!
          self._tic     = [self._tic[X]        for X in self._ttorsions]
          self._tsigma  = [self._tsigma[X]     for X in self._ttorsions]
          self._tlimit  = [int(round(360/s))   for s in self._tsigma   ]
          self._precond = [self._precond[X]    for X in self._ttorsions]
          self._tdomain = [self._tdomain[X]    for X in self._ttorsions]
          # fix according to limit
          self._precond = [ [angle for angle in precond if angle <= tlimit] \
                           for precond,tlimit in zip(self._precond,self._tlimit)]
          self._tdomain = [[(min(phi1,tlimit),min(phi2,tlimit)) for phi1,phi2 in tdomain] \
                           for tdomain,tlimit in zip(self._tdomain,self._tlimit)]
          self._tdomain = [[(phi1,phi2) for phi1,phi2 in tdomain if phi1 != phi2] \
                           for tdomain,tlimit in zip(self._tdomain,self._tlimit)]
          # convert domains in constraints
          self._hconstr = [(ic,tfh.str2interval(icdomain)) for ic,icdomain in self._hconstr]
          self._sconstr = [(ic,tfh.str2interval(icdomain)) for ic,icdomain in self._sconstr]


      def check(self):
          if self._zmatfile is None       : return -1, "zmatfile"
          if self._tmpll    is None       : return -1, "tmpll"
          if self._tmphl    is None       : return -1, "tmphl"
          if self._dirll    is None       : return -1, "dirll"
          if self._dirhl    is None       : return -1, "dirhl"
          if self._optmode  not in [0,1]  : return -1, "optmode"
          return 0, ""

      def yield_precond(self):
          # precond angles from file (and deftor if needed)
          if self._pcfile is not None and os.path.exists(self._pcfile):
             angles = []
             for values in rw.read_pcfile(self._pcfile):
                 # initialize vector
                 vec = [None for torsion in self._ttorsions]
                 # Add defined values in line
                 for idx,torsion in enumerate(self._ttorsions):
                     vec[idx] = values.get("torsion"+torsion,None)
                 # Vector is complete?
                 if None not in vec: angles.append( vec )
                 else:
                    nodefined = [idx for idx,vi in enumerate(vec) if vi is None]
                    # predefined values for those angles
                    lpredef = [self._precond[idx] for idx in nodefined]
                    # generate combinations using deftors
                    for subvec in itertools.product(*lpredef):
                        fvec = list(vec)
                        # update fvec
                        for idx,value in zip(nodefined,subvec): fvec[idx] = value
                        # add fvec to angles
                        angles.append( fvec )
          # precond angles using default values
          else:
             lpredef = tuple([precond for precond in self._precond])
             angles  = [vec for vec in itertools.product(*lpredef)]
          # Get closest to reference vector
          try:
            vecR = [int(ii) for ii in str(self._vecR).split("_")]
            midx,mindist = None,float("inf")
            for idx,vec in enumerate(angles):
                dist = 0.0
                for aa,bb in zip(vecR,vec):
                    dist += min((aa-bb)%360,(bb-aa)%360)
                    if dist > mindist: break
                if dist < mindist: midx,mindist = idx,dist
            if midx is not None: angles = angles[midx:] + angles[0:midx]
          except: pass
          #------------#
          # yield data #
          #------------#
          nblocks = self._prec[0]  # number of blocks
          bidx    = self._prec[1]  # block idx
          blocks  = {idx:[] for idx in range(1,nblocks+1)}
          # divide data into blocks
          count = 1
          for point in angles:
              blocks[count].append(point)
              count += 1
              if count > nblocks: count = 1
          # yield vectors in block
          for block in blocks[bidx]: yield block
#==================================================#

