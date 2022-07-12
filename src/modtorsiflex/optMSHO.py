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
| Sub-module :  optMSHO            |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''

#==================================================#
import os
import numpy  as     np
from   shutil import move
#--------------------------------------------------#
import common.fncs       as fncs
import common.physcons   as pc
import common.Exceptions as exc
from   common.files      import mkdir_recursive
from   common.internal   import zmat2xcc
from   common.Molecule   import Molecule
#--------------------------------------------------#
import modtorsiflex.tfgau     as tgau
import modtorsiflex.tfrw      as rw
import modtorsiflex.tfhelper  as tfh
import modtorsiflex.printing  as pp
from   modtorsiflex.tfvars    import NIBS2, IBS2
from   modtorsiflex.tfvars    import EPS_KCALMOL, DIREXCLUD
from   modtorsiflex.tfvars    import MINGTXT, ENERGYSUMLL, ALLCONFS
#==================================================#


REASONS  = "wrongconn,repeated,outofdomain,norestr,wimag"


#==================================================#
def read_dipole(log):
    with open(log,'r') as asdf: lines = asdf.readlines()
    text = "".join([line.strip() for line in lines])
    text = text.split("\\Dipole=")[1].split("\\DipoleDeri")[0]
    mux,muy,muz = text.split(",")
    aa = np.array( [mux,muy,muz] )
    # in debyes
    dipole = np.linalg.norm(aa)*pc.DEBYE
    return dipole
#--------------------------------------------------#
def get_rotcons(lzmat,zmatvals,symbols):
    xcc = zmat2xcc(lzmat,zmatvals)
    # correct symbols (just in case)
    symbols,atonums = fncs.symbols_and_atonums(symbols)
    # Data without dummies
    symbols_wo,xcc_wo = fncs.clean_dummies(symbols,xcc=list(xcc))
    # generate Molecule instance
    molecule = Molecule()
    molecule.setvar(xcc=xcc_wo)
    molecule.setvar(symbols=symbols_wo)
    molecule.setvar(V0=0)
    molecule.setvar(ch=0,mtp=1)
    molecule.prepare()
    molecule.setup()
    # in freq units (GHz)
    imoms   = [Ii * pc.KG*pc.METER**2 for Ii in molecule._imoms]
    rotcons = [pc.H_SI/(Ii*8*np.pi**2)/1E9 for Ii in imoms]
    return rotcons
#--------------------------------------------------#
def check_if_same_energy(vec,gts,V0,V1,idx,data,inpvars,EXCLUDED,dcorr):
    #-------------------------------------------#
    # Select those conformers with same V0 & V1 #
    #-------------------------------------------#
    idx_same = []
    for idx_j in range(idx-1,-1,-1):
        vec_j,V0j,V1j,weightj,rotconsj,ifreqj,Qrvj,Gj,lzmatj,zmatvalsj,sfccardsj,gtsj = data[idx_j]
        diffV0 = abs(V0j-V0)*pc.KCALMOL
        diffV1 = abs(V1j-V1)*pc.KCALMOL
        if diffV0 > EPS_KCALMOL: break
        if diffV1 > EPS_KCALMOL: break
        idx_same.append(idx_j)

    #----------------#
    # Compare angles #
    #----------------#
    same,enantio,VECj,GTSj = False,False,None,None
    for idx_j in idx_same:
        repeated = False
        # Data of previous points
        vec_j,V0j,V1j,weightj,rotconsj,ifreqj,Qrvj,Gj,lzmatj,zmatvalsj,sfccardsj,gtsj = data[idx_j]
        # check if same point
        if vec.is_same(vec_j,inpvars._epsdeg):
           same,VECj,GTSj = True,vec_j,gtsj
           break
        # System presents torsional enantiomerism?
        if not inpvars._enantio: continue
        # Try with enantiomer vector
        evec_j = tfh.enantiovec(lzmatj,zmatvalsj,dcorr,inpvars)
        if vec.is_same(evec_j,inpvars._epsdeg):
           enantio,VECj,GTSj = True,vec_j,gtsj
           break

    # Save if repeated
    if same or enantio:
       # if possible, keep the preconditiones one
       if gts.startswith("prec.") and GTSj.startswith("stoc."): EXCLUDED["repeated"].append(GTSj)
       else                                                   : EXCLUDED["repeated"].append(gts)
    return EXCLUDED,(same,enantio,VECj)
#--------------------------------------------------#
def check_point(gts,lzmat,zmatvals,cmatrix,inpvars,ifreq,ifreqconstr,EXCLUDED):

    #----------------#
    # Carry on tests #
    #----------------#
    bools = tfh.precheck_geom(lzmat,zmatvals,cmatrix,inpvars)[0]
    bool_imag = True
    if ifreqconstr != []:
       bool_imag = fncs.float_in_domain(fncs.afreq2cm(abs(ifreq)),ifreqconstr)

    #-----------------------#
    # If excluded, classify #
    #-----------------------#
    allok = False
    if   (inpvars._tests[1][0] == 1) and (bools[0] is False): EXCLUDED["wrongconn"  ].append(gts)
    elif (inpvars._tests[1][1] == 1) and (bools[1] is False): EXCLUDED["outofdomain"].append(gts)
    elif (inpvars._tests[1][2] == 1) and (bools[2] is False): EXCLUDED["norestr"    ].append(gts)
    elif (inpvars._tests[1][3] == 1) and (bools[3] is False): EXCLUDED["norestr"    ].append(gts)
    elif not bool_imag                                      : EXCLUDED["wimag"      ].append(gts)
    else                                                    : allok = True

    # Return dictionary and boolean
    return EXCLUDED, allok
#--------------------------------------------------#
def cleanup_conformers(data,inpvars,cmatrix,symbols,ifreqconstr,dcorr):

    # Analyze point by point
    EXCLUDED = {reason:[] for reason in REASONS.split(",")}

    count_pre = 0
    count_sto = 0

    # Now, one by one
    extra = {}
    for idx,conftuple in enumerate(data):

        vec,V0,V1,weight,rotcons,ifreq,Qrv,G,lzmat,zmatvals,fccards,gts = conftuple
        # extra data
        try   : rotcons = sorted(list(get_rotcons(lzmat,zmatvals,symbols)))
        except: rotcons = None
        # preconditioned or stochastic?
        if gts.startswith("stoc"): count_sto  += weight
        else                     : count_pre  += weight
        # Check if it must be removed
        EXCLUDED,allok = check_point(gts,lzmat,zmatvals,cmatrix,inpvars,ifreq,ifreqconstr,EXCLUDED)
        # If all ok, then check if there are previous conformers with same energy
        if allok:
           EXCLUDED,extra_j = check_if_same_energy(vec,gts,V0,V1,idx,data,inpvars,EXCLUDED,dcorr)
           extra[idx]       = extra_j

    return EXCLUDED,extra, count_pre, count_sto
#--------------------------------------------------#
def relocate_excluded(folder,EXCLUDED):
    nreloc  = 0
    direxcl = folder+DIREXCLUD

    n1,n2,n3,n4,n5 = [len(EXCLUDED[which]) for which in REASONS.split(",")]
    if n1+n2+n3+n4+n5 > 0:
       string  = "Conformers with problems (<--) will be moved to: %s\n"%direxcl
       mkdir_recursive(direxcl)
       for gts in [gts for which in REASONS.split(",") for gts in EXCLUDED[which]]:
           # gts file
           if not os.path.exists(folder+gts): continue
           move(folder+gts,direxcl+gts)
           nreloc += 1
           # zmat & molden files
           for ext in (".zmat",".molden"):
               filename = gts.replace(".gts",ext)
               if not os.path.exists(folder+filename): continue
               move(folder+filename,direxcl+filename)
       string += "Number of conformers moved: %i\n"%nreloc
       string += ""
       # repeat again all
       for line in string.split("\n"): print(IBS2+line)
    return nreloc
#--------------------------------------------------#
def classify_files(inpvars,cmatrix,dcorr,level="LL"):
    if level == "LL":
       folder      = inpvars._dirll
       fscal       = inpvars._freqscalLL
       ifreqconstr = inpvars._ifqrangeLL
    if level == "HL":
       folder      = inpvars._dirhl
       fscal       = inpvars._freqscalHL
       ifreqconstr = inpvars._ifqrangeHL
    # Folder exists?
    if not os.path.exists(folder):
       pp.print_dirnotfound(folder,NIBS2,1)
       return
   
    # print status of validation tests
    pp.print_tests(inpvars._tests,1)
    pp.print_ifreqconstr(ifreqconstr)
    pp.print_torsiangleclass()


    #------------------------------------#
    # Read data in folder and print info #
    #------------------------------------#
    nreloc = 1 # number of relocated files
    ntries = 0 # number of attempts
    while nreloc != 0:
       if ntries != 0: print("\n"+IBS2+"TRYING AGAIN...\n\n")
       # Read logs
       print(IBS2+"Reading files inside '%s'..."%folder)
       data, symbols = rw.folder_data(folder,inpvars,fscal)
       print(IBS2+"Number of gts files: %i\n"%len(data))
       if len(data) == 0: return

       if symbols is None: 
          print(IBS2+"ERROR! No data inside %s\n"%folder)
          raise exc.END

       # Get reference point
       data.sort(key=lambda x:x[1])
       vec0    = str(data[0][0])
       Eref    = data[0][1]
       # Clean-up conformers
       EXCLUDED,extra,npre,nsto = cleanup_conformers(data,inpvars,cmatrix,symbols,ifreqconstr,dcorr)
       pp.print_table_conformers(data,inpvars,fscal,EXCLUDED,extra,npre,nsto)
       nreloc = relocate_excluded(folder,EXCLUDED)
       # Just in case, break this loop if carried out 4 times
       ntries += 1
       if ntries == 4: break

    #----------------------#
    # Rotational constants #
    #----------------------#
    pp.print_table_rotcons(data)

    #------------------------------------------#
    # Partition functions for each temperature #
    #------------------------------------------#
    pp.print_partition_function(data,inpvars._temps,fscal,inpvars._enantio)

    #--------------------#
    # Create molden file #
    #--------------------#
    pp.sprint("Generating xyz file with all geometries: %s"%(folder+ALLCONFS),NIBS2)
    rw.write_molden_allconfs(data,folder)
    pp.sprint("Done!",NIBS2,1)
#==================================================#




