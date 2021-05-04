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
| Sub-module :  optMSHO            |
| Last Update:  2020/12/21 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''

#==================================================#
import os
import numpy  as     np
from   shutil import move
#--------------------------------------------------#
import common.fncs     as fncs
import common.physcons as pc
import common.Exceptions as exc
from   common.files    import mkdir_recursive
from   common.internal   import zmat2xcc
from   common.Molecule   import Molecule
#--------------------------------------------------#
import modtorsiflex.tfgau   as tgau
import modtorsiflex.tfrw    as rw
import modtorsiflex.tfhelper     as tfh
import modtorsiflex.printing   as pp
from   modtorsiflex.tfvars       import NIBS2, IBS2, EPS_KCALMOL, DIREXCLUD, MINGTXT
#==================================================#




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
def classify_files(inpvars,cmatrix,case,dcorr):
    if case == "LL":
       folder = inpvars._dirll
       fscal  = inpvars._freqscalLL
    if case == "HL":
       folder = inpvars._dirhl
       fscal  = inpvars._freqscalHL
    # Folder exists?
    if not os.path.exists(folder):
       pp.print_dirnotfound(folder,NIBS2,1)
       return
   
    # print status of validation tests
    pp.print_tests(inpvars._tests,1)

    # Read logs
    print(IBS2+"Reading log files inside '%s'..."%folder)
    # has to be done in this way due to the yield in folder_data function
    data      = [(ctuple,symbols) for ctuple,symbols in tgau.folder_data(folder,inpvars,fscal)]
    print(IBS2+"Number of (frq) log files: %i\n"%len(data))
    if len(data) == 0: return
    symbols   = data[0][1]
    dataconfs = [ctuple  for ctuple,symbols in data]
    dataconfs.sort(key=lambda x:x[1])
    del data
    minG  = min([ctuple[3] for ctuple in dataconfs])
    if case == "LL":
       with open(MINGTXT,'w') as asdf: asdf.write("minG %.8f %.2f\n"%(minG,inpvars._temp))

    if symbols is None: 
       print(IBS2+"ERROR! No data inside %s\n"%folder)
       raise exc.END

    # Get reference point
    vec0    = str(dataconfs[0][0])
    Eref    = dataconfs[0][1]
    print(IBS2+"Minimum total energy:")
    print(IBS2+"  * point: %s"%str(vec0))
    print(IBS2+"  * value: %.7f hartree"%Eref)
    print("")
    print(IBS2+"Relative energies (in kcal/mol)")
    print(IBS2+"     E     : total energy")
    print(IBS2+"     E+ZPE : total energy plus zero-point energy")
    print(IBS2+"     G     : Gibbs free energy at %.2f K"%inpvars._temp)
    print(IBS2+"     X     : Contribution to the total partition function at %.2f K"%inpvars._temp)
    print("")
    print(IBS2+"     frequencies scaled by: %.4f"%fscal)
    print("")
    # get minimum energies
    minV0 = min([conftuple[1] for conftuple in dataconfs])
    minV1 = min([conftuple[2] for conftuple in dataconfs])
    minG  = min([conftuple[3] for conftuple in dataconfs])

    # Analyze point by point
    wrongconn   = []
    repeated    = []
    outofdomain = []
    norestr = []
    iii = len(str(len(dataconfs)))
    jjj = max(len(str(vec0)),len("angles"))
    kkk = max(inpvars._ntorsions*2+2,len("conf"))
    if inpvars._ts:
       lengths = (iii,jjj,kkk,6,7,7,7,7,8)
       head = ("","angles","conf","weight","E","E+ZPE","Gibbs","X","ifreq")
    else:
       lengths = (iii,jjj,kkk,6,7,7,7,7)
       head = ("","angles","conf","weight","E","E+ZPE","Gibbs","X")
    table_head = " "+" | ".join(fncs.fill_string(col,length) for col,length in zip(head,lengths))+" "
    print(IBS2+"-"*len(table_head))
    print(IBS2+table_head)
    print(IBS2+"-"*len(table_head))
    nomenclatures = []
    QMSHO     = np.array([0.0 for temp in inpvars._temps])
    QMSHO_pre = np.array([0.0 for temp in inpvars._temps])
    QMSHO_sto = np.array([0.0 for temp in inpvars._temps])
    count_pre = 0
    count_sto = 0
    # calculate quotient of Xj
    xj_quo = 0.0
    beta   = 1.0/pc.KB/inpvars._temp
    for idx_j,conftuple_j in enumerate(dataconfs):
        vec_j,V0_j,V1_j,Gj,weight_j,Qrv_j,ifreq_j,lzmat_j,zmatvals_j,log_j = conftuple_j
        xj_quo += weight_j*np.exp(-(Gj-minG)*beta)
    # Now, one by one
    sumxj_pre = 0.0
    sumxj_sto = 0.0
    table2 = []
    for idx1,conftuple1 in enumerate(dataconfs):
        vec1,V0_1,V1_1,G1,weight1,Qrv1,ifreq1,lzmat1,zmatvals1,log1 = conftuple1
        # extra data
        try   : dipole  = read_dipole(folder+log1)
        except: dipole  = None
        try   : rotcons = sorted(list(get_rotcons(lzmat1,zmatvals1,symbols)))
        except: rotcons = None
        # Generate line
        str_number = "%%%ii"%iii%(idx1+1)
        # nomenclature
        nomenclature = vec1.nomenclature()
        idx = 1
        while True:
            idx += 1
            if nomenclature in nomenclatures: nomenclature = vec1.nomenclature()+"_%i"%idx
            else: break
        nomenclatures.append(nomenclature)
        # get relative energies
        relV0 = (V0_1-minV0)*pc.KCALMOL
        relV1 = (V1_1-minV1)*pc.KCALMOL
        relG  = (G1  -minG )*pc.KCALMOL
        xj    = weight1*np.exp(-(G1-minG)*beta)/xj_quo
        datainline  = [fncs.fill_string(str_number,iii)]
        datainline += [fncs.fill_string(str(vec1),jjj)]
        datainline += [fncs.fill_string(nomenclature,kkk)]
        datainline += ["%2i"%weight1]
        datainline += ["%7.3f"%relV0,"%7.3f"%relV1,"%7.3f"%relG,"%7.5f"%xj]
        if inpvars._ts:
            try   : sifreq = "%.1fi cm^-1"%abs(fncs.afreq2cm(ifreq1))
            except: sifreq = "  -  "
            datainline += [sifreq]
        table_line = " "+" | ".join(fncs.fill_string(col,length) \
                         for col,length in zip(datainline,lengths))
        table2.append( datainline[0:2]+[dipole,rotcons] )
        # Add to QMSHO
  ###   add_to_qmsho = weight1*Qrv1[str(vec1)]*np.exp(-(V1_1-minV1)/pc.KB/inpvars._temps)
        add_to_qmsho = weight1*Qrv1*np.exp(-(V1_1-minV1)/pc.KB/inpvars._temps)
        QMSHO += add_to_qmsho
        # preconditioned or stochastic?
        if log1.startswith("stoc"):
            table_line += "   "
            QMSHO_sto  += add_to_qmsho
            count_sto  += weight1
            sumxj_sto  += xj
        else:
            table_line += " * "
            QMSHO_pre  += add_to_qmsho
            count_pre  += weight1
            sumxj_pre  += xj
        # check geom
        bools = tfh.precheck_geom(lzmat1,zmatvals1,cmatrix,inpvars)
        if   (inpvars._tests[1][0] == 1) and (bools[0] is False):
           table_line += " <-- wrong connectivity"
           wrongconn.append(log1)
        elif (inpvars._tests[1][1] == 1) and (bools[1] is False):
           table_line += " <-- invalid domain"
           outofdomain.append(log1)
        elif (inpvars._tests[1][2] == 1) and (bools[2] is False):
           table_line += " <-- constraint not fulfilled"
           norestr.append(log1)
        elif (inpvars._tests[1][3] == 1) and (bools[3] is False):
           table_line += " <-- constraint not fulfilled"
           norestr.append(log1)
        else:
           same = False
           enantiomers = False
           # Append indices of point with same energy (within EPS_KCALMOL)
           idx_same = []
           for idx2 in range(idx1-1,-1,-1):
               vec2,V0_2,V1_2,G2,weight2,Qrv2,ifreq2,lzmat2,zmatvals2,log2 = dataconfs[idx2]
               diffE = abs(V0_2-V0_1)*pc.KCALMOL
               if diffE > EPS_KCALMOL: break
               idx_same.append(idx2)
           # Check angles
           for idx2 in idx_same:
               vec2,V0_2,V1_2,G2,weight2,Qrv2,ifreq2,lzmat2,zmatvals2,log2 = dataconfs[idx2]
               # check if same point or if enantiomer
               if vec1.is_same(vec2,inpvars._epsdeg)   : same = True       ; break
               if not inpvars._enantio                 : continue
               vec2_enantio = tfh.enantiovec(lzmat2,zmatvals2,dcorr,inpvars)
               if vec1.is_same(vec2_enantio,inpvars._epsdeg): enantiomers = True; break
           # indicate ir
           if same       : table_line += " <-- same as %s"%str(vec2)
           if enantiomers: table_line += " <-- enantio. of %s"%str(vec2)
           # save point to remove
           if same or enantiomers:
              # if possible, keep the preconditiones one
              if log1.startswith("prec.") and log2.startswith("stoc."): repeated.append(log2)
              else                                                    : repeated.append(log1)
        print(IBS2+table_line)
    print(IBS2+"-"*len(table_head))
    print(IBS2+"   * From preconditioned search")
    print("")
    nconfs = count_pre+count_sto
    print(IBS2+"   number of total          conformers: %3i"%nconfs)
    print(IBS2+"   number of preconditioned conformers: %3i (%6.2f%%); sum(X) = %.4f"%(count_pre,100.0*count_pre/nconfs,sumxj_pre))
    print(IBS2+"   number of stochastic     conformers: %3i (%6.2f%%); sum(X) = %.4f"%(count_sto,100.0*count_sto/nconfs,sumxj_sto))
    print("")
    pp.print_torsiangleclass()

    if len(repeated)+len(outofdomain)+len(norestr)+len(wrongconn) > 0:
       print(IBS2+"Number of wrong-connectivity points: %i"%len(wrongconn))
       print(IBS2+"Number of repeated           points: %i"%len(repeated))
       print(IBS2+"Number of out-of-domain      points: %i"%len(outofdomain))
       print(IBS2+"Number of restr-broken       points: %i"%len(norestr))
       folder2 = folder+DIREXCLUD
       print(IBS2+"  * these points will be moved to: %s"%folder2)
       mkdir_recursive(folder2)
       for log in repeated+outofdomain+norestr+wrongconn: move(folder+log,folder2+log)
       print(IBS2+"  * moved!")
       print("")
       # repeat again all
       del dataconfs
       classify_files(inpvars,cmatrix,case,dcorr)
    else:
       # Print MSHO table
       pp.print_mshotable(inpvars._temps,QMSHO,QMSHO_pre,QMSHO_sto,fscal)
       # Create molden file
       rw.write_molden_allconfs(dataconfs,Eref,symbols,folder)
       # table 2
       print(IBS2+"Total dipole moment (d.mom., debye) and rotational constants (Bi, GHz):")
       print("")
       lengths = (iii,jjj,7,7,7,7,7)
       head = ("","angles","d.mom.","B1","B2","B3")
       table_head = " "+" | ".join(fncs.fill_string(col,length) for col,length in zip(head,lengths))+" "
       print(IBS2+"-"*len(table_head))
       print(IBS2+table_head)
       print(IBS2+"-"*len(table_head))
       for idx,name,dipole,rotcons in table2:
           if dipole is None: dipole = "   -   "
           else             : dipole = "%7.3f"%dipole
           if rotcons is None: rotcons = "   -    |    -    |    -    "
           else : rotcons = "%7.3f | %7.3f | %7.3f "%tuple(rotcons)
           print("        "+idx+" | "+name+" | "+dipole+" | "+rotcons)
       print(IBS2+"-"*len(table_head))
       print("")
#==================================================#




