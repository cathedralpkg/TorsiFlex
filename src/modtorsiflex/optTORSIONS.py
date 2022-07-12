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
| Sub-module :  optTORSIONS        |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''

#==================================================#
import os
import readline
import shutil
#--------------------------------------------------#
import common.Exceptions       as exc
import common.gaussian        as gau
from   common.files       import read_zmat
from   common.MyCompleter import MyCompleter
#--------------------------------------------------#
from   modtorsiflex.inpvars      import InpVars
import modtorsiflex.tfgau      as tgau
import modtorsiflex.tfrw       as rw
import modtorsiflex.tfhelper   as tfh
from   modtorsiflex.tpespoint  import TorPESpoint, fvec2svec
from   modtorsiflex.printing   import sprint
import modtorsiflex.tfvars     as tvars
from   modtorsiflex.version    import PROGNAMEnopy
from   modtorsiflex.optREGEN   import regen_from_tmp
#==================================================#


#===============================================================#
# Activate autocomplete
#===============================================================#
def set_completer(options=[]):
    completer = MyCompleter(options)
    readline.set_completer(completer.complete)
    readline.parse_and_bind('tab: complete')
#===============================================================#


#==================================================#
def rewrite_input(inpvars,initial,final):
    if not os.path.exists(tvars.IFILE): return

    print("       - modifying file (%s)..."%tvars.IFILE)

    precond = []
    tdomain = []
    tsigma  = []

    with open(tvars.IFILE,'r') as asdf: lines = asdf.readlines()
    new = []
    for idx_new,torsion in enumerate(final):
        if torsion in initial:
            X_old = inpvars._ttorsions[initial.index(torsion)]
            idx_old = int(X_old)-1
        else:
            X_old = None
            idx_old = None
            
        if X_old in inpvars._ttorsions:
            precond_i = inpvars._precond[idx_old]
            tdomain_i = inpvars._tdomain[idx_old]
            tsigma_i  = inpvars._tsigma[idx_old]
        else:
            precond_i = [60,180,300]
            tdomain_i = [(0,360)]
            tsigma_i  = 1
        sdomain_i = "U".join("(%i,%i)"%(phi1,phi2) for phi1,phi2 in tdomain_i)

        new.append("torsion%s    %s"%(idx_new+1,torsion))
        new.append("precond%s    %s"%(idx_new+1," ".join([str(angle) for angle in precond_i])))
        new.append("tsigma%s     %i"%(idx_new+1,tsigma_i))
        new.append("tdomain%s    %s"%(idx_new+1,sdomain_i))

    # Add
    STRING = ""
    IDX = None
    for idx,line in enumerate(lines):
        if idx-1 == IDX:
           for line in new: STRING += line+"\n"
        line = line.strip()
        avoid = ["torsion","tsigma","precond","tdomain"]
        skip = False
        for key in avoid:
            if line.startswith(key): skip = True
        if skip: continue
        STRING += line+"\n"
        # locate where to put new info
        if line.startswith("#    Target t"):
            IDX = idx
            if lines[idx+1].startswith("#-----"): IDX += 1
            else: STRING += "#-----------------------#\n"
    if IDX is None:
       STRING += "#-----------------------#\n"
       STRING += "#    Target torsions    #\n"
       STRING += "#-----------------------#\n"
       for line in new: STRING += line+"\n"

    with open(tvars.IFILE,'w') as asdf: asdf.write(STRING)
#--------------------------------------------------#
def svec_from_zmatvals(zmatvals,final):
    vec = [zmatvals[ic]%360 for ic in final]
    return fvec2svec(vec)
#--------------------------------------------------#
def svec_from_gjf(gjf,final):
    with open(gjf,'r') as asdf: lines = asdf.readlines()
    ivalues = {}
    for line in lines:
        line = line.strip()
        line = line.split()
        if len(line) != 2: continue
        ic,value = line
        if ic not in final: continue
        ivalues[ic] = float(value)
    assert len(final) == len(ivalues)
    return svec_from_zmatvals(ivalues,final)
#--------------------------------------------------#
def rename_files_in_tmp(tmp_folder,final):
    if not os.path.exists(tmp_folder): return
    print("       - modifying files in folder (%s)..."%tmp_folder)
    new_files = []
    # Files in folder
    gjfs = [gjf for gjf in os.listdir(tmp_folder) if gjf.endswith(".gjf") and gjf.startswith("opt")]
    # modify files in scratch
    for gjf in gjfs:
        # -- opt files -- #
        gjf1 = gjf
        err1 = gjf1.replace(".gjf",".err")
        log1 = gjf1.replace(".gjf",".log")
        # -- frq files -- #
        gjf2 = gjf1.replace("opt","frq")
        err2 = err1.replace("opt","frq")
        log2 = log1.replace("opt","frq")
        # current and new svecs
        svecA = gjf.split(".")[1]
        svecB = svec_from_gjf(tmp_folder+gjf,final)
        # rename file(s) - adding "tmp" prefix
        for fnameA in [gjf1,err1,log1,gjf2,err2,log2]:
            fnameB = "tmp-"+fnameA.replace(svecA,svecB)
            if not os.path.exists(tmp_folder+fnameA): continue
            os.rename(tmp_folder+fnameA,tmp_folder+fnameB)
            new_files.append(fnameB)
    # (c) remove "tmp" prefix
    for fname in new_files:
        final_name = fname.replace("tmp-","")
        if os.path.exists(tmp_folder+fname):
           os.rename(tmp_folder+fname,tmp_folder+final_name)
#--------------------------------------------------#
def rename_files_in_files(files_folder,final):
    if not os.path.exists(files_folder): return
    print("       - modifying files in folder (%s)..."%files_folder)
    # Files in folder
    names = [".".join(fname.split(".")[:-1]) for fname in os.listdir(files_folder) \
                                             if  fname.endswith(".zmat") ]
    # Remove repetitions
    names = set(names)
    # modify files in scratch
    new_files = []
    for name in names:
        # svec from name
        svecA    = name.split(".")[1]
        # name of files
        zmatfile = name + ".zmat"
        molden   = name + ".molden"
        gts      = name + ".gts"
        # zmatvals from file
        if not os.path.exists(files_folder+zmatfile): continue
        zmatvals = read_zmat(files_folder+zmatfile)[0][1]
        # Final vector
        svecB = svec_from_zmatvals(zmatvals,final)
        #------------------------------------#
        # rename files - adding "tmp" prefix #
        #------------------------------------#
        for iname in [gts,zmatfile,molden]:
            if not os.path.exists(files_folder+iname): continue
            fname = files_folder+"tmp-"+iname.replace(svecA,svecB)
            iname = files_folder+iname
            os.rename(iname,fname)
            new_files.append(fname)
    # (c) remove "tmp" prefix
    for tmpname in new_files:
        if not os.path.exists(tmpname): continue
        fname = tmpname.replace("tmp-","")
        os.rename(tmpname,fname)
#==================================================#



#==================================================#
# 
#==================================================#
def redefine_torsions(inpvars,dcorr,zmat,cmatrix,symbols):

    lzmat, zmatvals, zmatatoms = zmat

    # List torsions in system
    all_torsions = {}
    for coord,atoms in zmatatoms.items():
        if len(atoms) != 4: continue
        if not tfh.is_proper_torsion(atoms,cmatrix): continue
        storsion = "-".join([symbols[atom]+"%i"%(atom+1) for atom in atoms])
        all_torsions[coord] = storsion

    # List selected torsions
    initial = list(inpvars._tic)
    final   = list(initial)

    # Print info
    set_completer(["remove","add","ok","exit"])
    while True:
        print("    Torsions in system:")
        for torsion,storsion in all_torsions.items():
            line = "       %s (%s)"%(torsion,storsion)
            if torsion in final:
               idx = final.index(torsion)
               line += " <-- torsion%i"%(idx+1)
            print(line)
        print("")
        print("    to remove a torsion    : remove name_of_torsion")
        print("    to add    a torsion    : add    name_of_torsion")
        print("    to confirm             : ok")
        print("    to exit without changes: exit")
        print("")
        answer = input("    > ").strip()
        print("")
        if answer.lower().startswith("exit"): raise exc.END
        if answer.lower().startswith("ok"  ): break
        action,target = answer.split()
        if target not in all_torsions: continue
        if action == "remove" and target     in final: final.remove(target)
        if action == "add"    and target not in final: final.append(target)

    #------------------#
    # CARRY ON CHANGES #
    #------------------#
    if initial == final:
        print("    Torsions remain the same!")
        print("")
        raise exc.END
    print("    Torsions have been modified!")
    print("")
    print("    The following files will be modified:")
    print("    - %s"%(tvars.IFILE))
    print("")
    print("    The following files will be renamed:")
    print("    - files inside '%s'"%(inpvars._tmpll))
    print("    - files inside '%s'"%(inpvars._tmphl))
    print("    - files inside '%s'"%(inpvars._dirll))
    print("    - files inside '%s'"%(inpvars._dirhl))
    print("")
    print("    The following files will be removed:")
    print("    - %s"%(inpvars._dirll+tvars.ALLCONFS))
    print("    - %s"%(inpvars._dirhl+tvars.ALLCONFS))
    print("    - %s"%(tvars.FDOMAINS))
    print("")
    print("    The following folders will be removed:")
    print("    - %s"%(inpvars._dirll+tvars.DIRMSTOR))
    print("    - %s"%(inpvars._dirhl+tvars.DIRMSTOR))
    print("")
    print("    Apply changes (y/n)?")
    print("")
    set_completer(["y","n"])
    answer = input("    > ").strip().lower()
    print("")
    if answer not in ["y","yes"]: raise exc.END

    # remove CONFORMERS.txt
    files2remove  = [inpvars._dirll+tvars.ALLCONFS,inpvars._dirhl+tvars.ALLCONFS]
    files2remove += [tvars.FDOMAINS,tvars.MINGTXT,tvars.ENERGYSUMLL,tvars.LLHLCORRFILE]
    for fname in files2remove:
        if not os.path.exists(fname): continue
        print("       - removing file (%s)"%fname)
        os.remove(fname)

    # remove MSTOR folders
    folders2remove = [inpvars._dirll+tvars.DIRMSTOR,inpvars._dirhl+tvars.DIRMSTOR]
    for folder in folders2remove:
        if not os.path.exists(folder): continue
        print("       - removing folder (%s)"%folder)
        shutil.rmtree(folder)

    # modify torsiflex.inp
    rewrite_input(inpvars,initial,final)

    # modify files in tmp folder
    rename_files_in_tmp(inpvars._tmpll,final)
    rename_files_in_tmp(inpvars._tmphl,final)

    # modify files in the conformer folders
    rename_files_in_files(inpvars._dirll,final)
    rename_files_in_files(inpvars._dirhl,final)

    # Now, regenerate from scratch
    if os.path.exists(tvars.IFILE):
       print("")
       print("       - regenerating domains...")
       inpvars = InpVars()
       inpvars.read_input(tvars.IFILE)
       inpvars.prepare_variables()
       for X,ic in zip(inpvars._ttorsions,inpvars._tic):
           inpvars._tatoms[X] = list(zmatatoms[ic])
       args = (inpvars,dcorr)
       regen_from_tmp(*args)
#==================================================#

