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
| Sub-module :  optREGEN           |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''

#==================================================#
import os
from  shutil import copyfile
#--------------------------------------------------#
from   common.files            import write_zmat,write_gtsfile, read_zmat
#--------------------------------------------------#
import modtorsiflex.tfgau      as tgau
from   modtorsiflex.tfhelper   import enantiovec
from   modtorsiflex.tfrw       import write_domains, gts2molden
from   modtorsiflex.tpespoint  import TorPESpoint
from   modtorsiflex.printing   import sprint
from   modtorsiflex.tfvars     import FDOMAINS, NIBS2, DOMAINS
#==================================================#

#==================================================#
def ddomains_from_dirll(inpvars,dcorr):
    '''
    Lists conformers in dirll and add them to the
    domains dictionary
    '''
    # initialize dictionary
    ddomains = {domain:set([]) for domain in DOMAINS}
    # No folder --> return empty variables
    if not os.path.exists(inpvars._dirll): return ddomains
    # Status depending of type of stationary point
    if inpvars._ts: statusFRQ = 1
    else          : statusFRQ = 0
    # list zmat files
    files = [fname for fname in os.listdir(inpvars._dirll) if fname.endswith(".zmat")]
    # Go one by one
    for zmatfile in sorted(files):
        vec      = TorPESpoint(zmatfile.split(".")[1])
        evec     = None
        ddomains = add_to_domain(ddomains,"conf",statusFRQ,vec,vec)
        # If enantio, read file and get enantio-vec
        if inpvars._enantio:
           lzmat,zmatvals,zmatatoms = read_zmat(inpvars._dirll+zmatfile)[0]
           evec     = enantiovec(lzmat,zmatvals,dcorr,inpvars)
           ddomains = add_to_domain(ddomains,"enan",statusFRQ,evec,evec)
    # Return data
    return ddomains
#--------------------------------------------------#
def ddomains_from_scratch(ddomains,inpvars,dcorr,confs=[]):
    # list of log files 
    opts = sorted([log for log in os.listdir(inpvars._tmpll) \
                        if log.startswith("opt") and log.endswith(".log")])
    sprint("Num opt log files in '%s': %i"%(inpvars._tmpll,len(opts)),NIBS2,1)
    # Loop over log list
    repe_or_excl = []
    for log in opts:
        vec1    = TorPESpoint(log.split(".")[1])
        log_frq = log.replace("opt","frq")
        if   log.startswith("optprec."): scase = "prec"
        elif log.startswith("optstoc."): scase = "stoc"
        else                           : continue

        # add folder
        log_opt  = inpvars._tmpll+log
        log_frq  = inpvars._tmpll+log_frq

        # Read log file and get status
        statusOPT,statusFRQ,vecs1,vecs2 = tgau.status_from_log(log_opt,inpvars,dcorr)
        evec1 = vecs1[1]
        vec2  = vecs2[0]
        evec2 = vecs2[1]

        # 1. Log file without Normal termination
        if statusOPT == -1:
           ddomains = add_to_domain(ddomains,"fail",-1,vec1,None,evec1,None)

        # 2. Log file with Normal termination & frq file exists
        elif os.path.exists(log_frq):
           statusFRQ = tgau.status_from_log(log_frq,inpvars,dcorr)[1]
           bool_min = (statusFRQ == 0) and (not inpvars._ts)
           bool_ts  = (statusFRQ == 1) and (    inpvars._ts)
           # 2.1. correct number of ifreqs
           if bool_min or bool_ts:
              bool_in,closest = vec2.is_in_list(confs,inpvars._epsdeg)
              if bool_in:
                 ddomains = add_to_domain(ddomains,"repe",statusFRQ,vec1,vec2,evec1,evec2)
              else:
                 ddomains = add_to_domain(ddomains,"conf",statusFRQ, vec1, vec2)
                 ddomains = add_to_domain(ddomains,"enan",statusFRQ,evec1,evec2)
                 # Convert to gts
                 gts_conf = inpvars._dirll+"%s.%s.gts"%(scase,str(vec2))
                 if not os.path.exists(gts_conf):
                    sprint("  %s  -->  %s"%(log_frq,gts_conf),NIBS2)
                    if not os.path.exists(inpvars._dirll): os.mkdir(inpvars._dirll)
                    copy_conformer(log_frq,gts_conf,inpvars)
                    sprint()
                 # update confs
                 confs.append( vec2)
                 confs.append(evec2)

           # 2.2. wrong number of ifreqs
           else:
               ddomains = add_to_domain(ddomains,"wimag",statusFRQ,vec1,vec2,evec1,evec2)

        # 3. Log file with Normal termination & frq file does not exist
        else:
           repe_or_excl.append( ( vec1, vec2) )
           repe_or_excl.append( (evec1,evec2) )

    #------------------------#
    # deal with repe or excl #
    #------------------------#
    confs = [conf for conf in confs if conf is not None]
    for vec1,vec2 in repe_or_excl:
        domain = "excl"
        if vec2.is_in_list(confs,inpvars._epsdeg): domain = "repe"
        ddomains = add_to_domain(ddomains,domain,-1,vec1,vec2)
    return ddomains
#--------------------------------------------------#
def add_to_domain(ddomains,which,statusFRQ,vec1,vec2,evec1=None,evec2=None):
    # add initial and final vectors
    if vec1 is None and vec2 is None: return ddomains
    ddomains[which].add( (vec1,vec2,statusFRQ) )
    # add initial and final vectors (enantiomers)
    if evec1 is None and evec2 is None: return ddomains
    ddomains[which].add( (evec1,evec2,statusFRQ) )
    # return dictionary
    return ddomains
#--------------------------------------------------#
def copy_conformer(ilog,ogts,inpvars):
    ozmat   = ogts[:-3]+"zmat"
    omolden = ogts[:-3]+"molden"

    statusFRQ,vec,lzmat,zmatvals,zmatatoms,basic = tgau.read_normal_log(ilog,inpvars,"LL")

    # create z-matrix file
    write_zmat(ozmat,lzmat, zmatvals)
    write_gtsfile(*tuple(list(basic)+[ogts]))
    gts2molden(ogts,omolden)
#==================================================#



#==================================================#
# Regenerate domains from dirll and tmpll          #
#==================================================#
def regen_from_tmp(inpvars,dcorr):
    '''
    Regenerate domains.txt file from scratch data
    '''

    # assert file does not exist
    if os.path.exists(FDOMAINS):
       sprint("File '%s' already exists!"%FDOMAINS,NIBS2)
       sprint()
       return

    #-------------------#
    # Use data in dirll #
    #-------------------#
    ddomains = ddomains_from_dirll(inpvars,dcorr)
    confs    = [fpoint for ipoint,fpoint,statusFRQ in ddomains["conf"]]
    confs   += [fpoint for ipoint,fpoint,statusFRQ in ddomains["enan"]]

    # print number of conformers
    if len(confs) != 0:
       sprint("Num conformers found in '%s': %i"%(inpvars._dirll,len(confs)),NIBS2)
       sprint()

    #-------------------#
    # Use data in tmpll #
    #-------------------#
    if os.path.exists(inpvars._tmpll):
       ddomains = ddomains_from_scratch(ddomains,inpvars,dcorr,confs)
       sprint()

    #---------------------------------#
    # Print num points in each domain #
    #---------------------------------#
    for domain in DOMAINS:
        np = len(ddomains[domain])
        str_domain = "'%s'"%domain
        sprint("Points in domain %-8s: %i"%(str_domain,np),NIBS2)
    sprint()

    #----------------------#
    # rewrite domains file #
    #----------------------#
    sprint("Writing file: %s"%FDOMAINS,NIBS2)
    sprint()
    write_domains(len(inpvars._ttorsions),ddomains)
    sprint()
#==================================================#

