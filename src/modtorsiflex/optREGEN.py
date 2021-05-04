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
| Sub-module :  optREGEN           |
| Last Update:  2020/12/21 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''

#==================================================#
import os
from  shutil import copyfile
#--------------------------------------------------#
import modtorsiflex.tfgau   as tgau
import modtorsiflex.tfrw    as rw
import modtorsiflex.tfhelper     as tfh
from   modtorsiflex.tpespoint  import TorPESpoint
from   modtorsiflex.printing   import sprint
from   modtorsiflex.tfvars       import FDOMAINS, NIBS2, DOMAINS
#==================================================#

#==================================================#
def dirll_conformers(dirll,ts):
    # initialize domains
    ddomains = {domain:set([]) for domain in DOMAINS}

    # check folder
    if os.path.exists(dirll):
       for fpoint in tfh.folder_points(dirll):
           if ts: statusFRQ = 1
           else : statusFRQ = 0
           ddomains["conf"].add( (None,fpoint,statusFRQ) )
    else: tfh.create_dir(dirll)

    # print conformers
    confs = [fpoint for ipoint,fpoint,statusFRQ in ddomains["conf"]]
    if len(confs) != 0:
       sprint("Num conformers found in '%s': %i"%(dirll,len(confs)),NIBS2)
       sprint()

    # return data
    return ddomains, confs
#==================================================#



#==================================================#
# Regenerate domains from dirll and tmpll          #
#==================================================#
def regen_from_tmp(inpvars,cmatrix):
    # assert file does not exist
    if os.path.exists(FDOMAINS):
       sprint("File '%s' already exists!"%FDOMAINS,NIBS2)
       sprint()
       return

    #-------------------------------------#
    # Already located conformers in dirll #
    #-------------------------------------#
    ddomains, confs = dirll_conformers(inpvars._dirll,inpvars._ts)

        
    #---------------#
    # Data in tmpll #
    #---------------#
    if os.path.exists(inpvars._tmpll):
       # list of log files 
       opts = [log for log in os.listdir(inpvars._tmpll) \
               if log.startswith("opt") and log.endswith(".log")]
       opts.sort()
       sprint("Num opt log files in '%s': %i"%(inpvars._tmpll,len(opts)),NIBS2)
       repe_or_excl = []
       for log in opts:
           vec1     = TorPESpoint(log.split(".")[1])
           if log.startswith("optprec."):
              case = "prec"
              log_frq  = log.replace("optprec.","frqprec.")
           elif log.startswith("optstoc."):
              case = "stoc"
              log_frq  = log.replace("optstoc.","frqstoc.")
           # add folder
           log_opt  = inpvars._tmpll+log
           log_frq  = inpvars._tmpll+log_frq
           # Read log file and get status
           statusOPT,statusFRQ,vec2 = tgau.status_from_log(log_opt,cmatrix,inpvars)

           # 1. Log file without Normal termination
           if statusOPT == -1:
              ddomains["fail"].add( (vec1,None,-1) )

           # 2. Log file with Normal termination & frq file exists
           elif os.path.exists(log_frq):
              statusFRQ = tgau.status_from_log(log_frq,cmatrix,inpvars)[1]
              bool_min = (statusFRQ == 0) and (not inpvars._ts)
              bool_ts  = (statusFRQ == 1) and (    inpvars._ts)
              # 2.1. correct number of ifreqs
              if bool_min or bool_ts:
                 bool_in,closest = vec2.is_in_list(confs,inpvars._epsdeg)
                 if bool_in: ddomains["repe"].add( (vec1,vec2,statusFRQ) )
                 else:
                     ddomains["conf"].add( (vec1,vec2,statusFRQ) )
                     # copy file to dirll
                     log_conf = inpvars._dirll+"%s.%s.log"%(case,str(vec2))
                     if not os.path.exists(log_conf):
                        sprint("  %s  -->  %s"%(log_frq,log_conf),NIBS2)
                        copyfile(log_frq,log_conf)
                     # update confs
                     confs.append(vec2)
              # 2.2. wrong number of ifreqs
              else: ddomains["wimag"].add( (vec1,vec2,statusFRQ) )

           # 3. Log file with Normal termination & frq file does not exist
           else:
              repe_or_excl.append( (vec1,vec2) )

       #------------------------#
       # deal with repe or excl #
       #------------------------#
       for vec1,vec2 in repe_or_excl:
           domain = "excl"
           if vec2.is_in_list(confs,inpvars._epsdeg): domain = "repe"
           ddomains[domain].add( (vec1,vec2,-1) )
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
    rw.write_domains(len(inpvars._ttorsions),ddomains)
    sprint()
#==================================================#

