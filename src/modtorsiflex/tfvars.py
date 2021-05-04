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
| Sub-module :  variables          |
| Last Update:  2021/02/25 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains different
values for some variables used in TorsiFlex
'''

#==========================================#
PROGNAME      = "torsiflex.py"
PROGNAMEnopy  = "TorsiFlex"
VERSION       = "2021.2"
DATE          = "2021-05-04"
PROGVER       = "v%s (%s)"%(VERSION,DATE)
#------------------------------------------#
# initial blank spaces to add to all print statements
NIBS  = 4 ; IBS  = " "*NIBS
NIBS2 = 7 ; IBS2 = " "*NIBS2
#------------------------------------------#
# name of folders of files
DIRTEMPL      = "GauTemplates/"
DIREXCLUD     = "Excluded/"
DIRMSTOR      = "MsTor/"
# name of files
IFILE         = "torsiflex.inp"
MINGTXT       = "minG.txt"
LRFILE        = "lastrun"
DOMAINS       = "conf,enan,repe,fail,excl,wimag".split(",")
FDOMAINS      = "domains.txt"
MSTORIF1      = DIRMSTOR+"mstor.dat"
MSTORIF2      = DIRMSTOR+"hess.dat"
# Gaussian templates
TEMPLMINOPTLL = DIRTEMPL+"min_optLL"
TEMPLMINOPTHL = DIRTEMPL+"min_optHL"
TEMPLMINFRQLL = DIRTEMPL+"min_frqLL"
TEMPLMINFRQHL = DIRTEMPL+"min_frqHL"
TEMPLTSOPTLL  = DIRTEMPL+"ts_optLL"
TEMPLTSOPTHL  = DIRTEMPL+"ts_optHL"
TEMPLTSFRQLL  = DIRTEMPL+"ts_frqLL"
TEMPLTSFRQHL  = DIRTEMPL+"ts_frqHL"
# File with all conformers
ALLCONFS      = "CONFORMERS.xyz"
# Correlation LL-HL
LLHLCORRFILE  = "corr_llhl.txt"
#==========================================#
EPS_KCALMOL   = 0.100 # energy difference to compare geoms (kcal/mol)
#==========================================#


