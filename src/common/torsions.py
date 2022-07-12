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
| Module     :  common             |
| Sub-module :  torsions           |
| Last Update:  2022/07/04 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

'''


#=============================================#
import os
import numpy as np
#---------------------------------------------#
import common.Exceptions as     Exc
from   common.fncs       import xyz
from   common.fncs       import clean_lines
from   common.fncs       import atonums2masses
from   common.files      import read_file
#=============================================#


#=============================================#
def torsion_of_cx3(atoms,symbols,cmatrix,out=0):
    at1,at2,at3,at4 = atoms
    # torsion (X)3-C-A-...
    for C,A in [(at2,at3),(at3,at2)]:
        # Check C is a C atom
        if symbols[C] != "C": continue
        # Check C atom is bonded to 3 atoms (excluding A)
        Xs = [X for X,bonded in enumerate(cmatrix[C]) if bonded and X != A]
        if len(Xs) != 3: continue
        # Check the 3 Xs are equal
        symbols_Xs = set([symbols[X] for X in Xs])
        if len(symbols_Xs) != 1: continue
        # Check each X is only bonded to C
        bonded_Xs = set([list(cmatrix[X]).count(True) for X in Xs])
        if len(bonded_Xs)  != 1: continue
        if bonded_Xs.pop() != 1: continue
        # WE DO HAVE A CX3!!
        if out == 0: return True
        else       : return Xs
    # No CX3 group involved in rotation
    if out == 0: return False
    else       : return None
#=============================================#


