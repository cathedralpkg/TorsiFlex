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
| Sub-module :  help               |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#==================================================#
from   modtorsiflex.version import PROGVER
from   modtorsiflex.version import PROGNAME
from   modtorsiflex.version import PROGNAMEnopy
from   modtorsiflex.version import VERSION
from   modtorsiflex.version import DATE
from   modtorsiflex.tfvars  import FDOMAINS
from   modtorsiflex.tfvars  import GAUTEMPL
#==================================================#

SPACE = " "*len(PROGNAME)
args  = (PROGNAMEnopy,PROGVER,PROGNAME,SPACE,SPACE,SPACE,SPACE,FDOMAINS,GAUTEMPL,PROGNAMEnopy)

HSTRING = '''
Current version: %s %s

Description:

  A program to seek the conformers of flexible molecules
  by means of a combined preconditioned-stochastic algorithm.
  The located conformers, calculated with a low-level (LL)
  electronic structure method, can be re-optimized using
  a high-level (HL) method.

Main authors:

  Dr.   David Ferro-Costas
  Prof. Antonio Fernandez-Ramos

  Centro Singular de Investigacion en Quimica Bioloxica
  e Materiais Moleculares (CIQUS), Universidade de
  Santiago de Compostela, Galicia, Spain

Execution:

  %s [--help /-h]  [--version/-v]
  %s [--smiles  ]  [--cartesian ]  [--inp/--input]
  %s [--prec    ]  [--stoc      ]  [--hlopt      ]
  %s [--msho    ]  [--mstor     ]
  %s [--torsions]  [--regen     ]

-------------------
  Program options  
-------------------

   --smiles smiles_code zmatfile.zmat

     Converts the SMILES code to the Z-matrix format (stored at zmatfile.zmat)
    
     Requires two arguments:

     - 1st argument --> smiles_code

       The smiles code of the molecule (given between quotation marks);

     - 2nd argument --> zmatfile.zmat

       The name of the file where the Z-matrix will be stored

     Example:

       --smiles "CCCCO" butanol.zmat


   --cartesian ccfile.xyz zmfile.zmat

     Converts Cartesian coordinates (stored in ccfile.xyz)
     to the Z-matrix format (stored at zmfile.zmat)

     Example:

       --cartesian butanol.xyz butanol.zmat


   --input [zmfile.zmat]

     Generates the input file and the templates for Gaussian.

     Admits one argument, zmfile.zmat, the name of the file containing
     the reference Z-matrix of the system.


   --prec [M m]

     Uses the preconditioned algorithm for the conformer location.

     To divide the preconditioned guesses into M groups and
     deal with the m-th group, use this option as follows:

          --hlopt M m

     For example:

          torsiflex --hlopt 10 2

     divides the guesses into 10 groups and only carries out
     the calculations associated to the 2nd group.

     --> Calculations with Gaussian are carried out <--


   --stoc [point(s)]

     Uses the stochastic algorithm for the conformer location.

     Specific points of the torsional space can be calculated by
     introducing them as arguments. For example, in a system with
     three target torsions, points (90,30,55) and (20,55,89), defined
     in degrees, can be calculated with:

        --stoc  90_33_55  20_55_89

     --> Calculations with Gaussian are carried out <--


   --hlopt [n1 [n2]] [nocalc]

     Re-optimizes LL conformers at HL.

     The n1-th LL conformer can be optimized with

          --hlopt n1
     
     whereas conformers n1-th to n2-th LL (both included) can
     be re-optimized at HL using

          --hlopt n1 n2

     For example, 

          torsiflex --hlopt 5 10

     considers conformers 5 to 10, both included. Notice that conformers
     are listed by increasing order of energy (see --msho option).

     When followed by 'nocalc':

          --hlopt nocalc

     it generates the gjf files (Gaussian inputs) without
     carrying out calculations. Useful to send the calculations
     on your own.

     --> Calculations with Gaussian are carried out <--


   --msho [ll/hl]

     Checks the located conformers, listing them by increasing
     order of energy and calculates the multi-structure
     harmonic-oscillator (MS-HO) partition functions.

     This option can be carried out exclusively for the
     low-level (ll) or the high-level (hl) conformers if
     followed by the corresponding abbreviation:

          --msho ll
          --msho hl


   --mstor [ll/hl]

     Generates the MsTor input files.

     This option can be carried out exclusively for the
     low-level (ll) or the high-level (hl) conformers if
     followed by the corresponding abbreviation:

          --mstor ll
          --mstor hl


   --torsions

     Adds or removes a target torsion (must be defined in
     the Z-matrix file).


   --regen

     Regenerates the %s file using the temporal files.


   --help (also -h)

     Prints this help message.


   --version (also -v)

     Prints the program version.

---------------------
  Extra Information  
---------------------

(a) Assert the path to the Gaussian executable
    is defined in your .bashrc file under the name
    'GauExe' and export it. Example:

    export GauExe='/home/programs/Gaussian/g09'

(b) Modify the Gaussian templates (%s file),
    taking into account that:

    * [nproc], [mem], [level], [optmode], [charge],
      [multipl], [zmat], [modred], and [fccards]
      are %s indications.
      They should not be removed.

    * Gaussian command line must start with '#p'

    * The command line includes 'iop(99/9=1,99/14=3)'.
      This is mandatory and should not be deleted:
         99/9=1  --> rotates to Z-matrix orientation first
         99/14=3 --> expresses final optimized structure
                     in terms of the input Z-matrix

    * 'scf=(incore)' is recommended with Hartree-Fock
      calculations. With it, Gaussian stores the full
      integral list in memory, speeding up the calculation.

(c) We highly recommend to firstly carry out the
    preconditioned search, which should found the
    chemically-intuitive conformers.
    After it, the stochastic algorithm should be used.
'''%args


