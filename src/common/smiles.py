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
| Sub-module :  smiles             |
| Last Update:  2022/04/23 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''

#==============================================================#
from importlib.util  import find_spec
#--------------------------------------------------------------#
from common.physcons import ANGSTROM
#--------------------------------------------------------------#
specs = []
for pkg in ["rdkit","rdkit.Chem"]:
    try   : specs.append( find_spec(pkg) )
    except: specs.append( None)
if None not in specs:
   from   rdkit      import Chem
   from   rdkit.Chem import AllChem
#==============================================================#

#==============================================================#
def smiles_to_cc(smiles):
    '''
    Converts smiles code to Cartesian coordinates

    inp1: smiles ; a string with the smiles code

    out1: xcc    ; a list of floats with the cartesian
                   coordinates in bohr
    out2: symbols; a list of strings with the atomic symbols
    '''
    # initialize lists
    xcc,symbols = [],[]
    # Create from code
    m = Chem.MolFromSmiles(smiles)
    # Add H atoms
    m = Chem.AddHs(m)
    # Now get Cartesian Coordinates
    AllChem.EmbedMolecule(m)
    sdata  = Chem.MolToMolBlock(m)
    for line in sdata.split("\n"):
        data = line.split()
        if len(data) != 16: continue
        x,y,z,symbol = data[0:4]
        symbols.append(symbol)
        xcc += [x,y,z]
    # From Ansgstrom to Bohr
    xcc = [float(xi)/ANGSTROM for xi in xcc]
    # Return data
    return xcc, symbols
#==============================================================#


