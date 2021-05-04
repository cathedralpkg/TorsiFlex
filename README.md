# _TorsiFlex_

## About _TorsiFlex_

    Name of the Program: TorsiFlex
    Program Version : 2021.2
    Program Version Date: May 04, 2021
    Manual  Version Date: May 04, 2021

_TorsiFlex_ is an user-friendly program written in Python 3.
It was designed to seek the conformers of a given molecule by
 using a combined low-level/high-level (LL/HL) methodology.



## How to cite

D. Ferro-Costas and A. Fernández-Ramos, Front. Chem., 2020, 8:16.


## Licensing and Distribution 

_TorsiFlex version 2021.2_

MIT LICENSE

Copyright (c) 2021, David Ferro Costas (david.ferro@usc.es) and Antonio Fernandez Ramos (qf.ramos@usc.es)

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


## Description of files

 Contents of the folders distributed in this version:
  - **src/**       : _TorsiFlex_ source files
  - **docs/**      : Manual of _TorsiFlex_
  - **tests/**     : All the files related to the tests set
        

## Installation

_TorsiFlex_ is a program written in Python 3. Consequently, it does not need any kind 
of compilation, as it would be the case with C or Fortran programs.
The user should install Python 3 in order to use _TorsiFlex_, 
as well as the following Python libraries:
   - cmath
   - fcntl
   - glob
   - math
   - matplotlib
   - multiprocessing
   - numpy
   - os
   - random
   - scipy
   - sys
   - time

WARNING: __do not__ use Python 2 to execute _TorsiFlex_.


## Setting up the program

Before using _TorsiFlex_, the user has to define the path to the executable(s) of the 
software for the electronic structure calculation (ESSO).
At the moment _TorsiFlex_ supports __Gaussian__ as ESSO.

In order to interact with the ESSO, _TorsiFlex_ needs to know the location of some executable files. 
Such information is obtained from the following environment variables, which have to be 
defined and exported by the user in their __.bashrc__ file:

For __Gaussian__, the environment variable is:

  - _GauExe_, the path to the Gaussian executable and

Example of paths are:

  - ```export Gauexe="/home/programs/G09_64D/g09/g09"``` 


## Execution

You can run _TorsiFlex_ by invoking the Python interpreter manually as follows:

```python3 torsiflex.py```

If you prefer to avoid invoking the Python interpreter, you have to follow these two simple steps:

Add as the first line in the torsiflex.py file the following:

```#!PATH_FOR_PYTHON python```

where PATH_FOR_PYTHON indicates the location of the Python interpreter.

Example:

```#!/usr/bin/python3```


In this example Python is located in `/usr/bin/`.

Make the main program torsiflex.py executable:

```chmod u+x torsiflex.py```

This allows you to run _TorsiFlex_ just using:

```torsiflex.py```

Before run _TorsiFlex_, we recommend to read the help menu. It can be displayed either by typing

```torsiflex.py --help```

or

```torsiflex.py -h```

## Tests set

Directory __tests/__ contains the z-matrix files associated to the three worked examples (see manual). 

                                                            

