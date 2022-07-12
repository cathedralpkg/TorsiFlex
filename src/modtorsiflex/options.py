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
| Sub-module :  options            |
| Last Update:  2022/07/12 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*
'''


#==================================================#
import os
import sys
#--------------------------------------------------#
import common.Exceptions         as exc
#--------------------------------------------------#
import modtorsiflex.printing     as pp
#==================================================#


#==================================================#
OPTIONS1  = "h,v"
#--------------------------------------------------#
OPTIONS2  = "help,version,"
OPTIONS2 += "smiles,cartesian,input,inp,"
OPTIONS2 += "prec,stoc,hlopt,"
OPTIONS2 += "msho,mstor,"
OPTIONS2 += "regen,torsions"
#==================================================#

OPTIONS = ["-" +opt for opt in OPTIONS1.split(",")] + \
          ["--"+opt for opt in OPTIONS2.split(",")]

#==================================================#
def prepare_smiles(arguments):
    # Evaluate arguments
    if   len(arguments) == 2: smiles,fname = arguments
    else                    : raise Exception
    # add extension
    if not fname.endswith(".zmat"): fname += ".zmat"
    # return data
    return smiles,fname
#--------------------------------------------------#
def prepare_cartesian(arguments):
    # Evaluate arguments
    if   len(arguments) == 2: file1,file2 = arguments
    else                    : raise Exception
    # add extension
    if not file2.endswith(".zmat"): file2 += ".zmat"
    # return data
    return file1,file2
#--------------------------------------------------#
def prepare_input(arguments):
    # Default value
    zmatfile = None
    # Evaluate arguments
    if   len(arguments) == 1: zmatfile = arguments[0]
    elif len(arguments) >= 2: raise Exception
    # zmatfile (get it from files)
    if zmatfile is None:
       zmats = [fname for fname in os.listdir(".") if fname.endswith(".zmat")]
       if   len(zmats) != 0: zmatfile = zmats[0]
       else                : zmatfile = "molecule.zmat"
    # return data
    return zmatfile
#--------------------------------------------------#
def prepare_prec(arguments):
    # Default value
    def_prec1 = 1
    def_prec2 = 1
    # Evaluate arguments
    if   len(arguments) == 0: prec1,prec2 = def_prec1,def_prec2
    elif len(arguments) == 2: prec1,prec2 = arguments
    else                    : raise Exception
    prec1,prec2 = int(prec1),int(prec2)
    assert prec1 >= def_prec1
    assert prec2 >= def_prec2
    prec1 , prec2 = max(prec1,prec2) , min(prec1,prec2)
    return prec1 , prec2
#--------------------------------------------------#
def prepare_stoc(arguments):
    # Default value
    if len(arguments) == 0: return []
    return arguments
#--------------------------------------------------#
def prepare_hlopt(arguments):
    # User does not want to perform calculations
    mode  = "0"
    if "nocalc" in arguments:
        mode = "1"
        arguments.remove("nocalc")
    # Get selected conformers
    lowerconf , upperconf = 0 , float("inf")
    if len(arguments) == 1:
       lowerconf = int(arguments[0])
       upperconf = int(arguments[0])
    elif len(arguments) == 2:
       lowerconf = int(arguments[0])
       upperconf = int(arguments[1])
    elif len(arguments) != 0: raise Exception
    # Return data
    return mode , min(lowerconf,upperconf) , max(lowerconf,upperconf)
#--------------------------------------------------#
def prepare_ms(arguments):
    mode = "all"
    if   "ll" in arguments: mode = "ll"
    elif "hl" in arguments: mode = "hl"
    elif len(arguments) != 0: raise Exception
    return mode
#==================================================#

#==================================================#
def get_dict_options():
    current = None
    options = {current:[]}
    for arg in sys.argv[1:]:
        if arg in OPTIONS:
           current = arg
           options[current] = []
        else:
            options[current].append(arg)
    return options
#--------------------------------------------------#
def prepare_options(options):
    for option,arguments in options.items():
        try:
           if   option == "--smiles"   : options[option] = prepare_smiles(arguments)
           elif option == "--cartesian": options[option] = prepare_cartesian(arguments)
           elif option == "--input"    : options[option] = prepare_input(arguments)
           elif option == "--prec"     : options[option] = prepare_prec(arguments)
           elif option == "--stoc"     : options[option] = prepare_stoc(arguments)
           elif option == "--hlopt"    : options[option] = prepare_hlopt(arguments)
           elif option == "--msho"     : options[option] = prepare_ms(arguments)
           elif option == "--mstor"    : options[option] = prepare_ms(arguments)
        except:
           pp.print_option_error(option)
           raise exc.END
    return options
#--------------------------------------------------#
def get_options_from_prompt():
    # Get user arguments (options)
    options = get_dict_options()

    # --inp --> --input
    if "--inp" in options:
        options["--input"] = options.get("--input",[]) + options["--inp"]
        options.pop("--inp")

    # Prepare options
    options = prepare_options(options)

    # Remove None from options
    options.pop(None)

    # Return
    return options
#==================================================#

