"""
compile_fortran.py: compile the fortran code in julia before running the script.
"""
import sh

from ... import get_julia


def compile_fortran():
    """
    check the compiling status.
    """
    # * firstly we get the path of the compile.jl
    compile_path = get_julia("specfem_gll.jl/src/fortran/compile.jl")
    sh.julia(compile_path)
