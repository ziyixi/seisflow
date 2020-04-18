"""
compile_julia.py: compile julia fortran dynamic libraries.
"""
from .. import get_julia
import sh


def compile_julia():
    path = "specfem_gll.jl/src/fortran/compile.jl"
    full_path = get_julia(path)
    sh.julia(full_path)


if __name__ == "__main__":
    compile_julia()
