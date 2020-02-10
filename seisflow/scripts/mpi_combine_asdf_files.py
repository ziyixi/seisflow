"""
mpi_combine_asdf_files.py: combine asdf files massively in two directories.
"""
from ..asdf.combine_asdf_files import combine_asdf
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
