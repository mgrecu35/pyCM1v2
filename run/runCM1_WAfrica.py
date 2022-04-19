import pyCM1 as cm1
import sys, io
import matplotlib.pyplot as plt
import numpy as np
from contextlib import redirect_stdout

f = io.StringIO()

#with redirect_stdout(f):
ib,ie,jb,je,kb,ke,\
    ibm,iem,jbm,jem,kbm,kem,numq,tq=cm1.pycm1_init()
from netCDF4 import Dataset

stop
