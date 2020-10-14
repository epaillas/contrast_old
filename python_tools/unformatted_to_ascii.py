import numpy as np
import sys
from scipy.io import FortranFile
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--input_filename', type=str, required=True)
parser.add_argument('--output_filename', type=str, required=True)
args = parser.parse_args()

print('Converting unformatted Fortran 90 file to ASCII.')
print('input_filename: {}'.format(args.input_filename))
print('output_filename: {}'.format(args.output_filename))

f = FortranFile(args.input_filename, 'r')
nrows = f.read_ints()[0]
ncols = f.read_ints()[0]
data = f.read_reals(dtype=np.float64).reshape(nrows, ncols)
f.close()

np.savetxt(args.output_filename, data)
