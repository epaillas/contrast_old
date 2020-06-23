import numpy as np
import argparse

parser = argparse.ArgumentParser(description='')

parser.add_argument('--data_filename', type=str)
args = parser.parse_args()  

fmt = 2*'%15.5f '

# infer handle from data filename
handle = args.data_filename.split('monopole')[0]

# read galaxy monopole data
data = np.genfromtxt(args.data_filename)

# xi_r
r_for_xi = data[:,0]
delta = data[:,1]
ext_out = 'xi_r'
fout = handle + ext_out
cout = np.c_[r_for_xi, delta]
np.savetxt(fout, cout, fmt=fmt)

# integrated xi_r
r_for_xi = data[:,0]
int_delta = data[:,2]
ext_out = 'int_xi_r'
fout = handle + ext_out
cout = np.c_[r_for_xi, int_delta]
np.savetxt(fout, cout, fmt=fmt)

# v_r
r_for_v = data[:,0]
v_r = data[:,3]
ext_out = 'v_r'
fout = handle + ext_out
cout = np.c_[r_for_v, v_r]
np.savetxt(fout, cout, fmt=fmt)

# sv_r
r_for_v = data[:,0]
sv_r = data[:,4]
ext_out = 'sv_r'
fout = handle + ext_out
cout = np.c_[r_for_v, sv_r]
np.savetxt(fout, cout, fmt=fmt)

# v_los
r_for_v = data[:,0]
v_los = data[:,5]
ext_out = 'v_los'
fout = handle + ext_out
cout = np.c_[r_for_v, v_los]
np.savetxt(fout, cout, fmt=fmt)


# sv_los
r_for_v = data[:,0]
sv_los = data[:,6]
ext_out = 'sv_los'
fout = handle + ext_out
cout = np.c_[r_for_v, sv_los]
np.savetxt(fout, cout, fmt=fmt)





