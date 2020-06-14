import sys
import subprocess
import os
from python_tools.utilities import Utilities

class Correlator:

    def __init__(self,
                 data_filename,
                 data_filename_2,
                 output_filename
                 box_size,
                 corr_type,
                 dim1_min,
                 dim1_max,
                 dim1_nbin,
                 dim2_min,
                 dim2_max,
                 dim2_nbin):

        # declare attributes
        self.output_filename = output_filename
        self.box_size = box_size
        self.corr_type = corr_type

        if os.path.isfile(data_filename):
            self.data_filename = data_filename
        else:
            sys.exit('{} not found.'.format(data_filename))

        if 'CCF' in self.corr_type:
            if os.path.isfile(data_filename_2):
                self.data_filename_2 = data_filename_2
            else:
                sys.exit('{} not found.'.format(data_filename_2))
        
        self.dim1_min = dim1_min
        self.dim2_min = dim2_min
        self.dim1_max = dim1_max
        self.dim2_max = dim2_max
        self.dim1_nbin = dim1_nbin
        self.dim2_nbin = dim2_nbin

        # need to check this
        self.ngrid = int(self.box_size / 15)

        print('Running contrast with the following parameters:\n')
        print('output_filename: {}'.format(self.output_filename))
        print('ngrid: {}'.format(self.ngrid))

        # run the desired correlation function
        getattr(self, corr_type)()



    def CF_monopole(self):
        '''
        Two point autocorrelation function
        in bins of r.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'CF_monopole.exe',
               self.data_filename,
               output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]
        
        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def CF_rmu(self):
        '''
        Two-point autocorrelation function
        in bins of r and mu.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'CCF_rmu.exe',
               self.data_filename,
               output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.dim2_nbin),
               str(self.ngrid)]
        
        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)


    def CCF_monopole(self):
        '''
        Two point cross-correlation function
        in bins of r.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'CCF_monopole.exe',
               self.data_filename,
               self.data_filename_2,
               output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]
        
        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)


    def CCF_rmu(self):
        '''
        Two point cross-correlation function
        in bins of r and mu.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'CCF_rmu.exe',
               self.data_filename,
               self.data_filename_2,
               output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.dim2_nbin),
               str(self.ngrid)]
        
        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def CCF_spi(self):
        '''
        Two point cross-correlation function
        in bins of s (perpendicular) and pi
        (parallel).
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'CCF_spi.exe',
               self.data_filename,
               self.data_filename_2,
               output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]
        
        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

def next_pow_two(n):
    '''
    Returns the largest power of two
    smaller than a given positive integer.
    '''
    i = 1
    while i < n:
        i = i << 1
    return i