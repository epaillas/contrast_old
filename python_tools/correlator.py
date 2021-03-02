import sys
import subprocess
import os
from python_tools.utilities import Utilities


class Correlator:

    def __init__(self,
                 data_filename,
                 data_filename_2,
                 output_filename,
                 box_size,
                 corr_type,
                 dim1_min,
                 dim1_max,
                 dim1_nbin,
                 dim2_min,
                 dim2_max,
                 dim2_nbin,
                 ngrid,
                 qperp,
                 qpara):

        # declare attributes
        self.output_filename = output_filename
        self.box_size = box_size
        self.corr_type = corr_type

        if os.path.isfile(data_filename):
            self.data_filename = data_filename
        else:
            sys.exit('{} not found.'.format(data_filename))

        if data_filename_2 is None:
            self.data_filename_2 = self.data_filename
        else:
            self.data_filename_2 = data_filename_2

        self.dim1_min = dim1_min
        self.dim2_min = dim2_min
        self.dim1_max = dim1_max
        self.dim2_max = dim2_max
        self.dim1_nbin = dim1_nbin
        self.dim2_nbin = dim2_nbin
        self.ngrid = ngrid
        self.qperp = qperp
        self.qpara = qpara

        print('Running contrast with the following parameters:\n')
        print('corr_type: {}'.format(self.corr_type))
        print('box_size: {}'.format(self.box_size))
        print('data_filename: {}'.format(self.data_filename))
        if self.data_filename_2 != self.data_filename:
            print('data_filename_2: {}'.format(self.data_filename_2))
        print('output_filename: {}'.format(self.output_filename))
        print('dim1_min: {}'.format(self.dim1_min))
        print('dim1_max: {}'.format(self.dim1_max))
        print('dim2_min: {}'.format(self.dim2_min))
        print('dim2_max: {}'.format(self.dim2_max))
        print('dim1_nbin: {}'.format(self.dim1_nbin))
        print('dim2_nbin: {}'.format(self.dim2_nbin))
        print('ngrid: {}'.format(self.ngrid))

        # run the desired correlation function
        getattr(self, corr_type)()

    def tpcf(self):
        '''
        Two point autocorrelation function
        in bins of r.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'tpcf.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid),
               str(self.qperp),
               str(self.qpara)]

        log = open(log_filename, 'w+')
        subprocess.run(cmd, stdout=log, stderr=log)

    def s_mu_tpcf(self):
        '''
        Two-point autocorrelation function
        in bins of r and mu.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 's_mu_tpcf_custombins.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.dim2_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def s_pi_tpcf(self):
        '''
        Two point cross-correlation function
        in bins of s (perpendicular) and pi
        (parallel).
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 's_pi_tpcf.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def mean_radial_velocity_vs_r(self):
        '''
        Line-of-sight mean pairwise velocity
        as a function of r and mu.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'mean_radial_velocity_vs_r.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def mean_velocity_sphere(self):
        '''
        Average velocity of particles within
        a sphere for each input centre.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'mean_velocity_sphere.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def mean_transverse_velocity_vs_r(self):
        '''
        Line-of-sight mean pairwise velocity
        as a function of r and mu.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'mean_transverse_velocity_vs_r.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def std_radial_velocity_vs_r(self):
        '''
        Radial pairwise velocity dispersion
        as a function of r.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'std_radial_velocity_vs_r.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def std_transverse_velocity_vs_r(self):
        '''
        Transverse pairwise velocity dispersion
        as a function of r.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'std_transverse_velocity_vs_r.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def std_los_velocity_vs_rmu(self):
        '''
        Line-of-sight pairwise velocity dispersion
        as a function of r and mu.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'std_los_velocity_vs_rmu.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.dim2_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.call(cmd, stdout=log, stderr=log)

    def omp_tpcf(self):
        '''
        Two point autocorrelation function
        in bins of r.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'omp_tpcf.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.run(cmd, stdout=log, stderr=log)

    def neighbour_search(self):
        '''
        Search for the number of neighbours at 
        a given distance.
        '''
        log_filename = self.output_filename + '.log'

        binpath = sys.path[0] + '/bin/'
        cmd = [binpath + 'neighbour_search.exe',
               self.data_filename,
               self.data_filename_2,
               self.output_filename,
               str(self.box_size),
               str(self.dim1_min),
               str(self.dim1_max),
               str(self.dim1_nbin),
               str(self.ngrid)]

        log = open(log_filename, 'w+')
        subprocess.run(cmd, stdout=log, stderr=log)
