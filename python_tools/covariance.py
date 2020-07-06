import numpy as np
import sys
import glob
import click
from utilities import Utilities
from scipy.integrate import quad, simps
from scipy.interpolate import InterpolatedUnivariateSpline


def CovarianceMatrix(data, norm=False):
    """
    Assumes rows are observations,
    columns are variables
    """
    nobs, nbins = np.shape(data)
    mean = np.mean(data, axis=0)
    cov = np.zeros([nbins, nbins])

    for k in range(nobs):
        for i in range(nbins):
            for j in range(nbins):
                cov[i, j] += (data[k, i] - mean[i])*(data[k, j] - mean[j])

    cov /= nobs - 1
    
    if norm:
        corr = np.zeros_like(cov)
        for i in range(nbins):
            for j in range(nbins):
                corr[i, j] = cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])
        return corr
    else:
        return cov

def CrossCovarianceMatrix(data1, data2, norm=False):
    """
    Assumes rows are observations,
    columns are variables
    """
    nobs, nbins = np.shape(data1)
    mean1 = np.mean(data1, axis=0)
    mean2 = np.mean(data2, axis=0)
    cov = np.zeros([nbins, nbins])

    for k in range(nobs):
        for i in range(nbins):
            for j in range(nbins):
                cov[i, j] += (data1[k, i] - mean1[i])*(data2[k, j] - mean2[j])

    cov /= nobs - 1
    
    if norm:
        corr = np.zeros_like(cov)
        for i in range(nbins):
            for j in range(nbins):
                corr[i, j] = cov[i, j] / np.sqrt(cov[i, i] * cov[j, j])
        return corr
    else:
        return cov

def MultipoleCovariance(handle_mocks, smin, smax, multipoles):
        smin = float(smin)
        smax = float(smax)
        files_mocks = sorted(glob.glob(handle_mocks))
        mock_datavec = []
        nmocks = len(files_mocks)

        print('nmocks: {}'.format(nmocks))

        for fname in files_mocks:
            s, mu, xi_smu_mock = Utilities.ReadData_TwoDims(fname)
            s, xi0 = Utilities.getMultipole(0, s, mu, xi_smu_mock)
            s, xi2 = Utilities.getMultipole(2, s, mu, xi_smu_mock)
            s, xi4 = Utilities.getMultipole(4, s, mu, xi_smu_mock)

            # only keep scales to fit later
            idx = (s >= smin) & (s <= smax)
            s = s[idx]
            xi0 = xi0[idx]
            xi2 = xi2[idx]
            xi4 = xi4[idx]

            if multipoles == '0+2':
                datavec = np.concatenate((xi0, xi2))
            elif multipoles == '0+2+4':
                datavec = np.concatenate((xi0, xi2, xi4))
            elif multipoles == '0':
                datavec = xi0
            elif multipoles == '2':
                datavec = xi2
            elif multipoles == '4':
                datavec = xi4
            else:
                sys.exit('Fit type not recognized.')

            mock_datavec.append(datavec)

        mock_datavec = np.asarray(mock_datavec)
        cov = CovarianceMatrix(mock_datavec)
        print('np.shape(cov): {}'.format(np.shape(cov)))

        return cov

def JointMultipoleCovariance(handle_mocks, smins, smaxs, multipoles):
        handle_mocks = handle_mocks.split(',')
        smins = [float(i) for i in smins.split(',')]
        smaxs = [float(i) for i in smaxs.split(',')]
        ndenbins = len(handle_mocks)

        files_mocks = {}
        smin = {}
        smax = {}
        for j in range(ndenbins):
            denbin = 'den{}'.format(j)
            files_mocks[denbin] = sorted(glob.glob(handle_mocks[j]))
            smin[denbin] = smins[j]
            smax[denbin] = smaxs[j]

        nmocks = len(files_mocks['den0'])
        mock_datavec = []

        print('nmocks: {}'.format(nmocks))
        print('ndenbins: {}'.format(ndenbins))

        for i in range(nmocks):
            datavec = np.array([])
            for j in range(ndenbins):
                denbin = 'den{}'.format(j)
                fname = files_mocks[denbin][i]
                s, mu, xi_smu_mock = Utilities.ReadData_TwoDims(fname)
                s, xi0 = Utilities.getMultipole(0, s, mu, xi_smu_mock)
                s, xi2 = Utilities.getMultipole(2, s, mu, xi_smu_mock)
                s, xi4 = Utilities.getMultipole(4, s, mu, xi_smu_mock)

                # only keep scales to fit later
                fitscale = (s >= smin[denbin]) & (s <= smax[denbin])
                s = s[fitscale]
                xi0 = xi0[fitscale]
                xi2 = xi2[fitscale]
                xi4 = xi4[fitscale]

                if multipoles == '0+2':
                    datavec = np.concatenate((datavec, xi0, xi2))
                elif multipoles == '0+2+4':
                    datavec = np.concatenate((datavec, xi0, xi2, xi4))
                elif multipoles == '0':
                    datavec = np.concatenate((datavec, xi0))
                elif multipoles == '2':
                    datavec = np.concatenate((datavec, xi2))
                elif multipoles == '4':
                    datavec = np.concatenate((datavec, xi4))
                else:
                    sys.exit('Fit type not recognized.')

            mock_datavec.append(datavec)

        mock_datavec = np.asarray(mock_datavec)
        cov = CovarianceMatrix(mock_datavec)
        print('np.shape(cov): {}'.format(np.shape(cov)))
        
        return cov

@click.command()
@click.option('--handle_in', type=str, required=True, help='Handle from mocks')
@click.option('--handle_out', type=str, required=True, help='Handle for the mean')
@click.option('--smin', type=str, default=0.0, help='Minimum scale to fit (in Mpc/h)')
@click.option('--smax', type=str, default=100.0, help='Maximum scale to fit (in Mpc/h)')
@click.option('--multipoles', type=str, default='0+2', help='Which multipoles to fit?')

def get_covariance(handle_in,
                   handle_out,
                   smin,
                   smax,
                   multipoles):

    # figure out if single or joint fit
    ndenbins = len(handle_in.split(','))

    if ndenbins > 1:
        print('Calculating covariance for joint fit.)
        cov = JointMultipoleCovariance(handle_mocks=handle_in,
                                smins=smin,
                                smaxs=smax,
                                multipoles=multipoles)
    elif ndenbins == 1:
        print('Calculating covariance for single fit.')
        cov = MultipoleCovariance(handle_mocks=handle_in,
                              smin=smin,
                              smax=smax,
                              multipoles=multipoles)
    else:
        sys.exit('Could not figure out number of density bins. Check handle.')

    np.save(handle_out, cov)

if __name__ == '__main__':
    get_covariance()