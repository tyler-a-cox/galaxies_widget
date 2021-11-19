import numpy as np
import pylab as plt
import shutil
import tempfile
from ftplib import FTP
import subprocess
import os, os.path, sys
from astropy.io import ascii

from scipy.optimize import curve_fit

from scipy.ndimage import gaussian_filter1d

_ERASESTR = ""
_CACHE_BASEDIR = os.path.join(os.getenv("HOME"), "cache")
_CACHE_VIZIER_DIR = os.path.join(_CACHE_BASEDIR, "vizier")
_ADIBEKYAN_VIZIER_NAME = "J/A+A/545/A32"


def vizier(cat, filePath, ReadMePath, catalogname="catalog.dat", readmename="ReadMe"):
    """
    """
    _download_file_vizier(cat, filePath, catalogname=catalogname)
    _download_file_vizier(cat, ReadMePath, catalogname=readmename)
    catfilename = os.path.basename(catalogname)
    with open(ReadMePath, "r") as readmefile:
        fullreadme = "".join(readmefile.readlines())
        if catfilename.endswith(".gz") and not catfilename in fullreadme:
            # Need to gunzip the catalog
            try:
                subprocess.check_call(["gunzip", filePath])
            except subprocess.CalledProcessError as e:
                print("Could not unzip catalog %s" % filePath)
                raise
    return None


def _download_file_vizier(cat, filePath, catalogname="catalog.dat"):
    sys.stdout.write(
        "\r" + "Downloading file %s from Vizier ...\r" % (os.path.basename(filePath))
    )
    sys.stdout.flush()
    try:
        # make all intermediate directories
        os.makedirs(os.path.dirname(filePath))
    except OSError:
        pass
    # Safe way of downloading
    downloading = True
    interrupted = False
    file, tmp_savefilename = tempfile.mkstemp()
    os.close(file)  # Easier this way
    ntries = 1
    while downloading:
        try:
            ftp = FTP("cdsarc.u-strasbg.fr")
            ftp.login()
            ftp.cwd(os.path.join("pub", "cats", cat))
            with open(tmp_savefilename, "wb") as savefile:
                ftp.retrbinary("RETR %s" % catalogname, savefile.write)
            shutil.move(tmp_savefilename, filePath)
            downloading = False
            if interrupted:
                raise KeyboardInterrupt
        except:
            raise
            if not downloading:  # Assume KeyboardInterrupt
                raise
            elif ntries > _MAX_NTRIES:
                raise IOError(
                    "File %s does not appear to exist on the server ..."
                    % (os.path.basename(filePath))
                )
        finally:
            if os.path.exists(tmp_savefilename):
                os.remove(tmp_savefilename)
        ntries += 1
    sys.stdout.write("\r" + _ERASESTR + "\r")
    sys.stdout.flush()
    return None


def read_sn_abu(verbose=False):
    """
    """
    # Generate file path and name
    tPath = os.path.join(_CACHE_VIZIER_DIR, "cats", *_ADIBEKYAN_VIZIER_NAME.split("/"))
    filePath = os.path.join(tPath, "table4.dat")
    readmePath = os.path.join(tPath, "ReadMe")
    # download the file
    if not os.path.exists(filePath):
        vizier(
            _ADIBEKYAN_VIZIER_NAME,
            filePath,
            readmePath,
            catalogname="table4.dat",
            readmename="ReadMe",
        )
    # Read with astropy, w/o the .gz
    table = ascii.read(filePath, readme=readmePath, format="cds")
    return table.to_pandas()


def get_sn_xhdist(element, plotrange=[-3, 1], bins=31):
    """
    """
    sndata = read_sn_abu()
    c, b = np.histogram(
        sndata["[{}/H]".format(element.capitalize())], bins=bins, density=True
    )
    b = (b[1:] + b[:-1]) / 2
    return c, b


Zsolar = 0.014


def fgas(t, mode="", tau=11e9):
    """
    Gas fraction
    """
    if mode == "constant":
        C = 0.7 / 13e9
        fraction = 1 - C * t
        return fraction

    else:
        C = 0.7 / (tau - tau * np.exp(-13e9 / tau))
        fraction = 1 - C * (tau - tau * np.exp(-t / tau))
        return fraction


def calc_metallicity(t, mode="constant", tau=11e9, z_yield=0.5):
    """
    """
    return -z_yield * np.log(fgas(t, mode=mode, tau=tau))


def metallicity_distribution(Z, z_yield=0.01):
    """
    """
    Z *= Zsolar
    return Z / z_yield * np.exp(-Z / z_yield)


def fit_closed(logZ, z_yield=0.01):
    """
    """
    dN = metallicity_distribution(10 ** logZ, z_yield=z_yield)
    return dN / np.trapz(dN, logZ)


def calc_metallicity_leaky(t, mode="constant", tau=11e9, z_yield=0.01, eta=3):
    """
    """
    return -z_yield / (1 + eta) * np.log(fgas(t, mode=mode, tau=tau))


def metallicity_distribution_leaky(Z, z_yield=0.01, eta=3, r=0.2):
    """
    """
    Z *= Zsolar
    pp = z_yield / (1 + eta - r)
    return Z / pp * np.exp(-Z / pp)


def fit_leaky(logZ, z_yield=0.01, eta=3):
    """
    """
    dN = metallicity_distribution_leaky(10 ** logZ, z_yield=z_yield, eta=eta)
    return dN / np.trapz(dN, logZ)


def calc_metallicity_accreting(t, mode="constant", tau=11e9, z_yield=0.01, eta=3):
    """
    """
    return calc_metallicity_accreting_leaky(
        t, mode=mode, tau=tau, z_yield=z_yield, eta=eta, r=0
    )


def calc_metallicity_accreting_leaky(
    t, mode="constant", tau=11e9, z_yield=0.01, eta=3, r=0.4
):
    """
    """
    return (
        z_yield
        / (1 + eta - r)
        * (1 - np.exp((1 + eta - r) / (1 - r) * (1 - 1 / fgas(t, mode=mode, tau=tau))))
    )


def dNdZ_accreting(Z):
    """
    """


def metallicity_distribution_accreting(Z, z_yield=0.01):
    """
    """
    return Z / z_yield * np.exp(-Z / z_yield)


def metallicity_distribution_accreting_leaky(Z, z_yield=0.01, eta=3, r=0.5):
    """
    """
    Z *= Zsolar
    pp = z_yield / (1 + eta - r)
    m = Z / (pp - Z)
    m[Z > pp] = 0
    return m


def fit_accreting(logZ, z_yield=0.01, eta=3):
    """
    """
    dN = metallicity_distribution_accreting_leaky(
        10 ** logZ, z_yield=z_yield, eta=eta, r=0
    )
    k = gaussian_filter1d(np.abs(dN), 3)
    dk = k / np.trapz(k, logZ)
    return dk


def fit_accreting_leaky(logZ, z_yield=0.01, eta=3, r=0.5):
    """
    """
    dN = metallicity_distribution_accreting_leaky(
        10 ** logZ, z_yield=z_yield, eta=eta, r=r
    )
    k = gaussian_filter1d(np.abs(dN), 3)
    dk = k / np.trapz(k, logZ)
    return dk
