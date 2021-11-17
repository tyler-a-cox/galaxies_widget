import numpy as np
import pylab as plt
import shutil
import tempfile
from ftplib import FTP
import subprocess
import os, os.path, sys
from astropy.io import ascii

_ERASESTR = ""
_CACHE_BASEDIR = os.path.join(os.getenv("HOME"), ".galdyncourse", "cache")
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
