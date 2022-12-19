from suncasa.tasks import task_calibeovsa as calibeovsa
from suncasa.tasks import task_importeovsa as timporteovsa
from split_cli import split_cli as split
from util import Time
import dump_tsys
import numpy as np
import os
import glob
from eovsapy import util
import datetime as dt
import subprocess

os.system("rm -rf inputs.pyc")
from inputs import *

os.chdir(workpath)
workdir = os.getcwd() + '/'  # main working directory. Using current directory in this example

dirs_ = [workdir, slfcaldir, imagedir, caltbdir]
for d in dirs_:
    if not os.path.exists(d):
        os.makedirs(d)

if not os.path.exists(workpath):
    os.makedirs(workpath)

os.chdir(workpath)

if endtime>starttime:
    trange = Time([starttime, endtime])  ## M4.4
else:
    raise RuntimeError("Ordering ot times wrong")

namesuffix = '_' + trange[0].to_datetime().strftime('%H%M') + '-' + trange[1].to_datetime().strftime('%H%M')
idbdir = util.get_idbdir(trange[0])
info = dump_tsys.rd_fdb(Time(starttime.split(' ')[0]))
idxs, = np.where(info['SOURCEID'] == '')
for k, v in info.iteritems():
    info[k] = info[k][~(info[k] == '')]
sidx = np.where(
    np.logical_and(info['SOURCEID'] == 'Sun', info['PROJECTID'] == 'NormalObserving') & np.logical_and(
        info['ST_TS'].astype(np.float) >= trange[0].lv,
        info['ST_TS'].astype(np.float) <= trange[1].lv))

filelist = info['FILE'][sidx]

outpath = 'msdata/'
if not os.path.exists(outpath):
    os.makedirs(outpath)
inpath = idbdir + '{}/'.format(trange[0].datetime.strftime("%Y%m%d"))
ncpu = 1

print ("Starting to import IDB MS")
msfiles = timporteovsa.importeovsa(idbfiles=[inpath + ll for ll in filelist], ncpu=ncpu, timebin="0s", width=1,\
                                 visprefix=outpath, nocreatms=False, doconcat=False,\
                              modelms="", doscaling=False, keep_nsclms=False, udb_corr=True,\
                            use_exist_udbcorr=True)
print ("IDB MS imported")
msfiles = [outpath + ll + '.ms' for ll in filelist]
print(msfiles)
concatvis = os.path.basename(msfiles[0])[:11] + namesuffix + '.ms'
print ("Running calibeovsa, concatting to ", concatvis)
vis = calibeovsa.calibeovsa(msfiles, caltype=['refpha', 'phacal'], interp='nearest', doflag=True,\
			 flagant='13~15', doimage=False, doconcat=True,\
	                    concatvis=outpath + concatvis, keep_orig_ms=True)
print("Initial calibration done")


