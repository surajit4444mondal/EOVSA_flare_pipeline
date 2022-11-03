

# =========== task handlers =============
cal_disk = 0  ## apply calibration tables from full disc imaging
fullday_cal_loc=''## Location of the full day calibration tables
identify_data_gap=1  ### identify data gaps
doslfcal = 1  # main cycle of doing selfcalibration
doapply = 1 # apply the results
doclean=1  ### produce final image
selfcal_spws=[3,8,15,24,35,42]  ### selfcal \
                        ### and final imaging would be only in these spws 

##===================== final imaging parameters ===========================
total_duration=240 ### seconds to image; time will be centred on detected flare peak
final_image_cadence=10 ### seconds
final_image_int=10 ### integration time of final images
min_restoring_beam=6 ## arcsec
beam_1GHz='89.7arcsec'### in arcseconds
cell_size=2 ## arcsec
final_imsize=512 ### pixel

# ============ declaring the working directories ============
### remember / is necessary in all the folder names

workpath = '/home/surajit/Downloads/flare_pipeline_results/20220930/' ## / is needed at the end of all paths
slfcaldir = workpath+ 'slfcal_v5/'  # place to put all selfcalibration products
imagedir = slfcaldir + 'images/'  # place to put all selfcalibration images
caltbdir = slfcaldir+'caltables/'  # place to put calibration tables
rawcal_ms='IDB20220930_1600-1700.ms'

# ============= time to image =================
starttime='2022-09-30 16:00:00'   ### has strict formating rules
endtime='2022-09-30 17:00:00'


# ============ selfcal parameters ===============

refantenna = '0'
calc_cell=True ### If set to False use the value in beam given below
cell=[10]  ### size needs to be same as the number of spws
calc_imsize=True   ### is False uses the value given below


max_frac_freq_avg=0.5  ### I will average at most this much fractional bandwidth

maxiter=10  ### maximum selfcal iterations
uvlim=25
avg_spw_max=5
flag_antennas = '' ###anything except 13~15. Those antennas are always flagged. 
phasecenter=''

# ========== end of input parameters =================
