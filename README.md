This is the newly developed EOVSA flare pipeline. This is intended to produce some quicklook images at some user defined frequencies. THE IMAGES ARE NOT INTENDED TO BE USED FOR ANY ANALYSIS WHICH REQUIRES ACCURATE FLUX MEASUREMENTS OR ACCURATE POSITION OR SHAPE MEASUREMENTS. The code is written in a manner such that the user can run it with minimal inputs. 

The essential parameters:

1. workpath: The directory where all the work will be done. Please leave rawcal_ms='' . Usecases when it should 
		not be empty is given later. NOTE THAT ALL PATHS PROVIDED IN INPUTS SHOULD END WITH /
		
2. rawcal_ms: ''   ### FILL THIS IF YOU ARE CHANGING MACHINES IN BETWEEN CODES. DESCRIBED BELOW.
3. starttime: Provide the starttime of the required data
4. endtime:   Provide the endtime of the required data.

Please ensure that the total duration is ~1 hour. The code automatically dteects the flare, and that requires a baseline. Baseline calculation becomes highly erroneous if the observation duration is very small. Additionally please note that format of the starttime and endtime is strict. The format is 
YYYY-MM-DD HH:MM:SS

5. selfcal_spws:  Give a list of spws where you want the images
6. total_duration: The total seconds to image. The center will be at the location of the flare
7. final_image_cadence: The cadence of the final images
8. final_image_int:     Integration time of final images

First please run gen_IDB_MS.py after going inside a CASA environment. This code is tested for CASA 5.4 only. 

Next please come out of CASA. The next code is run in a normal python environment and uses the modular CASA 6. Please run IDB_selfcal_pipeline_version.py from the command line or inside Ipython.

Advanced parameters:

1. identify_data_gap: Often EOVSA data shows data gaps, and setting this to 1, identifies those gaps and flags them
2. doslfcal:          1 means slfcal will be done. 0 means otherwise
3. doapply:           1 means will do the following action. Apply the caltables, and write the self-calibrated dataset. Note that the MS, if exists, will not be overwritten.
4. doclean:            1 means will produce the final cleaned images using the self-calibrated dataset.
5. cell:               Setting calc_cell  to False, will mean that you want to provide cell yourself. Note that you need to give cell for each spw in the MS in the form of a list.


THINGS to NOTE:

1. Sometimes, we need to run the selfcal code in a different machine than where the gen_IDB_MS was run. If you are in a similar situation, you can copy just the MS to the new machine. Then supply the name of the MS in rawcal_MS. Please ensure that it is accessible from workpath .

2. Always ensure that you have a working installation of SUNCASA

3. Please ensure that EOVSA is added to the observatories list of CASA. If not, please consider adding it. As a hack, you can run listobs in the pipeline machine and write the ouput in a file. If the rawcal_ms="temp.ms", then the output of listobs should be dumped into "temp.listobs". Copy the output of listobs to the workpath where you want to run "IDB_selfcal_pipeline_version.py" .


