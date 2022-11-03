This is the newly developed EOVSA flare pipeline. This is intended to produce some quicklook images at some user defined frequencies. THE IMAGES ARE NOT INTENDED TO BE USED FOR ANY ANALYSIS WHICH REQUIRES ACCURATE FLUX MEASUREMENTS OR ACCURATE POSITION OR SHAPE MEASUREMENTS. The code is written in a manner such that the user can run it with minimal inputs. 

The essential parameters:

1. workpath: The directory where all the work will be done. Please leave rawcal_ms='' . Usecases when it should 
		not be empty is given later. NOTE THAT ALL PATHS PROVIDED IN INPUTS SHOULD END WITH /
		
2. rawcal_ms: ''
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

identify_data_gap: Often EOVSA data shows data gaps, and setting this to 1, identifies those gaps and flags them
doslfcal:          1 means slfcal will be done. 0 means otherwise
doapply:           1 means will do the following action. Apply the caltables, and write the self-calibrated 		    dataset. Note that the MS, if exists, will not be overwritten.
doclean:            1 means will produce the final cleaned images using the self-calibrated dataset.

Sometimes, we need to run the selfcal code in a different machine than where the gen_IDB_MS was run. If you are in a similar situation, you can copy just the MS to the new machine. Then supply the name of the MS in rawcal_MS. Please ensure that it is accessible from workpath .