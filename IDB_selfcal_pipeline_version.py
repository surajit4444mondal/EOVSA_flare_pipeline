import numpy as np
import os
import glob
import time
from astropy.io import fits
import datetime as dt
import subprocess
from casatasks import *
from casatools import image,ms,msmetadata,table
from suncasa.utils import qlookplot
from suncasa.utils import helioimage2fits as hf
import timeit
from sunpy.time import parse_time
import matplotlib
matplotlib.use('Agg')

start=timeit.default_timer()
ia=image()
ms=ms()
msmd=msmetadata()
tb=table()

os.system("rm -rf imports.pyc")
from inputs import *

maximum_spw=max(selfcal_spws)
minimum_spw=min(selfcal_spws)
selfcal_spw='0~'+str(maximum_spw)

def find_sidelobe_level(image):
	head=fits.getheader(image)
	imsize=head['naxis1']
	cellsize=float(abs(head['CDELT1']))*3600
	
	psf_image=image[:-5]+".psf"
	ia.open(psf_image)
	psf_data=ia.getchunk()[:,:,0,0]
	ia.close()


	x0=head['CRPIX2']-1
	y0=head['CRPIX1']-1
			
	posang=head['BPA']*np.pi/180 # in degrees --> radian
	bmaj=head['BMAJ']*3600/cellsize # in degrees --> pix Full maj axis extent
	bmin=head['BMIN']*3600/cellsize # ''
	R=np.array([[np.cos(posang),-np.sin(posang)],[np.sin(posang),np.cos(posang)]])
	a=bmaj/2. # Chosing nsig region for calculating flux density
	b=bmin/2.
	
	low_limit_x=int(x0-a)
	low_limit_y=int(y0-a)
	upper_limit_x=int(x0+a)
	upper_limit_y=int(y0+a)
	psf_region_values=[]
			
	for k in range(low_limit_y,upper_limit_y):
		for j in range(low_limit_x,upper_limit_x):
			x_v=np.matrix([[k-x0],[j-y0]])	
			R=np.matrix(R)
			xp_v=R*x_v
			xp_v=np.array(xp_v.transpose()).reshape(2)
			if xp_v[0]**2/a**2 + xp_v[1]**2/b**2<1:
				psf_region_values.append(psf_data[k,j])
				psf_data[k,j]=np.nan
	max_val=np.nanmax(psf_data)
	min_lev=np.nanmin(np.array(psf_region_values))
	return max_val,min_lev	
						
def flag_data_gap(msname,sp):
	ms.open(msname)
	ms.selectinit(datadescid=0,reset=True)
	ms.selectinit(datadescid=sp)
	data=ms.getdata('amplitude')['amplitude']
	ms.close()
	ms.open(msname,nomodify=False)
	ms.selectinit(datadescid=0,reset=True)
	ms.selectinit(datadescid=sp)
	flag=ms.getdata('flag')
	pos=np.where(data<1e-3)
	flag['flag'][pos]=True
	ms.putdata(flag)
	ms.close()
	return

def check_shift(image,shift,cell):
	header=imhead(image)
	major=header['restoringbeam']['major']['value']
	minor=header['restoringbeam']['minor']['value']
	pa=header['restoringbeam']['positionangle']['value']*np.pi/180
	v2=np.array([np.sin(pa),np.cos(pa)])
	s=np.dot(v2,shift)/(np.sqrt(v2[0]**2+v2[1]**2)*np.sqrt(shift[0]**2+shift[1]**2))
	if shift[0]**2+shift[1]**2<1e-3:
		return True
	theta=np.arccos(s)
	abs_shift=np.sqrt(shift[0]**2+shift[1]**2)
	major_shift=abs(abs_shift*np.cos(theta))
	minor_shift=abs(abs_shift*np.cos(theta))
	print (major_shift, minor_shift, major/cell, minor/cell)
	if major_shift>0.75*major/cell:
		return False
	if minor_shift>0.75*minor/cell:
		return False
	return True
	
def grow_mask(image,mask,thres):
	image_data=fits.getdata(image)[0,0,:,:]
	ia.open(mask)
	mask_data=ia.getchunk()[:,:,0,0].T
	ia.close()
	
	shape=np.shape(mask_data)
	rows=shape[0]
	cols=shape[1]
	
	#print (np.shape(image_data))
	
	pos=np.where(mask_data==1)
	if len(pos)==0:
		print("Mask blank. First run gen_mask. Exiting")
		return False
	
	
	xpos1=pos[1]
	ypos1=pos[0]
	
	while np.size(xpos1)!=0:
		for x,y in zip(xpos1,ypos1):
			j=x
			i=y+1
			while i<rows:
				if mask_data[i,j]==1 or mask_data[i,j]==2 or image_data[i,j]<thres:
					break
				else:
					mask_data[i,j]=2
					i+=1
			i=y
			j=x+1
			while j<cols:
				if mask_data[i,j]==1 or mask_data[i,j]==2 or image_data[i,j]<thres:
					break
				else:
					mask_data[i,j]=2
					j+=1
			
		del xpos1
		del ypos1
		del pos
		pos=np.where(mask_data==2)
		xpos1=pos[1]
		ypos1=pos[0]
		mask_data[pos]=1
	
	pos=np.where(mask_data==1)	
	xpos1=pos[1]
	ypos1=pos[0]
	
	while np.size(xpos1)!=0:
		for x,y in zip(xpos1,ypos1):
			j=x
			i=y-1
			while i<rows:
				if mask_data[i,j]==1 or mask_data[i,j]==2 or image_data[i,j]<thres:
					break
				else:
					mask_data[i,j]=2
					i-=1
			i=y
			j=x-1
			while j<cols:
				if mask_data[i,j]==1 or mask_data[i,j]==2 or image_data[i,j]<thres:
					break
				else:
					mask_data[i,j]=2
					j-=1
			
		del xpos1
		del ypos1
		del pos
		pos=np.where(mask_data==2)
		xpos1=pos[1]
		ypos1=pos[0]
		mask_data[pos]=1
	
	pos=np.where(mask_data==1)	
	xpos1=pos[1]
	ypos1=pos[0]
	
	while np.size(xpos1)!=0:
		for x,y in zip(xpos1,ypos1):
			j=x
			i=y-1
			while i<rows:
				if mask_data[i,j]==1 or mask_data[i,j]==2 or image_data[i,j]<thres:
					break
				else:
					mask_data[i,j]=2
					i-=1
			i=y
			j=x-1
			while j<cols:
				if mask_data[i,j]==1 or mask_data[i,j]==2 or image_data[i,j]<thres:
					break
				else:
					mask_data[i,j]=2
					j+=1
			
		del xpos1
		del ypos1
		del pos
		pos=np.where(mask_data==2)
		xpos1=pos[1]
		ypos1=pos[0]
		mask_data[pos]=1
		
	pos=np.where(mask_data==1)
	xpos1=pos[1]
	ypos1=pos[0]
	
	while np.size(xpos1)!=0:
		for x,y in zip(xpos1,ypos1):
			j=x
			i=y+1
			while i<rows:
				if mask_data[i,j]==1 or mask_data[i,j]==2 or image_data[i,j]<thres:
					break
				else:
					mask_data[i,j]=2
					i+=1
			i=y
			j=x-1
			while j<cols:
				if mask_data[i,j]==1 or mask_data[i,j]==2 or image_data[i,j]<thres:
					break
				else:
					mask_data[i,j]=2
					j-=1
			
		del xpos1
		del ypos1
		del pos
		pos=np.where(mask_data==2)
		xpos1=pos[1]
		ypos1=pos[0]
		mask_data[pos]=1
		
	pos=np.where(mask_data==1)
	xpos1=pos[1]
	ypos1=pos[0]
	
	
	del mask_data
	ia.open(mask)
	mask_data=ia.getchunk()
	
	for x,y in zip(pos[1],pos[0]): 
		mask_data[x,y,0,0]=1
	ia.putchunk(mask_data)
	ia.close()
	
	return True

def get_selfcal_times(msname,nspw, starttime):
    times=[]
    flare=[]
    flare_peak=[]
    start_time=[]
    end_time=[]
    start_datetime_obj=[]
    end_datetime_obj=[]

    if os.path.isfile(msname[:-3]+".listobs")==False:
        listobs(msname, listfile=msname[:-3]+".listobs")
    
    t=subprocess.check_output("grep Observed "+msname[:-3]+".listobs",shell=True).decode('utf8').strip().split(' ')[4]
    t=parse_time(t.replace('/',' ')).value
    hmy=t.split('T')[1].split(':')
    hour=int(hmy[0])
    minute=int(hmy[1])
    second=int(float(hmy[2]))

    ymd=starttime.split(' ')[0].split('-')
    year=int(ymd[0])
    month=int(ymd[1])
    day=int(ymd[2])

    start=dt.datetime(year,month,day,hour,minute,second)

    peak=0
    peak_time=0
    
    ms.open(msname)
    for sp in range(minimum_spw,maximum_spw+1):
        ms.selectinit(datadescid=0,reset=True)
        ms.selectinit(datadescid=sp)
        ms.select({'Antenna1':[0,1]})
        data=ms.getdata(['amplitude','axis_info'],ifraxis=True)
        startmjd=data['axis_info']['time_axis']['MJDseconds'][0]
        endmjd=data['axis_info']['time_axis']['MJDseconds'][-1]
        power=data['amplitude'][0,3,0,:]
        #print (data['axis_info']['time_axis'].keys())
        print(sp, np.shape(data['amplitude'])[1])
        median_power=np.median(power)
        smoothed=np.convolve(power, np.ones(10),mode='same')*0.1
        peak_pow=np.max(smoothed)
        mad=np.median(abs(power-median_power))
        thres=5.0
        y=(smoothed-median_power)/mad
        pos=np.argmax(y)
        tot_times=np.size(y)
        if y[pos]>peak:
                peak=y[pos]
                peak_time=data['axis_info']['time_axis']['MJDseconds'][pos]
        del pos
        pos=np.where(y>thres)[0]
        duration=60
        if len(pos)==0:
            flare.append(False)
        else:
            flare.append(True)
            temp=np.convolve(pos,np.array([1,-1]),mode='same')
            #### Couting the number of continuous zeros in temp. This will give the event duration
            counts=[]
            count=1
            for i in temp:
                if i==1:
                    count+=1
                else:
                    if count!=1:
                        counts.append(count)
                        count=0
            counts.append(count)
            
            counts.sort()
            duration=counts[-1]
            
            
            #### max integration time=60s#####
            if duration/3>=60:
                duration=60
            elif duration/3>20 and duration/3<60:
                duration=40
            else:
                duration=20

        peak_pos=np.argmax(smoothed-median_power)
        flare_peak.append(smoothed[peak_pos]-median_power)
        #flare_peak.append(np.max(smoothed-median_power))
        ### assume that the duration is maximum when the peak is highest
        #peak_pos=np.where(abs(smoothed-flare_peak[-1]-median_power)<0.1)[0][0]



        ### check if peak_pos can be kept at middle of integration window.
        start_fine=True
        end_fine=True
               
        try:
            if smoothed[int(peak_pos-duration/2)]<=median_power+thres*mad:
                start_fine=False
        except IndexError:
            start_fine=False
        try:
            if smoothed[int(peak_pos+duration/2)]<=median_power+thres*mad:
                    end_fine=False
        except IndexError:
            end_fine=False

        if (start_fine==True and end_fine==True) or (start_fine==False and end_fine==False):
            start_time.append(data['axis_info']['time_axis']['MJDseconds'][int(max(0,peak_pos-duration/2))])
            end_time.append(data['axis_info']['time_axis']['MJDseconds'][int(min(peak_pos+duration/2,np.size(power)-1))])
        elif start_fine==False:
            diff=peak_pos-duration/2+1
            end_found=False
            while end_found==False:
                try:
                    if smoothed[int(diff)]>median_power+thres*mad:
                        end_found=True
                    else:
                        diff+=1
                except IndexError:
                    diff+=1
                

            distance=peak_pos-diff+1
            distance_left=duration-distance
            start_time.append(data['axis_info']['time_axis']['MJDseconds'][int(diff)])
            end_time.append(data['axis_info']['time_axis']['MJDseconds'][int(min(peak_pos+distance_left,np.size(power)-1))])
                
        else:
            diff=peak_pos+duration/2-1
            end_found=False
            while end_found==False:
                try:
                    if smoothed[int(diff)]>median_power+thres*mad:
                        end_found=True
                    else:
                        diff-=1
                except IndexError:
                    diff-=1
            distance=diff-peak_pos+1
            distance_left=duration-distance
            start_time.append(data['axis_info']['time_axis']['MJDseconds'][int(max(0,peak_pos-distance_left))])
            end_time.append(data['axis_info']['time_axis']['MJDseconds'][int(diff)])
        time_diff=start_time[-1]-startmjd
        startstr=(start+dt.timedelta(seconds=time_diff)).strftime('%Y/%m/%d/%H:%M:%S')
        start_datetime_obj.append(start+dt.timedelta(seconds=time_diff))
        time_diff=end_time[-1]-startmjd
        endstr=(start+dt.timedelta(seconds=time_diff)).strftime('%Y/%m/%d/%H:%M:%S')
        times.append(startstr+'~'+endstr)
        end_datetime_obj.append(start+dt.timedelta(seconds=time_diff))
        
    ms.close()
    flare_peak_time=start+dt.timedelta(seconds=peak_time-startmjd)
    peak_time_mjd=peak_time
    return times,flare,start_datetime_obj, end_datetime_obj,peak_time_mjd,startmjd,endmjd,flare_peak_time

def find_previous_image(imagedir, spw):
	imgprefix=imagedir+"selfcal"
	imagename = imgprefix + '_spw_*.image'
	images=glob.glob(imagename)
	num_img=len(images)
	if num_img==0:
		return False, ''
	s=np.zeros(num_img)
	for i in range(num_img):
		chunks=images[i].split('_')
		for n,c in enumerate(chunks):
			if c=='spw':
				break
		
		if '~' in chunks[n+1]:
			temp=chunks[n+1].split('~')
			freq1=int(temp[0])
			freq2=int(temp[1])
			s[i]=0.5*(freq1+freq2)
		else:
			s[i]=int(chunks[n+1])  ## the way I have named the images, the spw number will
						## always be just after the str "spw"
	diff=np.abs(-s+spw)
	ind=np.argsort(diff)
	non_zero=True
	for index in ind:
		if diff[index]!=0:
			return True, images[index]
	return False,''
	
			

def get_spw_num(msname):
    msmd.open(msname)
    nspws=msmd.nspw()
    msmd.done()
    return nspws

def calc_cellsize(msname):
    nspws=get_spw_num(msname)
    beams=[]
    ms.open(msname)
    for i in range(nspws):
        ms.selectinit(datadescid=0,reset=True)
        ms.selectinit(datadescid=i)
        data=ms.getdata(['u','v','axis_info'],ifraxis=True)
        fghz=data['axis_info']['freq_axis']['chan_freq'][0,0]/1e9   ### taking the first frequency. EOVSA spws bandwidths are small. Hence this is fine
        uvdist=np.sqrt(data['u']**2+data['v']**2)
        max_uvlambda=np.max(uvdist/299792458.0*fghz*1e9)
        beam_val=1.0/max_uvlambda*180/3.14159*3600/4.0  ### we take that the beam is resolved by 4 pixels by default
        beams.append(beam_val)
    ms.close()
    return beams

def get_ref_freqlist(msname):
    tb.open(msname + '/SPECTRAL_WINDOW')
    reffreqs = tb.getcol('REF_FREQUENCY')
    bdwds = tb.getcol('TOTAL_BANDWIDTH')
    cfreqs = reffreqs + bdwds / 2.
    tb.close()
    return cfreqs

def get_img_stat(imagename):
    ia.open(imagename+".image")
    data=ia.getchunk()[:,:,0,0]
    max_pix=np.nanmax(data)
    min_pix=np.nanmin(data)
    data[int(imsize*0.25):int(imsize*0.75),int(imsize*0.5):int(imsize*0.75)]=np.nan
    rms=np.nanstd(data)
    ia.close()
    return max_pix,min_pix,rms

def read_bandpass(CALTABLE,NSPW,NANT=16):	
    BCAL=tb.open(CALTABLE,nomodify=True) # Open the bandpass table 
    bpass=tb.getcol('CPARAM') # Extract the CPARAM=GAIN column from the bandpass table
    flag=tb.getcol('FLAG') # Extract FLAG column from bandpass table
    snr=tb.getcol('SNR');
    dims = bpass.shape	
    NPOL = dims[0]
    NCHAN = dims[1]
    NTIME = int(dims[2]/NANT)
    #if NCHAN!=NSPW:
     #   raise RuntimeError("Problem in caltable shape. Number of spws and channels in caltable does not match")
    flag = np.reshape(flag,[NPOL,NCHAN,NTIME,NANT])
    bpass=	np.reshape(bpass,[NPOL,NCHAN,NTIME,NANT])
    success=True
    tb.close()
    del BCAL
    del snr
    del NPOL
    del NCHAN
    del NTIME
    return bpass,flag, success

def combine_groups(group,pos):
    groups_com=[]
    deleted_groups=[]


    for g1 in range(len(group)):
        if g1 in deleted_groups:
            continue
        for g2 in range(g1+1,len(group)):
            combine=False
            for mem1 in group[g1]:
                for mem2 in group[g2]:
                    x1=pos[1][mem1]
                    x2=pos[1][mem2]
                    y1=pos[0][mem1]
                    y2=pos[0][mem2]
                    #if x1==x2 and y1==y2:
                    if np.sqrt((x1-x2)**2+(y1-y2)**2)<3:
                        combine=True
                        break
                if combine==True:
                    break
            if combine==True:
                for mem2 in group[g2]:
                    group[g1].append(mem2)
                deleted_groups.append(g2)
                continue
        temp=[]
        for mem in group[g1]:
            temp.append(mem)
        groups_com.append(temp)
    return groups_com

def gen_fof_groups(data3,thres):
    pos=np.where(data3>thres)
    y=pos[0]
    x=pos[1]
    size=len(x)

    groups=[[0]]
    for i in range(1,size):
        x0=pos[1][i]
        y0=pos[0][i]
        neighbour=False
        for g in groups:
            for mem in g:
                x1=pos[1][mem]
                y1=pos[0][mem]

                if np.sqrt((x1-x0)**2+(y1-y0)**2)<3:
                    g.append(i)
                    neighbour=True
                    break
            if neighbour==True:
                break
        if neighbour==False:
            groups.append([i])
    #print("Combining groups")
   # print(groups)
    group_com=combine_groups(groups,pos)
    #print("Groups combined")
    final_groups_x=[]
    final_groups_y=[]
    for g1 in group_com:
       # print (len(g1))
        if len(g1)>3:
            for mem in g1:
                final_groups_x.append(pos[1][mem])
                final_groups_y.append(pos[0][mem])
    return [final_groups_y,final_groups_x]

def gen_mask(image1,image2,mask1,mask2,threshold,imsize,filename,s,\
				make_shifted_mask=False, grow_threshold=0.5):  ### threshold corresponds to image1
	'''
	We will allow for small shifts and small change of size ere
	'''
	data2=fits.getdata(image2)[0,0,:,:]
	sidelobe_level,min_lev=find_sidelobe_level(image2)
	if image1!='':
		data1=fits.getdata(image1)[0,0,:,:]
		head1=fits.getheader(image1)
		data2=fits.getdata(image2)[0,0,:,:]
		pos1=np.where(data2==np.nanmax(data2))
		head2=fits.getheader(image2)
		#print (mask1)
		ia.open(mask1)
		mask1_data=ia.getchunk()[:,:,0,0].T
		ia.close()
		pos=np.where(mask1_data==1)
		#print (pos)
		if np.min(pos[0])==0 and np.min(pos[1])==0 and np.max(pos[1])==imsize-1 and np.max(pos[1])==imsize-1:
		    del pos
		    pos=np.where(data1>threshold)
		if np.size(pos)==0:
			return False
		x_min=np.min(pos[1])
		x_max=np.max(pos[1])
		y_min=np.min(pos[0])
		y_max=np.max(pos[0])
		box_x=x_max-x_min
		box_y=y_max-y_min
		#print(x_min,x_max,y_min,y_max)
		filename.write("spw="+str(s).zfill(2)+" Previous image supplied\n")
		filename.write("spw="+str(s).zfill(2)+" Mask boundaries in previous image:"+\
				str(x_min)+","+str(x_max)+","+str(y_min)+","+str(y_max))

		cell1=head1['CDELT1']
		cell2=head2['CDELT1']
		#print (cell1,cell2)
		ref_ra1=head1['CRVAL1']
		ref_dec1=head1['CRVAL2']
		ref_ra2=head2['CRVAL1']
		ref_dec2=head2['CRVAL2']
		ref_xpix1=head1['CRPIX1']
		ref_ypix1=head1['CRPIX2']
		ref_xpix2=head2['CRPIX1']
		ref_ypix2=head2['CRPIX2']

		ra_min=(x_min-ref_xpix1)*cell1+ref_ra1
		ra_max=(x_max-ref_xpix1)*cell1+ref_ra1
		dec_min=(y_min-ref_ypix1)*cell1+ref_dec1
		dec_max=(y_max-ref_ypix1)*cell1+ref_dec1

		x2_min=max(int((ra_min-ref_ra2)/cell2+ref_xpix2),0)
		x2_max=min(int((ra_max-ref_ra2)/cell2+ref_xpix2),imsize-1)
		y2_min=max(int((dec_min-ref_dec2)/cell2+ref_ypix2),0)
		y2_max=min(int((dec_max-ref_dec2)/cell2+ref_ypix2),imsize-1)
		#print (x2_min,x2_max,y2_min,y2_max)
		if x2_min>x2_max or y2_min>y2_max:  #### I do not know what to do here.
							#### Using the full image. Hopefull
			x2_min=0
			x2_max=imsize-1
			y2_min=0
			y2_max=imsize-1
		filename.write("spw="+str(s).zfill(2)+" Somehow lower boundary is greater than "+\
				"the outer boundary. Masing the full image\n")
		data3=data2[y2_min:y2_max,x2_min:x2_max]
		pos2=np.where(data3==np.nanmax(data3))
		if pos2[0][0]+y2_min!=pos1[0][0] or pos2[1][0]+x2_min!=pos1[1][0]:
			filename.write("spw="+str(s).zfill(2)+ "Peak pixel outside the mask."\
				" Check for position shifts\n")
			#filename.write("Trying to phase up at the location of the previous mask.\n")
			filename.write("spw="+str(s).zfill(2)+ "Difference in x:" +\
					str(int(pos1[1][0]-pos2[1][0]))+"pixels\n")
			filename.write("spw="+str(s).zfill(2)+"Difference in y:"+ \
					str(int((pos1[0][0]-pos2[0][0])))+"pixels\n")
		    ##### for testing #####
			if make_shifted_mask==False:
				filename.write("spw="+str(s).zfill(2)+ "Mske_shifted_mask is False. Will "+\
					"mask the entire region within the region of previous image\n")
				mask_y=np.arange(y2_min,y2_max+1)
				mask_x=np.arange(x2_min,x2_max+1)
				ia.open(mask2)
				data=ia.getchunk()
				data*=0
			    #for i,j in zip(mask_y,mask_x):
				for j in mask_x:
					for i in mask_y:
						data[j,i,0,0]=1                                           
				ia.putchunk(data)
				ia.close() 
				del data3                                 
				return False
		   
	else:
		y2_min=0
		y2_max=imsize
		x2_min=0
		x2_max=imsize
	data3=data2[y2_min:y2_max,x2_min:x2_max]
	#print (y2_min,y2_max,x2_min,x2_max)
	filename.write("spw="+str(s).zfill(2)+" Mask boundaries in previous image:"+\
				"{}\n".format(x2_min,x2_max,y2_min,y2_max))
	#print ("shift successfull")
	max_data3=np.nanmax(data3)
	chosen_thres=sidelobe_level+0.01 ### some more caution
	mask_y,mask_x=gen_fof_groups(data3,chosen_thres*max_data3)
	#print (mask_y, mask_x)
	if len(mask_x)==0:
		filename.write("spw="+str(s).zfill(2)+" Mask not found. Masking the region"+\
				" covered by the previous image\n")
		mask_y=np.arange(y2_min,y2_max+1)
		mask_x=np.arange(x2_min,x2_max+1)
		ia.open(mask2)
		data=ia.getchunk()
		data*=0
		#for i,j in zip(mask_y,mask_x):
		for j in mask_x:
			for i in mask_y:
			    data[j,i,0,0]=1                                           
		ia.putchunk(data)
		ia.close()      
		return False
	for i in range(len(mask_y)):
		mask_y[i]+=y2_min
		mask_x[i]+=x2_min
	ia.open(mask2)
	data=ia.getchunk()
	data*=0
	for i,j in zip(mask_y,mask_x):
		data[j,i,0,0]=1
	#print (i,j)

	ia.putchunk(data)
	ia.close()
	filename.write("spw="+str(s).zfill(2)+" Mask generated. Growing the mask"+\
			" with threshold {}\n".format(grow_threshold))
	success=grow_mask(image2,mask2,max(grow_threshold,min_lev-0.1)*max_data3)  ### again some safety here.
	filename.write("spw="+str(s).zfill(2)+" Mask generated\n")
	return True

def confirm_maximum_pixel(imagename,mask,freq_list,spwran,msname, uvrange, \
						imsize, cell, phasecenter,filename,s):
	try:
		first_freq=int(spwran.split('~')[0])
		last_freq=int(spwran.split('~')[1])
	except IndexError:
		first_freq=int(spwran)
		last_freq=first_freq
	#print (first_freq, last_freq)

	filename.write("spw="+str(s).zfill(2)+" confirm_maximum_pixel called with {}\n".format(spwran))
	maxpix=imstat(imagename,mask="mask(\'"+mask+"\')")["maxpos"]
	
	first_freq=max(minimum_spw,first_freq-1)
	last_freq=min(last_freq+1,maximum_spw)

	#if abs(freq_list[last_freq]-freq_list[first_freq])>0.75*freq_list[0.5*(first_freq+last_freq)]:
	#	print ("Suitable frequency averaging not possible")
	#	return True
	spwran=str(first_freq)+"~"+str(last_freq)
	filename.write("spw="+str(s).zfill(2)+" Checking against "+str(first_freq)+"~"+str(last_freq)+"\n")
	
	tclean(vis=msname, imagename="check_maxpix_robust",uvrange=uvrange,\
			 spw=spwran, imsize=imsize,cell=cell, niter=1, gain=0.05,\
			 cyclefactor=10, weighting='briggs', robust=0.0,savemodel='none',\
			 pbcor=False,pblimit=0.001, phasecenter=phasecenter,stokes='XX')

	maxpix1=imstat("check_maxpix_robust.image",mask="mask(\'"+mask+"\')")["maxpos"]
	
	shift=np.array([maxpix[0]-maxpix1[0],maxpix[1]-maxpix1[1]])
	filename.write("spw="+str(s).zfill(2)+" Shift is {}\n".format(shift))
	cellsize=float(cell.split('arcsec')[0])
	success=check_shift("check_maxpix_robust.image",shift,cellsize)
	filename.write("spw="+str(s).zfill(2)+" Shift ok? {}\n".format(success))
	os.system("rm -rf check_maxpix_robust.*")
	return success

def restore_previous_condition(imagename):
	temp=imagename.split('_')
	prefix='_'.join(temp[:-1])
	iteration_num=int(temp[-1])
	psf=glob.glob(prefix+"*.psf")
	num_psf=len(psf)
	for i in range(iteration_num+1,num_psf):
		image=prefix+'_'+str(i).zfill(2)
		for str1 in ['.image','.mask','.residual','.psf',\
					'.model','.sumwt','.pb','.fits']:
			
			os.system("rm -rf "+image+str1)
		os.system("rm -rf "+image+".gcal")
	return
		


def do_selfcal(msname, sp,spwran, imagedir,caltbdir, imsize, cell, maxiter,refantenna,nspw, uvrange, \
		filename, phasecenter, freq_list, ref_image='',make_shifted_mask=False, combine_spws=False):
	clearcal(msname)
	calprefix = caltbdir + 'selfcal'
	imsize=2048
	imgprefix = imagedir + 'selfcal'
	os.system('rm -rf '+caltbdir+'*_spw_'+sp.zfill(2)+'*')
	os.system('rm -rf '+imagedir+'*spw_'+sp.zfill(2)+'*')
	iteration_num=0
	found_mask=False
	change_grow_mask_threshold=False
	grow_mask_threshold=0.7
	image_to_go_back_to=''
	filename.write("spw="+str(sp).zfill(2)+" grow_mask_threshold={}\n".format(grow_mask_threshold))
	while iteration_num<maxiter:	
		#print("starting round "+str(iteration_num))
		caltable = calprefix + '_spw_'+sp.zfill(2)+'_'+str(iteration_num).zfill(2)+'.gcal'
		imagename = imgprefix + '_spw_' + sp.zfill(2) + '_' +str(iteration_num).zfill(2)
		
		filename.write("spw="+str(sp).zfill(2)+" Creating the first diry map\n") 
		if iteration_num==0:		
			tclean(vis=msname, imagename=imagename,uvrange=uvrange, spw=spwran, imsize=imsize,\
			cell=cell, niter=1, gain=0.05, cyclefactor=10, weighting='briggs',stokes='XX',\
			robust=0.0,savemodel='none', pbcor=False,pblimit=0.001, phasecenter=phasecenter)
			
			#print ("Check if image produced or not")
			if os.path.isdir(imagename+".image")==False:
				#print ("image not produced")
				filename.write("spw="+str(sp).zfill(2)+" Image not produced\n")
				os.system('rm -rf '+caltbdir+'*_spw_'+sp.zfill(2)+'*')
				os.system('rm -rf '+imagedir+'*spw_'+sp.zfill(2)+'*') 
				return False, False
			#print ("Getting image statistics")
			filename.write("spw="+str(sp).zfill(2)+" Getting image statistics\n")
			max_pix,min_pix,rms=get_img_stat(imagename)
			#print (max_pix)
			threshold=10*rms
			if max_pix<15*rms and len(ref_image)==0:  
				filename.write("spw="+str(sp).zfill(2)+" Image SNR too low for selfcal with no mask\n")
				filename.write("spw="+str(sp).zfill(2)+" Trying to find past image\n")
				#print("Trying to find past image")
				past_image_found, past_image=find_previous_image(imagedir, int(sp))
				if past_image_found==True:
					os.system('rm -rf '+caltbdir+'*_spw_'+sp.zfill(2)+'*')
					os.system('rm -rf '+imagedir+'*spw_'+sp.zfill(2)+'*')   
					filename.write("spw="+str(s).zfill(2)+" Previous mask files found."+\
						"Trying to use the previous mask\n")
					return False, True
				filename.write("spw="+str(sp).zfill(2)+" Previous mask does not exist. Autogenerating mask\n")
				threshold=0.5*max_pix
				exportfits(imagename=imagename+".image",fitsimage=imagename+".fits")
				found_mask=gen_mask('',imagename+".fits",'',imagename+".mask",10*rms,\
						imsize, filename,sp,make_shifted_mask,grow_mask_threshold)
				if found_mask==False:
					filename.write("spw="+str(sp).zfill(2)+" Mask not produced. Exiting\n")
					#print("Mask not produced. Exitting.")
					os.system('rm -rf '+caltbdir+'*_spw_'+sp.zfill(2)+'*')
					os.system('rm -rf '+imagedir+'*spw_'+sp.zfill(2)+'*')
					return False, True
		    
		    #elif max_pix<10*rms and len(ref_image)==0:
		     #   print ("Source not detected")
		      #  return False, False
		else:
		    threshold=10*rms
		    if combine_spws==False:
		    	spwran=sp
		    
		#print ("Starting actual clean")
		mask=''
		if len(ref_image)!=0 and iteration_num==0:
			#print (len(ref_image),ref_image)
			filename.write("spw="+str(sp).zfill(2)+" Generating mask using the previous mask\n")
			#print ("Generating mask using the previous mask\n")
			threshold=0.5*max_pix
			exportfits(imagename=imagename+".image", fitsimage=imagename+".fits")
			max_pix,min_pix,rms=get_img_stat(ref_image[:-5])
			found_mask=gen_mask(ref_image,imagename+".fits",ref_image[:-5]+".mask",\
				imagename+".mask",10*rms,imsize,filename,sp, make_shifted_mask,\
				grow_mask_threshold)

			if found_mask==False: 
				filename.write("spw="+str(sp).zfill(2)+" Mask not generated."+\
					" Trying to do more averaging in frequency\n ")
				#print("Mask  not generated. Trying to do more averaging in frequency")
				return False, True

			#print(imagename+".image")
			filename.write("spw="+str(sp).zfill(2)+" Checking if the maximum pixel inside the "+\
					"mask is in the correct position\n")
			success=confirm_maximum_pixel(imagename+".image",imagename+".mask",freq_list,spwran\
										,msname, uvrange, imsize, cell, phasecenter,filename,sp)

			if success==False:
				filename.write("spw="+str(sp).zfill(2)+" Probable shift in image plane."+\
					" Do more averaging in frequency.\n ")
				print("Probable shift in image plane. Do more averaging in frequency.")
				return False, True
			
		    

		elif len(ref_image)!=0:
		    #max_pix,min_pix,rms=get_img_stat(imagename+'.image')
		    threshold=0.2*max_pix
		    max_pix,min_pix,rms=get_img_stat(ref_image[:-5])
		    exportfits(imagename= imagename[:-2]+str(iteration_num-1).zfill(2)+".image", fitsimage=\
		            imagename[:-2]+str(iteration_num-1).zfill(2)+".fits",overwrite=True)
		    #print("Fits file generated")
		    #print (imagename, iteration_num)
		    found_mask=gen_mask(ref_image, imagename[:-2]+str(iteration_num-1).zfill(2)+".fits",\
					ref_image[:-5]+'.mask',imagename[:-2]+"00.mask",10*rms,imsize,filename,\
					sp,make_shifted_mask,grow_mask_threshold)  ### if found_mask=False, \
													#### the earlier mask was \
		                                            #### not updated and hence can be used.
		    #print (found_mask)
		    
		    mask=imagename[:-2]+"00.mask"
		    #print (max_pix)

		elif iteration_num!=0 and len(ref_image)==0 and max_pix<15*rms:
			threshold=0.2*max_pix
			exportfits(imagename=imagename[:-2]+str(iteration_num-1).zfill(2)+".image",fitsimage=\
		            imagename[:-2]+str(iteration_num-1).zfill(2)+".fits",overwrite=True)
			found_mask=gen_mask('',imagename[:-2]+str(iteration_num-1).zfill(2)+".fits",\
					'',imagename[:-2]+"00.mask",10*rms,imsize,filename,sp,make_shifted_mask,\
					grow_mask_threshold)
			mask=imagename[:-2]+"00.mask"
		#print (str(threshold)+"Jy")
		filename.write("spw="+str(sp).zfill(2)+" Cleaning threshold:{}\n".format(threshold))	
		#if no_shift_mask==True and iteration_num>3:
		 #   raise RuntimeError("shift not in image")
		if found_mask==True or max_pix>15*rms:
			tclean(vis=msname, imagename=imagename, spw=spwran, uvrange=uvrange,\
				 imsize=imsize,threshold=str(threshold)+"Jy",cell=cell, niter=10000,\
				 gain=0.05, cyclefactor=10, weighting='briggs',robust=0.0,\
				savemodel='modelcolumn', pbcor=False,pblimit=0.01,mask=mask,\
				phasecenter=phasecenter,stokes='XX')
			
			### sometimes it happens that due to some error, a blank model is generated.
			### But this should not be the case, as at least the maximum pixel should be
			### picked. I saw this happening in the case of 20210507 spw 45. CASA said 
			### and I am copying it " Caught Exception: NonLinearFitLM: error in loop
			### solution" and then it restored a blank image. In these situations, gaincal
			### will essentially use a point source at the phase center as a model. But, this
			### should not be the case. Hence, I will check if the model is blank or not. 
			### If yes, I will assume that some problem has happened and I will not proceed
					
			ia.open(imagename+".model")
			modeldata=ia.getchunk()
			ia.close()
			if np.nansum(modeldata)<1e-6:
				#print("Blank model. Checking if this is because of very high threshold.")
				filename.write("spw="+str(sp).zfill(2)+" Blank model. Checking if this"+\
					" is because of very high threshold.\n")
				max_pix=imstat(imagename=imagename+".image",mask="mask(\'"+mask+"\')")["max"]
				if max_pix<threshold:
					threshold=max_pix*0.5
					filename.write("spw="+str(sp).zfill(2)+" Threshold was too high."+\
							"new threshold={}\n".format(threshold))
					tclean(vis=msname, imagename=imagename, spw=spwran, uvrange=uvrange,\
				 		imsize=imsize,threshold=str(threshold)+"Jy",cell=cell, niter=10000,\
				 		gain=0.05, cyclefactor=10, weighting='briggs',robust=0.0,\
						savemodel='modelcolumn', pbcor=False,pblimit=0.01,mask=mask,\
						phasecenter=phasecenter,stokes='XX')
					ia.open(imagename+".model")
					modeldata=ia.getchunk()
					ia.close()
					if np.nansum(modeldata)<1e-6:
						filename.write("spw="+str(sp).zfill(2)+" Model still blank. Exiting\n")
						return False, True
				else:
					return False, True  ### I am returning True here, as I am not really sure
									### if it cannot recover ever. It is better to let it
									### try, rather than leaving it False. However, my gut
									### feeling is that, it would not be able to recover.

		else:
			#print("Mask not found or rms is very high")
			return False, True


		max_pix,min_pix,rms=get_img_stat(imagename)
		dyn_range=max_pix/rms
		print ("dyn_range=",dyn_range, iteration_num)
		filename.write("spw="+str(sp).zfill(2)+" dynamic range:"+str(dyn_range)+", Iteration"+\
				" number:"+str(iteration_num)+"\n")
		if iteration_num<=1:
			dyn_range1=dyn_range
			min_pix1=min_pix
			max_pix1=max_pix
		else:
			if dyn_range<10 and dyn_range<dyn_range1*0.95:##len(ref_image)==0 and mask=='':
				filename.write("spw="+str(sp).zfill(2)+" SNR too low for selfcal and SNR decreasing\n")
				print("SNR too low for selfcal")
				if change_grow_mask_threshold==True:
					print ("spw="+str(sp).zfill(2)+" Reverting back and exiting with success")
					restore_previous_condition(image_to_go_back_to)
					return True, True
				return False, True
			elif dyn_range/dyn_range1>0.95 and dyn_range/dyn_range1<1.05 and iteration_num>2: ### otherwise, this means no more improvement
				print (dyn_range, max_pix/abs(min_pix))
				
				if change_grow_mask_threshold==False:
					change_grow_mask_threshold=True
					grow_mask_threshold=0.2
					filename.write("spw="+str(sp).zfill(2)+" Changing grow_mask_threshold to"+\
							"{}\n".format(grow_mask_threshold))
					image_to_go_back_to=imagename
					print (imagename)
				else:
					filename.write("spw="+str(sp).zfill(2)+" Converged\n")
					return True, True
			elif dyn_range/dyn_range1<1.1 and iteration_num==5 and spwran==sp:
				filename.write("spw="+str(sp).zfill(2)+" Not enough improvement.\n")
				print("Not enough improvement")
				if change_grow_mask_threshold==False:
					change_grow_mask_threshold=True
					grow_mask_threshold=0.5
					image_to_go_back_to=imagename
					filename.write("spw="+str(sp).zfill(2)+" Changing grow_mask_threshold to"+\
							"{}\n".format(grow_mask_threshold))
				else:
					return True, True

		#if iteration_num>=1:
	#	flagmanager(vis=slfcalms,mode='restore',versionname='applycal_1')
	#	flagmanager(vis=slfcalms,mode='delete',versionname='applycal_1')
		if os.path.isdir(caltable):
		    os.system("rm -rf "+caltable)
		gaincal(vis=msname, refant=refantenna, spw=sp,caltable=caltable,uvrange=uvrange, \
			calmode='p', combine='scan', minblperant=4, minsnr=3, append=False,\
			 solnorm=True,solmode='L1R',rmsthresh=[10,7],normtype='median')

		if os.path.isdir(caltable)==False:
			filename.write("spw="+str(sp).zfill(2)+" Caltable not produced\n")
			print ("Solution not found")
			if change_grow_mask_threshold==True:
				restore_previous_condition(image_to_go_back_to)
				print ("spw="+str(sp).zfill(2)+" Reverting back and exiting with success")
				return True, True
			return False,True
		
		bpass,flag,success=read_bandpass(caltable,nspw)
		pos=np.where(flag[0,0,0,:]==True)[0]
		num_flagged_ant=np.size(pos)
		if iteration_num==0:
		    num_flagged_ant1=num_flagged_ant
		filename.write("Antennas flagged in caltable:" +str(pos)+'\n')
		if num_flagged_ant>num_flagged_ant1+3:
			filename.write("spw="+str(sp).zfill(2)+" Flagged antennas have increased with iteration. "+\
						"Possible divergence. Exiting\n")
			print("Flagged anttenae number increased a lot")
			if change_grow_mask_threshold==True:
				print ("spw="+str(sp).zfill(2)+" Reverting back and exiting with success")
				restore_previous_condition(image_to_go_back_to)
				return True, True
			return False, True

		
		if 15-len(pos)<5: ### will not proceed with this selfcal if 5 or more antennas are flaggged
			filename.write("spw="+str(sp).zfill(2)+" Too many flagged antennas. Exiting\n")
			print("Too many flagged antennas")
			if change_grow_mask_threshold==True:
				print ("spw="+str(sp).zfill(2)+" Reverting back and exiting with success")
				restore_previous_condition(image_to_go_back_to)
				return True, True
			return False, True
		
		### TODO Need to check if this helps the process or makes it worse. Initial guess is it makes it worse.
		applycal(vis=msname,spw=sp,gaintable=caltable,applymode='calonly')
		if iteration_num!=maxiter-1:
		    delmod(msname,scr=True)
		iteration_num+=1
		dyn_range1=dyn_range
		min_pix1=min_pix
		max_pix1=max_pix
		num_flagged_ant1=num_flagged_ant
	filename.write("spw="+str(sp).zfill(2)+" Successfull exit\n")
	return True, True 
        
def gen_blank_cal(msname, spw, caltbdir):
    calprefix=caltbdir+'selfcal'
    iteration_num=0 
    caltable = calprefix + '_spw_'+spw.zfill(2)+'_'+str(iteration_num).zfill(2)+'.gcal'
    os.system("rm -rf "+calprefix+'*')
    gencal(vis=msname, caltable=caltable,caltype='amp', spw=spw,antenna='', pol='', parameter=[1.0])
    return
    
    
def calling_do_selfcal(slfcalms, s, imagedir,caltbdir, \
			casa_imsize, cell_val,maxiter,refantenna,nspw,\
			uvrange,filename,phasecenter,freqs,avg_spw_max):
						
	sp=str(s)
	ref_image=''
	print('processing spw: ' + sp)

	filename.write("spw="+str(s).zfill(2)+" Cell size="+str(cell_val)+"\n")
	#f.write("Imsize="+str(imsize)+'\n')
	filename.write("spw="+str(s).zfill(2)+" Uvrange= "+uvrange+'\n')


	success,signal=do_selfcal(slfcalms, sp,sp, imagedir,caltbdir, \
					casa_imsize, cell_val,maxiter,refantenna,nspw,\
					uvrange=uvrange,filename=filename,phasecenter=phasecenter,freq_list=freqs)

	if success==False and signal == False:
		return success  ### I would only put signal to be False if there is
							### something seriously wrong in the data

				
	if success==False:
		filename.write("spw="+str(s).zfill(2)+" Trying with previous mask\n")
		print ("Trying with previous mask \n \n \n \n")
		### here I will try to find the image at nearest spw. Assumption is 
		### that if it is at the nearest frequency, it is very likely that
		### the sources at both of them will be located close by.
		past_image_found, past_image=find_previous_image(imagedir, s)
		if past_image_found==True:
			filename.write("spw="+str(s).zfill(2)+" Past image found. Using already used masks\n")
			fitsimage=past_image[:-6]+".fits"
			if os.path.isfile(fitsimage)==False:
				exportfits(imagename=past_image, fitsimage=fitsimage)

			ref_image=fitsimage
			success,signal=do_selfcal(slfcalms, sp,sp, imagedir,caltbdir,\
					casa_imsize, cell_val,maxiter,refantenna,nspw,uvrange=uvrange,\
					ref_image=ref_image,filename=filename,phasecenter=phasecenter,freq_list=freqs)

	if success==False:  ### this means that we would need to do 
						### multi-frequency synthesis. To do that
						### first I would make sure that the maximum
						### pixel within the previous mask in the 
						### dirty map is a real feature and not an
						### artefact.
		print ("Trying to do mfs \n \n \n \n")
		filename.write("spw="+str(s).zfill(2)+" Trying to do mfs\n")
		imgprefix=imagedir+'selfcal'
		imagename = imgprefix + '_spw_' + sp.zfill(2) + '_00'
		mask=imagename+".mask"
		maxpix=imstat(imagename+".image",mask="mask(\'"+mask+"\')")["maxpos"]

	###TODO A big assumption is that the source is at the same location at all 
	### frequencies, in the sense that they lie in the box created by the lower 
	### frequency. If the box is very small and the true source location moves, 
	###  then there is a problem. Hence, we need to do another check. We shall 
	### take a box of size lets say 5 arcminutes centred on the original mask and
	### also caluclate the maximum pixel there. If the maximum pixel inside this 
	### box remains same even when we do mfs, then probably the true source is 
	### also shifted. 
		freq_int_found=False		
		for avg_spw in range(1,avg_spw_max):
			min_spw=max(minimum_spw,s-avg_spw)
			max_spw=min(maximum_spw,s+avg_spw)
			#### implementing a check where I ensure that I will never average more 
			#### more than 0.5 times the central frequency
			if abs(freqs[max_spw]-freqs[min_spw])>max_frac_freq_avg*freqs[s]:
				filename.write("Did not a suitable averaging range in frequency."+\
						"Trying with maximum possible frequency bandwidth.\n")
				break
			if s-avg_spw==min_spw-1 and s+avg_spw==max_spw:
				f.write("Did not a suitable averaging range in frequency. Exiting.")
				break
			if min_spw!=max_spw:
				spwran=str(min_spw)+'~'+str(max_spw)
			else:
				spwran=sp
			print('processling spw {0:s} using spw {1:s} as model'.format(sp, spwran))
			filename.write('Processing spw {0:s} using spw {1:s} as model\n'.format(sp,spwran))
	    		
			if os.path.isdir("check_maxpix.psf"):
	    			os.system("rm -rf check_maxpix*")
			tclean(vis=slfcalms, imagename='check_maxpix', spw=spwran, uvrange=uvrange, \
						imsize=2048,cell=cell_val, niter=0, gain=0.05, cyclefactor=10,\
						 weighting='briggs',stokes='XX',\
	 			robust=0.0,savemodel='none', pbcor=False,pblimit=0.01,\
						mask='user',phasecenter=phasecenter)
			maxpix1=imstat("check_maxpix.image",mask="mask(\'"+mask+"\')")["maxpos"]
			print (maxpix, maxpix1)
			if maxpix[0]!=maxpix1[0] and maxpix[1]!=maxpix1[1]:
				shift_ok=check_shift("check_maxpix.image",\
						np.array([maxpix[0]-maxpix1[0],maxpix[1]-maxpix1[1]]),\
						cell[s])  ## 2 is the cellsize
			else:
				shift_ok=True
			print ("not an artefact", shift_ok)
			filename.write("spw="+str(s).zfill(2)+" Not an artefact? {}\n".format(shift_ok))
			if shift_ok==True:
				freq_int_found=True													
				break
			else:
				maxpix=maxpix1
				os.system("rm -rf check_maxpix.*")
		if freq_int_found==False:
			filename.write("spw="+str(s).zfill(2)+" Doing the best we can get and trying it out\n")
			s-=1
			min_spw=max(minimum_spw,s-avg_spw)
			max_spw=min(maximum_spw,s+avg_spw)
			freq_int_found=True
		print (min_spw, max_spw)
		f.write("spw="+str(s).zfill(2)+":"+str(min_spw)+","+str(max_spw))
		os.system("rm -rf check_maxpix.*")		
		if freq_int_found==True:
			for avg_spw1 in range(avg_spw,avg_spw_max): 
				min_spw=max(minimum_spw,s-avg_spw1)
				max_spw=min(maximum_spw,s+avg_spw1)
				if s-avg_spw==min_spw-1 and s+avg_spw==max_spw:
					success=False
					break
				if min_spw!=max_spw:
					spwran=str(min_spw)+'~'+str(max_spw)
				else:
					spwran=sp
					
				print ("Calling do_selfcal with mfs")
				filename.write("spw="+str(s).zfill(2)+":"+str(min_spw)+","+str(max_spw)+"\n")
				success,signal=do_selfcal(slfcalms, sp,spwran, imagedir,caltbdir,\
						 	casa_imsize,cell_val,maxiter,refantenna,\
							nspw,uvrange=uvrange,ref_image=ref_image,\
							filename=filename, make_shifted_mask=True,phasecenter=phasecenter\
							,freq_list=freqs)    
				print("success=",success)
				f.write("spw="+str(s).zfill(2)+" success={}".format(success))
		            
				if success==True and signal==True:
					break
		else:
			success=False
		
	if success==False and freq_int_found==True:
		print ("Calling do_selfcal with mfs and solving gains at all frequencies")
		filename.write("spw="+str(s).zfill(2)+":"+str(min_spw)+","+str(max_spw)+"\n")
		success,signal=do_selfcal(slfcalms, sp,spwran, imagedir,caltbdir,\
				 	casa_imsize,cell_val,maxiter,refantenna,\
					nspw,uvrange=uvrange,ref_image=ref_image,\
					filename=filename, make_shifted_mask=True,phasecenter=phasecenter\
					,freq_list=freqs,combine_spws=True)    
		print("success=",success)
		f.write("spw="+str(s).zfill(2)+" success={}".format(success))
		
	return success

def get_img_center_heliocoords(images):
    num_img=len(images)
    x=np.zeros(num_img)
    y=np.zeros(num_img)
    for j,img in enumerate(images):
        head=fits.getheader(img)
        data=fits.getdata(img)[0,0,:,:]
        x0=head['CRVAL1']
        y0=head['CRVAL2']
        cell_x=head['CDELT1']
        cell_y=head['CDELT2']
        pos=np.where(data==np.nanmax(data))
        xmax=pos[1][0]
        ymax=pos[0][0]
        xcen=head['CRPIX1']
        ycen=head['CRPIX2']
        dx=xmax-xcen
        dy=ymax-ycen
        dx_asec=dx*cell_x+x0
        dy_asec=dy*cell_y+y0
        x[j]=dx_asec
        y[j]=dy_asec
    xycen=[np.median(x),np.median(y)]
    return xycen

	

os.chdir(workpath)
workdir=os.getcwd()+"/"
try:
	os.mkdir(slfcaldir)
except OSError:
	pass
	
try:
	os.mkdir(imagedir)
except OSError:
	pass
	
try:
	os.mkdir(caltbdir)
except OSError:
	pass


#outpath = 'msdata/'
#os.chdir(outpath)
#msfiles=glob.glob("*.ms")
#for file1 in msfiles:
 #   if '_' in file1:
  #     concatvis=file1

#print (concatvis)

rawcal_ms='IDB20220930_1600-1700.ms'
os.chdir(workpath)
if os.path.isdir(rawcal_ms)==False:
	print ("Splitting here")
	split(vis=outpath + concatvis, outputvis=rawcal_ms, correlation='XX,YY', datacolumn='data',
     antenna='')

# ============ Split a short time for self-calibration ===========
print ("Starting the self-calibration process")
f=open("selfcal.log","a")
ms_in = rawcal_ms
if cal_disk:
	ms_caldisk = concatvis[:-3]+'_fullday_calibrated.ms'
	slftbs = glob.glob('fullday_cals/*')
	if os.path.exists(ms_caldisk)==False:
		applycal(vis=ms_in, gaintable=slftbs, selectdata=True, \
	     		interp='linear', flagbackup=False, applymode='calonly', calwt=False)
		#os.system('rm -rf ' + ms_caldisk)
		split(vis=ms_in, outputvis=ms_caldisk, datacolumn='corrected', width=1)
	ms_in = ms_caldisk


#============ Starting the selfcal loop ================

if doslfcal:
	print ("Determining the selfcal timeranges")
	times,flares,start_time_obj, end_time_obj,peak_time_mjd,startmjd,endmjd,flare_peak_time=get_selfcal_times(ms_in,maximum_spw+1,starttime)

 # =========== Step 3, main step of selfcalibration =========

	

	if len(flag_antennas)!=0:
		 flagdata(vis=ms_in,mode='manual',antenna=flag_antennas)

	split(vis=ms_in,outputvis='temp_ms.ms',spw=selfcal_spw,datacolumn='data',timerange=times[0])
	slfcalms='temp_ms.ms'
	nspw=get_spw_num(slfcalms)
		
	freqs=get_ref_freqlist(slfcalms)
	cell_val=[]
	uvranges=[]
	if calc_cell==True:
		cell=calc_cellsize(slfcalms)
		for i in range(len(cell)):
		    cell_val.append(str(cell[i])+'arcsec')
		    uvranges.append('>'+str(uvlim*freqs[i]/freqs[0])+'lambda') 
	elif len(beam)<nspw:
		print("Number of beams provided does not match number of spw."+
			" Using the first beam value and then I will scale with frequency") 
		
		for i in range(len(freqs)):
		    cell_val.append(str(cell[0]*freqs[0]/freqs[i])+'arcsec')
		    uvranges.append('>'+str(uvlim*freqs[i]/freqs[0])+'lambda')
		### scaled the beam here with frequency
	else:
		for i in range(nspw):
		    cell_val.append(str(cell[i])+'arcsec')
		    uvranges.append('>'+str(uvlim*freqs[i]/freqs[0])+'lambda') 

	if calc_imsize==True:
		imsize=42*16*60  ### 42 solar radii 
	casa_imsize=imsize/cell[0]  ## the imsize stays same as we are scaling both FOV and cell size with frequency
	os.system("rm -rf temp_ms.ms")

if doslfcal:
	if phasecenter=='':
		### then I need to find the phase center. To do that is to first make a big dirty 
		### map starting from the lowest frequency. Then I shall check if the DR of
		### the dirty map is at least 8. If not, I shall proceed to the next frequency. 
		### Once, I find a source detection, I shall set the peak location to the phasecenter.
		f.write("Finding phasecenter\n")
		print("Finding phasecenter")
		imagename="find_phasecenter"
		find_phase_center=False
		ra=[]
		dec=[]
		j=0
		print (flares)
		for s in range(5,45):  ### I will find phasecenter only at these frequencies
			#if flares[s-minimum_spw]==False:
			#	continue
			if j==6:
				break
			
			### TODO This following line should be times[s]. I changed that for testing only.
			
			current_trange=times[s-minimum_spw]
			ms_slfcaled = workdir + os.path.basename(ms_in).replace('.ms', '.slfcaled_v2.ms') 
							# output, selfcaled, visibility
			slfcalms = slfcaldir + 'slfcalms.XX.slfcal'  
				# intermediate small visibility for selfcalbration 
			print ("Spliting from "+ms_in)
			if os.path.exists(slfcalms):
				os.system("rm -rf "+slfcalms+"*")
			split(vis=ms_in, outputvis=slfcalms, datacolumn='data', \
					timerange=current_trange,correlation='XX', spw=selfcal_spw)

			tclean(vis=slfcalms, imagename=imagename, spw=str(s), uvrange=uvranges[s], \
					imsize=4096,cell=cell_val[s], niter=0, gain=0.05, cyclefactor=10,\
					weighting='briggs',robust=0.0,savemodel='none', pbcor=False,
					pblimit=0.01,stokes='XX')
			max_pix,min_pix,rms=get_img_stat(imagename)
			if max_pix/rms>10:
				j+=1
				find_phase_center=True
				pos=imstat(imagename+".image")['maxposf']
				temp=pos.split(',')
				ra_str=temp[0].split(':')
				dec_str=temp[1].split('.')
				ra.append((abs(int(ra_str[0]))+int(ra_str[1])/60.0+float(ra_str[2])/3600.)*15)
				sign=1
				if '-' in ra_str[0]:
					sign=-1
				ra[-1]=ra[-1]*sign
				sign=1
				if '-' in dec_str[0]:
					sign=-1
				try:
					dec.append(abs(int(dec_str[0]))+int(dec_str[1])/60.0+float(dec_str[2]+'0.'+dec_str[3])/3600.)
				except IndexError:
					dec.append(abs(int(dec_str[0]))+int(dec_str[1])/60.0+float(dec_str[2])/3600.)
				dec[-1]=dec[-1]*sign
				os.system("rm -rf "+imagename+".*")
				os.system("rm -rf "+slfcalms+"*")
			else:
				os.system("rm -rf "+imagename+".*")
				os.system("rm -rf "+slfcalms+"*")

		if find_phase_center==True:
			ra_final=np.median(np.array(ra))
			dec_final=np.median(np.array(dec))
			phasecenter='J2000 '+str(ra_final)+"deg "+str(dec_final)+"deg"
			print (phasecenter)
			f.write("Phasecenter:{}\n".format(phasecenter))
			os.system("rm -rf "+imagename+".*")
		else:
			doslfcal=False
			f.write("Appropriate phase center could not be found."+\
					"Please restart after providing an appropriate one.\n")
			print("Appropriate phase center could not be found."+\
					"Please restart after providing an appropriate one.")

if doslfcal:
	for s in range(50):
		if s not in selfcal_spws:
			continue
		success=False
		start_time=start_time_obj[s-minimum_spw]
		end_time=end_time_obj[s-minimum_spw]
		time_delta=(end_time-start_time).seconds
		if flares[s-minimum_spw]==True:
			max_time_delta=180 ## 3 minutes
			if time_delta<20:
				time_increase=5
			elif time_delta>20 and time_delta<50:
				time_increase=10
			else:
				time_increase=15
		else:
			max_time_delta=300 ## 5 minutes
			time_increase=30
		while success!=True:
			time_delta=(end_time-start_time).seconds
			if time_delta>max_time_delta:
				break
			startstr=start_time.strftime('%Y/%m/%d/%H:%M:%S')
			endstr=end_time.strftime('%Y/%m/%d/%H:%M:%S')
			current_trange=startstr+"~"+endstr
			f.write("spw="+str(s).zfill(2)+" Self-calibrating spw "+str(s)+\
						" with time range "+current_trange+'\n')

			ms_slfcaled = workdir + os.path.basename(ms_in).replace('.ms', '.slfcaled_v2.ms') # output, selfcaled, visibility
			slfcalms = slfcaldir + 'slfcalms.XX.slfcal'  # intermediate small visibility for selfcalbration 
			slfcaledms = slfcaldir + 'slfcalms.XX.slfcaled'
			print ("Spliting from "+ms_in)
			f.write("spw="+str(s).zfill(2)+" Splitting..\n")
			if not os.path.exists(slfcalms):
				split(vis=ms_in, outputvis=slfcalms, datacolumn='data', timerange=current_trange, correlation='XX', spw=selfcal_spw)
			if identify_data_gap:
				flag_data_gap(slfcalms,s)
		
			print('Processing ' + current_trange)
			
			success=calling_do_selfcal(slfcalms, s, imagedir,caltbdir, \
				casa_imsize, cell_val[s],maxiter,refantenna,nspw,\
				uvranges[s],f,phasecenter,freqs,avg_spw_max)
			f.write("success="+str(success))
			if success!=True:
				start_time=start_time-dt.timedelta(seconds=time_increase)
				end_time=end_time+dt.timedelta(seconds=time_increase)
			os.system("rm -rf "+slfcalms+"*")
			break
		sp=str(s)
		f.write("spw="+str(s).zfill(2)+" success={}\n".format(success))        
		calprefix=caltbdir+'selfcal'
		caltable = calprefix + '_spw_'+sp.zfill(2)
		imgprefix=imagedir+"selfcal"
		image = imgprefix + '_spw_' + sp.zfill(2)
		
		if success==False:
			for str1 in ['.image','.mask','.residual','.psf',\
					'.model','.sumwt','.pb','.fits']:
				os.system("rm -rf "+image+"*"+str1)
			os.system("rm -rf "+caltable+"*.gcal")
			gen_blank_cal(ms_in,sp, caltbdir)

		caltables=glob.glob(caltable+"*.gcal")
		num_caltable=len(caltables)
		num_img=len(glob.glob(image+"*.image"))
		if num_caltable!=1:
			for i in range(0,num_caltable-1):
				os.system("rm -rf "+caltable+"_"+str(i).zfill(2)+".gcal")
		
		for i in range(0,num_img-1):
			for str1 in ['.image','.mask','.residual','.psf',\
					'.model','.sumwt','.pb','.fits']:
				os.system("rm -rf "+image+"_"+str(i).zfill(2)+str1)
		final_img=image+"_"+str(num_img-1).zfill(2)+".image"
		print (ms_in,current_trange)
		hf.imreg(vis=ms_in,imagefile=final_img,timerange=current_trange,fitsfile=final_img[:-6]+"_helio.fits",\
                            usephacenter=False,verbose=False) ### converting final image to heliocentric coordinates
		final_cal=caltbdir+"final_cal_spw_"+str(s).zfill(2)+".gcal"
		if os.path.isdir(final_cal)==False:
			os.system("mv "+caltable+"_"+str(num_caltable-1).zfill(2)+".gcal "+final_cal)
		else:
			gaincal(vis=msname, refant=refantenna, spw=sp,caltable=caltable, uvrange=uvrange, \
					calmode='p', combine='scan', minblperant=4, minsnr=3, append=True,solnorm=True)
		#print('Calibration done in {0:s}'.format(slfcalms))
		os.system("rm -rf "+slfcalms+"*")

selfcal_time=timeit.default_timer()
f.write("Time taken for selfcal in seconds is "+str(selfcal_time-start))
time1=selfcal_time		
if doapply==1:
	final_ms=ms_in.split('.')[0]+'_selfcalled.ms'
	print (final_ms, ms_in)
	print ("Applying caltables")
	if os.path.isdir(final_ms)==False:
            for s in selfcal_spws:
		#if s not in [3,6,12,24,27]:
		 #   continue
                caltable=caltbdir+'final_cal_spw_'+str(s).zfill(2)+'.gcal'
                if os.path.isdir(caltable)==False:
                    continue
                applycal(vis=ms_in,gaintable=caltable,spw=str(s),applymode='calonly',interp='nearest')
	#if os.path.isdir(final_ms)==True:
	#	os.system("rm -rf "+final_ms)
            split(vis=ms_in,outputvis=final_ms, correlation='XX',datacolumn='corrected')
            f.write("Calibrated MS split\n")
            time1=timeit.default_timer()
            f.write("Calibrated data split in seconds:"+str(time1-selfcal_time))


if doclean:
	import glob
	
	image_list=glob.glob(imgprefix+"*_helio.fits") 
	msname=final_ms
	xycen=get_img_center_heliocoords(image_list)
	
	imaging_start_mjd=peak_time_mjd-total_duration/2
	if imaging_start_mjd<startmjd:
		diff=peak_start_mjd-startmjd
		imaging_start_time=(flare_peak_time-dt.timedelta(seconds=diff)).strftime("%Y/%m/%d/%H:%M:%S")
	else:
		imaging_start_time=(flare_peak_time-dt.timedelta(seconds=total_duration/2)).strftime("%Y/%m/%d/%H:%M:%S")
	imaging_end_mjd=peak_time_mjd+total_duration/2
	if imaging_end_mjd>endmjd:
		diff=endmjd-peak_start_mjd
		imaging_end_time=(flare_peak_time+dt.timedelta(seconds=diff)).strftime("%Y/%m/%d/%H:%M:%S")
	else:
		imaging_end_time=(flare_peak_time+dt.timedelta(seconds=total_duration/2)).strftime("%Y/%m/%d/%H:%M:%S") 
	time_str=imaging_start_time+"~"+imaging_end_time
	
	specfile=msname[:-3]+"_dspec.npz"
	pol='XX'
	uvrange_ds='0.15~5.0km'


	spws=[]
	for s in selfcal_spws:
		spws.append(str(s))	
	cell_size1=[str(cell_size)+"arcsec"] 
	stokes='XX'
	mkmovie=True
	docompress=True	
	opencontour=False
	clevels=[0.8,1.0]
	plotaia=True
	aiawave=1600
	movieformat='mp4'
	overwrite=False
	subregion=''#box[[128pix,128pix],[284pix,384pix]]'
	qlookplot.qlookplot(vis=msname,specfile=specfile,timerange=time_str, spw=spws,\
                    ncpu=1,imsize=[final_imsize],cell=cell_size1, restoringbeam=[beam_1GHz],\
                    robust=0.5,opencontour=opencontour,clevels=clevels,plotaia=plotaia,\
                    aiawave=aiawave,mkmovie=mkmovie,twidth=final_image_int,docompress=docompress,\
                    stokes=stokes,movieformat=movieformat,uvrange='',\
                    niter=300,overwrite=overwrite,xycen=xycen,fov=[256,256])
	final_clean_time=timeit.default_timer()
	f.write("Final clean done in seconds:"+str(final_clean_time-time1))
f.close()
