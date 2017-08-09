import numpy as np, aipy as a, os
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import healpy as hp
from capo import omni

def quick_image(arr, cmap='jet'):
	##Quickly generate a 2D heatmap, akin to Bryna's quick_image IDL procedure.

	plt.imshow(arr,cmap=plt.get_cmap(cmap))
	plt.colorbar()
	plt.show()

def get_filelist(obsid):
        filelist=[
                'vis_data/'+obsid+"_vis_XX.sav",
                'vis_data/'+obsid+"_vis_YY.sav",
                'vis_data/'+obsid+"_autos.sav",
                'vis_data/'+obsid+"_flags.sav",
                'metadata/'+obsid+"_params.sav",
                obsid+"_settings.txt"
        ]
	return filelist

def quick_image_animate(arr, cmap='jet', axis=0, save=None, annotate='frame'):
	"""Animate a 2D heatmap along a given axis.
		axis = Which axis of the array is "time". The other two are x/y
		save = Save file name
                annotate = If string, 'frame' or 'none'. Frame = write frame number in top left. None = write nothing
                           If list, treat each entry as a string to annotate the corresponding frame
	"""
	if not len(arr.shape) == 3:
             print "Invalid array shape"
             return

	#Swap axes so that 0 is the time axis, if necessary.
        if not axis == 0:
		arr= np.swapaxes(arr,axis,0)
	fig = plt.figure()
        if isinstance(annotate, list):
            if not len(annotate) == arr.shape[0]: annotate='frame'    #Default to frame if incorrect shape
            else: frame_labels = map(str, annotate)
        if isinstance(annotate, str):
            if annotate=='frame': frame_labels=map(str,range(arr.shape[0]))
            if annotate=='none' : frame_labels=None

	im = plt.imshow(arr[0,:,:], cmap=cmap, interpolation='nearest')
        if not annotate is None:  text = plt.title(frame_labels[0], ha='center', va='center')
	def update_image(n):
		im.set_data(arr[n,:,:])
                if not annotate is None: text.set_text(frame_labels[n])

	ani = anim.FuncAnimation(fig, update_image, arr.shape[0], interval=1)
        if not save is None:
		writer = anim.writers['ffmpeg'](fps=20)
		ani.save(save,writer=writer, dpi=200)
	else:
		plt.show()


## Convert between healpix indices and ra/dec

def idx_to_radec(index,NSIDE):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)

def radec_to_index(ra,dec,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-dec+90.),np.radians(360.-RA))



def build_redundant(mode='sep', aa=None, cal=None, sep=-1, restore="nofile", length=None, min_count=2):
   if cal is None and restore is None:
        raise NameError("Either cal or restore, or both, must be defined.")
   if os.path.isfile(restore):
        cal = None
        fil=np.load(restore)
        reds=fil['reds']
        if 'lens' in fil.keys(): lens = fil['lens']
   if ((cal is not None) or (aa is not None)) and (not os.path.isfile(restore)):
        if aa is None: aa = a.cal.get_aa(cal, 1024/100e6, .1, 1024)
        info = omni.aa_to_info(aa)
        reds = info.get_reds()
        reds = [gp for gp in reds if len(gp) >= int(min_count)]
        lens = []
        for gp in reds:
                i,j = gp[0]
                bl = aa.get_baseline(i,j)   #baseline vector in light ns
                lens.append(np.sqrt(sum(i**2 for i in bl)) * a.const.len_ns / 100.)    #Convert to meters
        reds = np.array(reds)
        lens = np.array(lens)
        if mode == 'lengths':
             lens = np.unique(np.round(lens,3))
             return ','.join(lens)
        tolerance=0.50 if float(length)>0 else min(lens)+0.50     #50 cm good enough?
        if not length is None:
            ### Select only groups of baselines with a given length.
            reds=reds[np.where(np.abs(lens - float(length)) < tolerance)]
        if not restore is None and not restore == 'nofile':
              np.savez(restore,reds=reds,lens=lens)

   if mode == 'flatten':
	return ",".join([str(it[0])+"_"+str(it[1]) for sublist in reds for it in sublist])

   if mode == 'count':
	return len(reds)

   if mode == 'sep':
        if sep == -1: raise ValueError("Choose a valid redundant group index.")
	return [str(it[0])+"_"+str(it[1]) for it in reds[int(sep)]]
