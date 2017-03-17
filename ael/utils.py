import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

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

def quick_image_animate(arr, cmap='jet', axis=0):
	##Animate a 2D heatmap along a third axis, 
	return 0


## Convert between healpix indices and ra/dec

def idx_to_radec(index,NSIDE):
    theta,phi=hp.pixelfunc.pix2ang(NSIDE,index)
    return -np.degrees(theta-np.pi/2.),np.degrees(np.pi*2.-phi)

def radec_to_index(ra,dec,NSIDE):
    return hp.pixelfunc.ang2pix(NSIDE,np.radians(-dec+90.),np.radians(360.-RA))
