import numpy as np
import matplotlib.pyplot as plt

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
