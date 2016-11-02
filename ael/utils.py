import numpy as np
import matplotlib.pyplot as plt

def quick_image(arr, cmap='jet'):
	##Quickly generate a 2D heatmap, akin to Bryna's quick_image IDL procedure.

	plt.imshow(arr,cmap=plt.get_cmap(cmap))
	plt.colorbar()
	plt.show()
