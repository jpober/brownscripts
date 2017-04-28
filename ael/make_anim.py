#!/usr/bin/env python

import sys, re, glob
import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.animation as anim
#from moviepy.editor import VideoClip

### Load the latest ffmpeg module before using this!!!

files = sys.argv[1:]  #glob.glob('*[0-17]*Dirty_XX.png')

frames = []

for f in files:
     print f
     img = mpimg.imread(f)
     frames.append(img)


fig = plt.figure()
im = plt.imshow(frames[0], cmap='brg_r', interpolation='nearest')
def update_image(n):
	im.set_data(frames[n])
ani = anim.FuncAnimation(fig, update_image, len(frames), interval=5)
plt.show()

#writer = anim.writers['ffmpeg'](fps=20)
#ani.save('anim.mp4',writer=writer, dpi=200)

