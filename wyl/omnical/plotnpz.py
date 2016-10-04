#wenyang li, May 24, 2016
#usage:
#For gain solution: python plotnpz.py npz.npz index pol
#                   eg: python plotnpz.py npz.npz 0 x
#For model vis: python plotnpz.py npz.npz index1 index2 polpol
#               eg: python plotnpz.py npz.npz 0 1 xx


import numpy as np, sys, matplotlib.pyplot as plt

length = len(sys.argv)
if length > 2:
    npz = sys.argv[1]
    if length == 4:
        string = sys.argv[2]+sys.argv[3]
    else:
        string = '<' + sys.argv[2] + ',' + sys.argv[3] + '> ' + sys.argv[4]
    
    data = np.load(npz)
    dd = data[string]
    
    fig=plt.figure()
    p1=fig.add_subplot(1,2,1)
    p1.set_title("real")
    i1=p1.imshow(dd.real, interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    p2=fig.add_subplot(1,2,2)
    p2.set_title("imag")
    i2=p2.imshow(dd.imag, interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    plt.show()