#wenyang li, May 24, 2016
#usage:
#For gain solution: python comparenpz npz1.npz npz2.npz index pol
#                   eg: python comparenpz npz1.npz npz2.npz 0 x
#For model vis: python comparenpz npz1.npz npz2.npz index1 index2 polpol
#               eg: python comparenpz npz1.npz npz2.npz 0 1 xx

import numpy as np, sys, matplotlib.pyplot as plt

length = len(sys.argv)

if length > 4:
    npz0 = sys.argv[1]
    npz1 = sys.argv[2]
    if length == 5:
        string = sys.argv[3] + sys.argv[4]
    else:
        string = '<' + sys.argv[3] + ',' + sys.argv[4] + '> ' + sys.argv[5]
    
    data0 = np.load(npz0)
    data1 = np.load(npz1)
    
    dd0 = data0[string]
    dd1 = data1[string]
    diff = 2*(dd0 - dd1)/(np.abs(dd0)+np.abs(dd1))

    fig=plt.figure()
    p1=fig.add_subplot(3,2,1)
    p1.set_title("solution from npz1 (real)")
    i1=p1.imshow(dd0.real, interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    p2=fig.add_subplot(3,2,2)
    p2.set_title("solution from npz1 (imag)")
    i2=p2.imshow(dd0.imag, interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    p3=fig.add_subplot(3,2,3)
    p3.set_title("solution from npz2 (real)")
    i3=p3.imshow(dd1.real, interpolation='nearest',aspect='auto')
    fig.colorbar(i3)
    p4=fig.add_subplot(3,2,4)
    p4.set_title("solution from npz2 (imag)")
    i4=p4.imshow(dd1.imag, interpolation='nearest',aspect='auto')
    fig.colorbar(i4)
    p5=fig.add_subplot(3,2,5)
    p5.set_title("fractional difference (real)")
    i5=p5.imshow(diff.real, interpolation='nearest',aspect='auto')
    fig.colorbar(i5)
    p6=fig.add_subplot(3,2,6)
    p6.set_title("fractional difference (imag)")
    i6=p6.imshow(diff.imag, interpolation='nearest',aspect='auto')
    fig.colorbar(i6)


    plt.show()





