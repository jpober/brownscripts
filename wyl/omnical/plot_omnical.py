import numpy as np
import sys, matplotlib.pyplot as plt

# usage: python plot_omnical obsid
obs=sys.argv[1]
name = { 0: '11', 1: '12', 2: '13', 3: '14', 4: '15', 5: '16', 6: '17', 7: '18', 8: '21', 9: '22', 10: '23',
11: '24', 12: '25', 13: '26', 14: '27', 15: '28', 16: '31', 17: '32', 18: '33', 19: '34', 20: '35', 21: '36',
22: '37', 23: '38', 24: '41', 25: '42', 26: '43', 27: '44', 28: '45', 29: '46', 30: '47', 31: '48', 32: '61',
33: '62', 34: '63', 35: '64', 36: '65', 37: '66', 38: '67', 39: '68', 40: '81', 41: '82', 42: '83', 43: '84',
44: '85', 45: '86', 46: '87', 47: '88', 48: '91', 49: '92', 50: '93', 51: '94', 52: '95', 53: '96', 54: '97',
55: '98', 56: '1001', 57: '1002', 58: '1003', 59: '1004', 60: '1005', 61: '1006', 62: '1007', 63: '1008',
64: '1009', 65: '1010', 66: '1011', 67: '1012', 68: '1013', 69: '1014', 70: '1015', 71: '1016', 72: '1017',
73: '1018', 74: '1019', 75: '1020', 76: '1021', 77: '1022', 78: '1023', 79: '1024', 80: '1025', 81: '1026',
82: '1027', 83: '1028', 84: '1029', 85: '1030', 86: '1031', 87: '1032', 88: '1033', 89: '1034', 90: '1035',
91: '1036', 92: '1037', 93: '1038', 94: '1039', 95: '1040', 96: '1041', 97: '1042', 98: '1043', 99: '1044',
100: '1045', 101: '1046', 102: '1047', 103: '1048', 104: '1049', 105: '1050', 106: '1051', 107: '1052',
108: '1053', 109: '1054', 110: '1055', 111: '1056', 112: '1057', 113: '1058', 114: '1059', 115: '1060',
116: '1061', 117: '1062', 118: '1063', 119: '1064', 120: '1065', 121: '1066', 122: '1067', 123: '1068',
124: '1069', 125: '1070', 126: '1071', 127: '1072'}

dx=np.load(obs+'.xx.npz')
dy=np.load(obs+'.yy.npz')
fm=np.zeros((56,384),dtype=bool)
badf = [0,1,2,13,14,15]
for nn in range(0,384):
    if nn%16 in badf: fm[:,nn]=True
freq=dx['freqs']
SH=fm.shape
sol={}
ampmax=1
ampmin=1
count=0
refphsx = dx['56x']/np.abs(dx['56x'])
refphsy = dy['56y']/np.abs(dy['56y'])
for ii in range(0,128):
    sol[ii]={}
    try: 
        x=dx[str(ii)+'x']/refphsx
        y=dy[str(ii)+'y']/refphsy
        exist=True
    except(KeyError):
        x=np.ones(SH)
        y=np.ones(SH)
        x*=np.nan
        y*=np.nan
        exist=False
    mx=np.ma.masked_array(x,mask=fm)
    mx=np.mean(mx,axis=0)
    my=np.ma.masked_array(y,mask=fm)
    my=np.mean(my,axis=0)
    ddx=mx.data
    ddy=my.data
    ff=mx.mask
    for jj in range(0,ff.size):
        if ff[jj]:
            ddx[jj]=np.nan
            ddy[jj]=np.nan
    if exist:
        dxmax=np.nanmax((np.abs(ddx)))
        dymax=np.nanmax((np.abs(ddy)))
        dxmin=np.nanmin((np.abs(ddx)))
        dymin=np.nanmin((np.abs(ddy)))
#        print dxmax,dymax,dxmin,dymin
        if count==0:
            ampmax=dxmax
            ampmin=dxmin
            count+=1
        if dxmax>ampmax: ampmax=dxmax
        if dymax>ampmax: ampmax=dymax
        if dxmin<ampmin: ampmin=dxmin
        if dymin<ampmin: ampmin=dymin
    sol[ii]['x']=ddx
    sol[ii]['y']=ddy

#plot and save amplitude
ly=[ampmin,(ampmax+ampmin)/2,ampmax]
lay=['%.2f'%(ampmin),'%.2f'%((ampmax+ampmin)/2),'%.2f'%(ampmax)]
lx=[37,287]
lax=[int(round(freq[37]/1e6)),int(round(freq[287]/1e6))]
fx=np.arange(0,384)
fig=plt.figure()
plt.suptitle(obs.split('/')[-1],y=0.99,size=15.0)
for ii in range(0,8):
    for jj in range(0,16):
        ind=ii*16+jj
        p=fig.add_subplot(8,16,ind+1)
        p.scatter(fx,np.abs(sol[ind]['x']),color='blue',s=0.01)
        p.scatter(fx,np.abs(sol[ind]['y']),color='red',s=0.01)
        plt.ylim((ampmin,ampmax))
        plt.xlim((0,383))
        p.set_xticks(lx)
        p.set_yticks(ly)
        if jj==0: p.yaxis.set_ticklabels(lay,size=6.5)
        else: p.yaxis.set_ticklabels([])
        if ii==7: p.xaxis.set_ticklabels(lax,size=6.5)
        else: p.xaxis.set_ticklabels([])
        p.set_title(name[ind],size=6.5,y=0.9)
plt.subplots_adjust(top=0.93,bottom=0.03,left=0.04,right=0.98)
fig.savefig(obs+'_amp_omnical.png')

#plot and save phase
pi=3.14159265358979323846
ly=[-pi,0,pi]
lay=['-3.14','0','3.14']
fig=plt.figure()
plt.suptitle(obs.split('/')[-1],y=0.99,size=15.0)
for ii in range(0,8):
    for jj in range(0,16):
        ind=ii*16+jj
        p=fig.add_subplot(8,16,ind+1)
        p.scatter(fx,np.angle(sol[ind]['x']),color='blue',s=0.01)
        p.scatter(fx,np.angle(sol[ind]['y']),color='red',s=0.01)
        plt.ylim((-pi,pi))
        plt.xlim((0,383))
        p.set_xticks(lx)
        p.set_yticks(ly)
        if jj==0: p.yaxis.set_ticklabels(lay,size=6.5)
        else: p.yaxis.set_ticklabels([])
        if ii==7: p.xaxis.set_ticklabels(lax,size=6.5)
        else: p.xaxis.set_ticklabels([])
        p.set_title(name[ind],size=6.5,y=0.9)
plt.subplots_adjust(top=0.93,bottom=0.03,left=0.04,right=0.98)
fig.savefig(obs+'_phs_omnical.png')
 
        

