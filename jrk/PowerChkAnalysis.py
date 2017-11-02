import numpy as n
from glob import glob
import pylab as pl

pwrtxts = glob('/users/jkerriga/data/jkerriga/AllSepsLST/even/*npy')
#print pwrtxts
#ptxt1 = glob('./even/Pzen*txt')
#ptxt2 = glob('./odd/Pzen*txt')
#pwrtxts=ptxt1+ptxt2
dirtyP = {}
resP = {}
lst = []
#chans = ['CH03','CH14','CH25']

blines = (30.0,40.0,60.0,70.0,90.0,120.0,150.0,180.0,210.0)
for init in blines:
    dirtyP[init] = []
    resP[init] = []
#subP = []
## Ch1:53_73 Ch2:73_93 Ch3:93_113
outpwr = []
suffix = 'uvcRREcACOTUc'
for i in pwrtxts:
    print i
    obsP = n.load(i)
    for j in blines:
        dirtyP[j].append(obsP.item()['dirty'][j])
        resP[j].append(obsP.item()['res'][j])
    #dirtyP.append(obsP[0])
    #resP.append(obsP[1])
    lst.append(obsP.item()['dirty']['time'])
    #outpwr.append('Pzen.'+i.split('npy')[0]+suffix)
#print outpwr
numGood = 0
maxGood = []
#for r in blines:
    #subP = n.subtract(n.abs(dirtyP[r]),n.abs(resP[r]))
    #divP = resP[r]/dirtyP[r]
    #if (subP>0).sum() > numGood:
    #    maxGood = subP
    #    numGood = (subP>0).sum()
    #print (subP>0).sum()
    #print n.size(subP)
   # print r
#print 'Power Removed: ',subP[subP>0].sum()
#print 'Power Added: ', subP[subP<0].sum()
lst = n.array(lst)
#dirtyP = n.array(dirtyP)
#resP = n.array(resP)
#pwrtxts = n.array(pwrtxts)
#badIdx= n.where(subP<0.)
#print n.shape(badIdx[0])
#print pwrtxts[subP<0]
x = pl.histogram(lst,32)
#summed = n.zeros(len(x[1])-1)

lstMean = []
totPower=n.zeros((len(x[1])-1,len(blines)),dtype=complex)
resPower=n.zeros((len(x[1])-1,len(blines)),dtype=complex)
divPower=n.zeros((len(x[1])-1,len(blines)))
for i in range(0,len(x[1])-1):
    lstMean.append((x[1][i]+x[1][i+1])/2.0)
    for j in range(len(blines)):
        b = blines[j]
        arr1=lst>x[1][i]
        arr2=lst<x[1][i+1]
        msk = ((arr1.astype(int) +  arr2.astype(int))/2).astype(bool)
        msk = n.array(msk)
        dirtyP[b] = n.array(dirtyP[b])
        resP[b] = n.array(resP[b])
        totPower[i,j]=dirtyP[b][msk].sum()
        resPower[i,j]=resP[b][msk].sum()
        print dirtyP[b][msk].sum(),resP[b][msk].sum()
        divPower[i,j]=n.divide(resP[b][msk].sum(),dirtyP[b][msk].sum())
    #summed[i]=(subP[msk]).sum()
maxGood = n.array(maxGood)
outpwr = n.array(outpwr)
#print outpwr.shape,subP.shape
#n.savetxt('GoodObs.txt',outpwr[maxGood[:,0]>0],fmt='%s')
#n.savetxt('BadObs.txt',outpwr[maxGood[:,0]<0],fmt='%s')
#print n.mean(subP)
#pl.plot(subP,'.')
#print n.shape(lstMean)
#print resPower
resPower = n.array(n.absolute(resPower))
#resPower = n.absolute(resPower-totPower)
totPower = n.array(n.absolute(totPower))
divPower = n.array(divPower)
#### Fornax is at LST 
lstMean = n.array(lstMean)
print lstMean[lstMean<6.24]
#print divPower
ax1 = pl.subplot(111)
ax1.plot(lstMean,1.-divPower[:,1])
#ax1.close()
for i in range(9):

    ax = pl.subplot(3,3,i+1,sharex=ax1)
    ax.plot(lstMean,1.-divPower[:,i],color='k')
    ### Fornax center of FoV
    ax.axvline(3.67,0,1.1,color='r',linewidth=2.0,linestyle='dashed')
    ### Pictor center FoV
    ax.axvline(5.316,0,1.1,color='b',linewidth=2.0,linestyle='dashed')
    fornax = lstMean[lstMean<5.82]
    fornax = fornax[fornax>1.26]
    ax.fill_between(fornax,0,1.1,facecolor='red',alpha=0.4)
    pictor = lstMean[lstMean<7.2]
    pictor = pictor[pictor>3.9]
    ax.fill_between(pictor,0,1.1,facecolor='blue',alpha=0.4)

    ax.set_ylim(0.0,1.05)
    ax.set_xlim(0.5,8.0)
    pl.title(str(blines[i])+'m')
pl.tight_layout()
#pl.show()
pl.savefig('RatioPerBaseline.png')   
