#run the following in command line:
#python npz2txt npzfiles 
#eg: python npz2txt zen.jds.xx.npz zen.jds.yy.npz
#It outputs a txt file zen.jds.txt

import numpy, sys

length=len(sys.argv)

p2pol={'x': 'EE', 'y':'NN'}

if length>1:
    
    #create output file
    args=sys.argv[1:]
    fn0=sys.argv[1].split('.')
    lenfn0=len(fn0)
    outfn=''
    for ii in range(0,lenfn0-3):
        outfn+=(fn0[ii]+'.')
    outfn+='npz.txt'
    outfile=open(outfn,'w') 
    outfile.write("# Program of origin: RTS\n")
    outfile.write("# Convention: Divide uncalibrated data by these gains to obtain calibrated data.\n")
    outfile.write("# ANT NAME, ANT INDEX, FREQ (MHZ), POL, TIME (JD), RE(GAIN), IM(GAIN), FLAG\n")
    
    #read gain solutions from npz
    npzdict={}
    for f,filename in enumerate(args):
        data=numpy.load(filename)
        ant=[]
        for ii, ss in enumerate(data):
            if ss[0].isdigit():
                ant.append(ss)
        time=data['jds']
        freq=data['freqs']
        pol=p2pol[ant[0][-1]]
        for ti, tt in enumerate(time):
            if not npzdict.has_key(tt):
                npzdict[tt]={}
            if not npzdict.has_key(pol):
                npzdict[tt][pol]={}
            for fi, ff in enumerate(freq):
                if not npzdict[tt][pol].has_key(ff):
                    npzdict[tt][pol][ff]={}
                for ai, aa in enumerate(ant):
                    aa0=int(aa[:-1])
                    if not npzdict[tt][pol][ff].has_key(aa0):
                        npzdict[tt][pol][ff][aa0]=data[aa][ti][fi]
    nullpol=['EN', 'NE']
    if not npzdict[tt].has_key('NN'):
        nullpol.append('NN')
    if not npzdict[tt].has_key('EE'):
        nullpol.append('EE')
    for pi,pp in enumerate(nullpol):
        for ti, tt in enumerate(time):
            npzdict[tt][pp]={}
            for fi, ff in enumerate(freq):
                npzdict[tt][pp][ff]={}
                for ai, aa in enumerate(ant):
                    aa0=int(aa[:-1])
                    npzdict[tt][pp][ff][aa0]=complex(1.0,0.0)
    
    #write to txt
    for tt in npzdict:
        for pp in npzdict[tt]:
            for ff in npzdict[tt][pp]:
                for aa0 in npzdict[tt][pp][ff]:
                    string='ant'+str(aa0)+', '+str(aa0)+', '+str(ff)+', '+pp+', '+str(tt)+', '+str(npzdict[tt][pp][ff][aa0].real)+', '+str(npzdict[tt][pp][ff][aa0].imag)+', 0\n'
                    outfile.write(string)
