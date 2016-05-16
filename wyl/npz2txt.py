#moved this function to wenyang-li/capo/src/omni.py

import numpy as np

def writetxt(npzfiles):
    
    p2pol = {'EE': 'x','NN': 'y','EN': 'cross', 'NE': 'cross'}  #check the convension
    
    #create output file
    fn0 = npzfiles[0].split('.')
    fn0[-1] = 'txt'
    outfn = '.'.join(fn0)
    outfile = open(outfn,'w')
    outfile.write("# Program of origin: Omnical\n")
    outfile.write("# Convention: Divide uncalibrated data by these gains to obtain calibrated data.\n")
    outfile.write("# ANT NAME, ANT INDEX, FREQ (MHZ), POL, TIME (JD), RE(GAIN), IM(GAIN), FLAG\n")
    
    #read gain solutions from npz
    
    #npzdict = {}
    for f,filename in enumerate(npzfiles):
        data = np.load(filename)
        ant = []
        for ii, ss in enumerate(data):
            if ss[0].isdigit():
                intss = int(ss[0:-1])
                if not intss in ant:
                    ant.append(intss)
        time = data['jds']
        freq = data['freqs']/1e6
        pol = ['EE', 'NN', 'EN', 'NE']
        nt = time.shape[0]
        nf = freq.shape[0]
        na = len(ant)
        for tt in range(0, nt):
            for pp in range(0, 4):
                for ff in range(0, nf):
                    for iaa in range(0, na):
                        aa = ant[iaa]
                        dt = time[tt]
                        dp = pol[pp]
                        df = freq[ff]
                        stkey = str(aa) + p2pol[pol[pp]]
                        try: da = data[stkey][tt][ff]
                        except: da = 1.0
                        string = 'ant'+str(aa)+', '+str(aa)+', '+str(df)+', '+dp+', '+str(dt)+', '+str(da.real)+', '+str(da.imag)+', 0\n'
                        outfile.write(string)
    outfile.close()


