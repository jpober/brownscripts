import aipy as a, numpy as n


class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
        self.antpos_ideal = kwargs.pop('antpos_ideal')
#        self.ant_layout = kwargs.pop('ant_layout')
    def update(self):
        a.pol.AntennaArray.update(self)
    def get_params(self, ant_prms={'*':'*'}):
        try: prms = a.pol.AntennaArray.get_params(self, ant_prms)
        except(IndexError): return {}
        for k in ant_prms:
            if k == 'aa':
                if not prms.has_key('aa'): prms['aa'] = {}
                for val in ant_prms[k]:
                    if   val == 'tau_ns': prms['aa']['tau_ns'] = self.tau_ns
                    elif val == 'tau_ew': prms['aa']['tau_ew'] = self.tau_ew
                    elif val == 'gain': prms['aa']['gain'] = self.gain
            else:
                try: top_pos = n.dot(self._eq2zen, self[int(k)].pos)
                # XXX should multiply this by len_ns to match set_params.
                except(ValueError): continue
                if ant_prms[k] == '*':
                    prms[k].update({'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]})
                else:
                    for val in ant_prms[k]:
                        if   val == 'top_x': prms[k]['top_x'] = top_pos[0]
                        elif val == 'top_y': prms[k]['top_y'] = top_pos[1]
                        elif val == 'top_z': prms[k]['top_z'] = top_pos[2]
        return prms
    def set_params(self, prms):
        changed = a.pol.AntennaArray.set_params(self, prms)
        for i, ant in enumerate(self):
            ant_changed = False
            top_pos = n.dot(self._eq2zen, ant.pos)
            try:
                top_pos[0] = prms[str(i)]['top_x']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[1] = prms[str(i)]['top_y']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[2] = prms[str(i)]['top_z']
                ant_changed = True
            except(KeyError): pass
            if ant_changed: ant.pos = n.dot(n.linalg.inv(self._eq2zen), top_pos) / a.const.len_ns * 100.
            changed |= ant_changed
        if changed: self.update()
        return changed



### The following line gets found and replaced 
INSERT_PRMS_HERE




def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    phsoff = {'x':[0.,0.], 'y':[0.,0.]}
    if not 'antpos' in prms:
        prms['antpos'] = prms['antpos_ideal']
    nants = len(prms['antpos_ideal'])
    antpos_ideal = n.zeros(shape=(nants,3),dtype=float)
    tops = {'top_x':0, 'top_y':1, 'top_z':2}
    for k in prms['antpos_ideal'].keys():
        for i,m in enumerate(prms['antpos_ideal'][k]):
            antpos_ideal[k,tops[m]] = prms['antpos_ideal'][k][m]
    for i in range(nants):
        bp_r = {'x':n.array([0.0]*9), 'y':n.array([0.0]*9)}
        bp_i = {'x':n.array([0.0]*3), 'y':n.array([0.0]*3)}
        antennas.append(a.pol.Antenna(0,0,0, a.fit.Beam(freqs),
             amp={'x':1.,'y':1.},bp_r=bp_r, bp_i=bp_i, phsoff=phsoff))
    aa = AntennaArray(prms['loc'], antennas, antpos_ideal=antpos_ideal)#, ant_layout=prms['ant_layout'])
#    for i in range(nants):
#        pos = prms['antpos-top'][i]
#        i = str(i)
#        aa.set_params({i:{'top_x':pos[0], 'top_y':pos[1], 'top_z':pos[2]}})
    for i in range(nants):
        pos = prms['antpos_ideal'][i]
        i = str(i)
        aa.set_params({i:pos})
    return aa

