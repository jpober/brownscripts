import aipy as a, numpy as n


class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
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


prms = {

#    'loc': ('-30:43:17.5', '21:25:41.9'), # Karoo, ZAR, GPS. elevation=1085m
    'loc': ('-26:42:11.916', '116:40:15.06'), # MWA,  (GPS)  -- For now, pretend it's at the MWA site
    'antpos': {
        0: {'top_y': -113.79573805727523, 'top_x': -80.3, 'top_z': 0.0},
        1: {'top_y': -113.79573805727523, 'top_x': -65.7, 'top_z': 0.0},
        2: {'top_y': -113.79573805727523, 'top_x': -51.1, 'top_z': 0.0},
        3: {'top_y': -113.79573805727523, 'top_x': -36.5, 'top_z': 0.0},
        4: {'top_y': -101.15176716202242, 'top_x': -87.6, 'top_z': 0.0},
        5: {'top_y': -101.15176716202242, 'top_x': -73.0, 'top_z': 0.0},
        6: {'top_y': -101.15176716202242, 'top_x': -58.4, 'top_z': 0.0},
        7: {'top_y': -101.15176716202242, 'top_x': -43.8, 'top_z': 0.0},
        8: {'top_y': -101.15176716202242, 'top_x': -29.2, 'top_z': 0.0},
        9: {'top_y': -88.50779626676963, 'top_x': -94.9, 'top_z': 0.0},
        10: {'top_y': -88.50779626676963, 'top_x': -80.3, 'top_z': 0.0},
        11: {'top_y': -88.50779626676963, 'top_x': -65.7, 'top_z': 0.0},
        12: {'top_y': -88.50779626676963, 'top_x': -51.1, 'top_z': 0.0},
        13: {'top_y': -88.50779626676963, 'top_x': -36.5, 'top_z': 0.0},
        14: {'top_y': -88.50779626676963, 'top_x': -21.9, 'top_z': 0.0},
        15: {'top_y': -75.86382537151681, 'top_x': -102.2, 'top_z': 0.0},
        16: {'top_y': -75.86382537151681, 'top_x': -87.6, 'top_z': 0.0},
        17: {'top_y': -75.86382537151681, 'top_x': -73.0, 'top_z': 0.0},
        18: {'top_y': -75.86382537151681, 'top_x': -58.4, 'top_z': 0.0},
        19: {'top_y': -75.86382537151681, 'top_x': -43.8, 'top_z': 0.0},
        20: {'top_y': -75.86382537151681, 'top_x': -29.2, 'top_z': 0.0},
        21: {'top_y': -75.86382537151681, 'top_x': -14.6, 'top_z': 0.0},
        22: {'top_y': -63.21985447626402, 'top_x': -94.9, 'top_z': 0.0},
        23: {'top_y': -63.21985447626402, 'top_x': -80.3, 'top_z': 0.0},
        24: {'top_y': -63.21985447626402, 'top_x': -65.7, 'top_z': 0.0},
        25: {'top_y': -63.21985447626402, 'top_x': -51.1, 'top_z': 0.0},
        26: {'top_y': -63.21985447626402, 'top_x': -36.5, 'top_z': 0.0},
        27: {'top_y': -63.21985447626402, 'top_x': -21.9, 'top_z': 0.0},
        28: {'top_y': -50.57588358101121, 'top_x': -87.6, 'top_z': 0.0},
        29: {'top_y': -50.57588358101121, 'top_x': -73.0, 'top_z': 0.0},
        30: {'top_y': -50.57588358101121, 'top_x': -58.4, 'top_z': 0.0},
        31: {'top_y': -50.57588358101121, 'top_x': -43.8, 'top_z': 0.0},
        32: {'top_y': -50.57588358101121, 'top_x': -29.2, 'top_z': 0.0},
        33: {'top_y': -37.931912685758405, 'top_x': -80.3, 'top_z': 0.0},
        34: {'top_y': -37.931912685758405, 'top_x': -65.7, 'top_z': 0.0},
        35: {'top_y': -37.931912685758405, 'top_x': -51.1, 'top_z': 0.0},
        36: {'top_y': -37.931912685758405, 'top_x': -36.5, 'top_z': 0.0}
	}
}



def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    phsoff = {'x':[0.,0.], 'y':[0.,0.]}
    nants = len(prms['antpos'])
    for i in range(nants):
        bp_r = {'x':n.array([0.0]*9), 'y':n.array([0.0]*9)}
        bp_i = {'x':n.array([0.0]*3), 'y':n.array([0.0]*3)}
        antennas.append(a.pol.Antenna(0,0,0, a.fit.Beam(freqs),
             amp={'x':1.,'y':1.},bp_r=bp_r, bp_i=bp_i, phsoff=phsoff))
    aa = AntennaArray(prms['loc'], antennas)#, ant_layout=prms['ant_layout'])
#    for i in range(nants):
#        pos = prms['antpos-top'][i]
#        i = str(i)
#        aa.set_params({i:{'top_x':pos[0], 'top_y':pos[1], 'top_z':pos[2]}})
    for i in range(nants):
        pos = prms['antpos'][i]
        i = str(i)
        aa.set_params({i:pos})
    return aa

