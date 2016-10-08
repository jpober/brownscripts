### pulled some code from https://github.com/HERA-Team/pyuvdata.git to read in data only ###

import uvdata
import numpy as np
from astropy.io import fits
import aipy as a
from scipy.io.idl import readsav
import warnings

def find_ex_ant(uvdata):
    ex_ant = []
    for ii in uvdata.antenna_numbers:
        if not ii in uvdata.ant_1_array and not ii in uvdata.ant_2_array:
            ex_ant.append(ii)
    return ex_ant
    
class data_uvfits(uvdata.uvfits.UVFITS):
    def read_data_only(self, filename):
        F = fits.open(filename)
        D = F[0]
        hdunames = self._indexhdus(F)
        time0_array = D.data['DATE']
        try:
            time1_array = D.data['_DATE']
            self.time_array = (time0_array.astype(np.double) + time1_array.astype(np.double))
        except(KeyError):
            self.time_array = time0_array
        
        self.Ntimes = len(np.unique(self.time_array))
        
        try:
            self.ant_1_array = np.int32(D.data.field('ANTENNA1')) - 1
            self.ant_2_array = np.int32(D.data.field('ANTENNA2')) - 1
        except(KeyError):
            bl_input_array = np.int64(D.data.field('BASELINE'))
            self.ant_1_array, self.ant_2_array = self.baseline_to_antnums(bl_input_array)
        self.Nblts = len(self.time_array)
        self.Nbls = self.Nblts / self.Ntimes
        
        if D.header['NAXIS'] == 7:
            self.data_array = (D.data.field('DATA')[:, 0, 0, :, :, :, 0] + 1j * D.data.field('DATA')[:, 0, 0, :, :, :, 1])
            self.flag_array = (D.data.field('DATA')[:, 0, 0, :, :, :, 2] <= 0)
        else:
            self.data_array = (D.data.field('DATA')[:, 0, 0, :, :, 0] + 1j * D.data.field('DATA')[:, 0, 0, :, :, 1])
            self.data_array = self.data_array[:, np.newaxis, :, :]
            self.flag_array = (D.data.field('DATA')[:, 0, 0, :, :, 2] <= 0)
            self.flag_array = self.flag_array[:, np.newaxis, :, :]
        self.Nfreqs = D.header['NAXIS4']
        self.Npols = D.header['NAXIS3']
        self.freq_array = self._gethduaxis(D, 4)
        self.polarization_array = np.int32(self._gethduaxis(D, 3))
        
        ant_hdu = F[hdunames['AIPS AN']]
        self.antenna_names = ant_hdu.data.field('ANNAME').tolist()
        self.antenna_numbers = ant_hdu.data.field('NOSTA') - 1
        self.Nants_telescope = len(self.antenna_numbers)
        
class data_miriad(uvdata.miriad.Miriad):
    def read_data_only(self, filename):
        uv = a.miriad.UV(filename)
        miriad_header_data = {'Nfreqs': 'nchan', 'Npols': 'npol', 'Nants_telescope':'nants', 'channel_width': 'sdf',}
        for item in miriad_header_data:
            if isinstance(uv[miriad_header_data[item]], str):
                header_value = uv[miriad_header_data[item]].replace('\x00', '')
            else:
                header_value = uv[miriad_header_data[item]]
            setattr(self, item, header_value)
        data_accumulator = {}
        for (uvw, t, (i,j)), d, f in uv.all(raw=True):
            try: data_accumulator[uv['pol']].append([t,i,j,d,f])
            except(KeyError): data_accumulator[uv['pol']]=[[t,i,j,d,f]]
        self.polarization_array = np.sort(data_accumulator.keys())
        times = list(set(np.concatenate([[k[0] for k in d] for d in data_accumulator.values()])))
        times = np.sort(times)
        ant_i_unique = list(set(np.concatenate([[k[1] for k in d] for d in data_accumulator.values()])))
        ant_j_unique = list(set(np.concatenate([[k[2] for k in d] for d in data_accumulator.values()])))
        self.antenna_numbers = np.arange(self.Nants_telescope)
        self.antenna_names = self.antenna_numbers.astype(str).tolist()
        t_grid = []
        ant_i_grid = []
        ant_j_grid = []
        for ant_i in ant_i_unique:
            for ant_j in ant_j_unique:
                if ant_i > ant_j: continue
                ant_i_grid.append(ant_i)
                ant_j_grid.append(ant_j)
        self.Nbls = len(ant_i_grid)
        self.Ntimes = len(times)
        ant_i_grid *= self.Ntimes
        ant_j_grid *= self.Ntimes
        for t in times:
            t_grid += [t] * self.Nbls
        ant_i_grid = np.array(ant_i_grid)
        ant_j_grid = np.array(ant_j_grid)
        t_grid = np.array(t_grid)        
        self.Nblts = len(t_grid)
        self.time_array = t_grid
        self.ant_1_array = ant_i_grid
        self.ant_2_array = ant_j_grid
        self.data_array = np.zeros((self.Nblts, 1, self.Nfreqs, self.Npols), dtype=np.complex64)
        self.flag_array = np.ones(self.data_array.shape, dtype=np.bool)
        self.freq_array = (np.arange(self.Nfreqs) * self.channel_width + uv['sfreq'] * 1e9)
        self.freq_array = np.tile(self.freq_array, (1, 1))
        for pol, data in data_accumulator.iteritems():
            pol_ind = self._pol_to_ind(pol)
            for ind, d in enumerate(data):
                t, ant_i, ant_j = d[0], d[1], d[2]
                blt_index = np.where(np.logical_and(np.logical_and(t == t_grid, ant_i == ant_i_grid),ant_j == ant_j_grid))[0].squeeze()
                self.data_array[blt_index, :, :, pol_ind] = d[3]
                self.flag_array[blt_index, :, :, pol_ind] = d[4]
                
class data_fhd(uvdata.fhd.FHD):
    def read_data_only(self,filelist,use_model=False):
        datafiles = {}
        params_file = None
        flags_file = None
        settings_file = None
        data_name = '_vis_'
        if use_model: data_name = '_vis_model_'
        for file in filelist:
            if file.lower().endswith(data_name + 'xx.sav'):
                datafiles['xx'] = file
            elif file.lower().endswith(data_name + 'yy.sav'):
                datafiles['yy'] = file
            elif file.lower().endswith(data_name + 'xy.sav'):
                datafiles['xy'] = file
            elif file.lower().endswith(data_name + 'yx.sav'):
                datafiles['yx'] = file
            elif file.lower().endswith('_params.sav'):
                params_file = file
            elif file.lower().endswith('_flags.sav'):
                flags_file = file
            elif file.lower().endswith('_settings.txt'):
                settings_file = file
            else:
                continue
        if len(datafiles) < 1:
            raise StandardError('No data files included in file list')
        if params_file is None:
            raise StandardError('No params file included in file list')
        if flags_file is None:
            raise StandardError('No flags file included in file list')
        if settings_file is None:
            warnings.warn('No settings file included in file list')
        vis_data = {}
        for pol, file in datafiles.iteritems():
            this_dict = readsav(file, python_dict=True)
            if use_model:
                vis_data[pol] = this_dict['vis_model_ptr']
            else:
                vis_data[pol] = this_dict['vis_ptr']
            this_obs = this_dict['obs']
            data_shape = vis_data[pol].shape
        obs = this_obs
        bl_info = obs['BASELINE_INFO'][0]
        self.Ntimes = int(obs['N_TIME'][0])
        self.Nbls = int(obs['NBASELINES'][0])
        self.Nblts = data_shape[0]
        self.Nfreqs = int(obs['N_FREQ'][0])
        self.Npols = len(vis_data.keys())
        fhd_pol_list = []
        for pol in obs['POL_NAMES'][0]:
            fhd_pol_list.append(pol.decode("utf-8").lower())
        
        #params_dict = readsav(params_file, python_dict=True)
        #params = params_dict['params']
        flag_file_dict = readsav(flags_file, python_dict=True)
        
        vis_weights_data = {}
        if 'flag_arr' in flag_file_dict:
            weights_key = 'flag_arr'
        elif 'vis_weights' in flag_file_dict:
            weights_key = 'vis_weights'
        else:
            raise ValueError('No recognized key for visibility weights in flags_file.')
        for index, w in enumerate(flag_file_dict[weights_key]):
            vis_weights_data[fhd_pol_list[index]] = w
        
        lin_pol_order = ['xx', 'yy', 'xy', 'yx']
        linear_pol_dict = dict(zip(lin_pol_order, np.arange(5, 9) * -1))
        pol_list = []
        for pol in lin_pol_order:
            if pol in vis_data:
                pol_list.append(linear_pol_dict[pol])
        self.polarization_array = np.asarray(pol_list)
        self.data_array = np.zeros((self.Nblts, 1, self.Nfreqs, self.Npols), dtype=np.complex_)
        self.flag_array = np.zeros((self.Nblts, 1, self.Nfreqs, self.Npols), dtype=np.bool_)
        
        for pol, vis in vis_data.iteritems():
            pol_i = pol_list.index(linear_pol_dict[pol])
            self.data_array[:, 0, :, pol_i] = vis
            self.flag_array[:, 0, :, pol_i] = vis_weights_data[pol] <= 0
            
        int_times = bl_info['JDATE'][0]
        bin_offset = bl_info['BIN_OFFSET'][0]
        self.time_array = np.zeros(self.Nblts)
        if self.Ntimes == 1:
            self.time_array.fill(int_times)
        else:
            for ii in range(0, self.Ntimes):
                if ii < (self.Ntimes - 1):
                    self.time_array[bin_offset[ii]:bin_offset[ii + 1]] = int_times[ii]
                else:
                    self.time_array[bin_offset[ii]:] = int_times[ii]
        self.ant_1_array = bl_info['TILE_A'][0] - 1
        self.ant_2_array = bl_info['TILE_B'][0] - 1
        self.Nants_data = np.max([len(np.unique(self.ant_1_array)),len(np.unique(self.ant_2_array))])
        self.antenna_names = bl_info['TILE_NAMES'][0].tolist()
        self.Nants_telescope = len(self.antenna_names)
        self.antenna_numbers = np.arange(self.Nants_telescope)
        self.freq_array = np.zeros((1, self.Nfreqs), dtype=np.float_)
        self.freq_array[0, :] = bl_info['FREQ'][0]
        
