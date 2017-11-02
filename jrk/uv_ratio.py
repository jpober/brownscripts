import uvdata
import numpy as n

f1 = uvdata.miriad.Miriad()
f2 = uvdata.miriad.Miriad()

f1.read_miriad('/users/jkerriga/data/jkerriga/128DayOutput/fhd_0/vis_data/Pzen.2456242.30605.uvcRREcACOTUcHP')
f2.read_miriad('/users/jkerriga/zen.2456242.30605.uvcRREcACOTUcP')

f2.flag_array = f1.flag_array

f2.write_miriad('OriginalwithFHDflags')
