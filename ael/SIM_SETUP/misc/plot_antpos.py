#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--aspect_neq', action='store_true', help='Do not force equal aspect in x/y plot.')
o.add_option('--nonumbers', action='store_true', help='Do not plot antenna numbers next to symbols.')
opts, args = o.parse_args(sys.argv[1:])

th = n.arange(0, 2*n.pi, .01)
r = 5.

aa = a.cal.get_aa(opts.cal, .1, .1, 1)
antpos = [aa.get_baseline(0,i,src='z') for i in range(len(aa.ants))]
antpos = n.array(antpos) * a.const.len_ns / 100.
x,y,z = antpos[:,0], antpos[:,1], antpos[:,2]
x -= n.average(x)
y -= n.average(y)
p.plot(x,y, 'k.')
if False:
    im = a.img.Img(size=300, res=30)
    DIM = 300./30
    im.put((x,y,z), z)
    _z = a.img.recenter(im.uv, (DIM/2,DIM/2))
    print _z
    _z = n.ma.array(_z, mask=n.where(_z == 0, 1, 0))
    _x,_y = im.get_uv()
    _x = a.img.recenter(_x, (DIM/2,DIM/2))
    _y = a.img.recenter(_y, (DIM/2,DIM/2))
    p.contourf(_x,_y,_z,n.arange(-5,5,.5))
for ant,(xa,ya,za) in enumerate(zip(x,y,z)):
    hx,hy = r*za*n.cos(th)+xa, r*za*n.sin(th)+ya
    if za > 0: fmt = '#eeeeee'
    else: fmt = '#a0a0a0'
    #p.fill(hx,hy, fmt)
    if not opts.nonumbers: p.text(xa,ya, str(ant))
p.grid()
#p.xlim(-100,100)
p.xlabel("East-West Antenna Position (m)")
p.ylabel("North-South Antenna Position (m)")
#p.ylim(-100,100)
a = p.gca()
if not opts.aspect_neq: a.set_aspect('equal')
p.show()
