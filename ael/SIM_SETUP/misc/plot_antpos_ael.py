#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--aspect_neq', action='store_true', help='Do not force equal aspect in x/y plot.')
o.add_option('--nonumbers', action='store_true', help='Do not plot antenna numbers next to symbols.')
o.add_option('-b', dest='highlight', help='Highlight listed antennas.')
opts, args = o.parse_args(sys.argv[1:])


th = n.arange(0, 2*n.pi, .01)
r = 5.

if not opts.highlight == None:
	mark = map(int,opts.highlight.split(","))
else:
	mark = []


aa = a.cal.get_aa(opts.cal, .1, .1, 1)
antpos = [aa.get_baseline(0,i,src='z') for i in range(len(aa.ants))]
antpos = n.array(antpos) * a.const.len_ns / 100.

x,y,z = antpos[:,0], antpos[:,1], antpos[:,2]
#x -= n.average(x); y -= n.average(y)

other = [j for j in range(len(aa.ants)) if j not in mark]
p.plot(x[mark],y[mark], 'ro',markersize=5)
p.plot(x[other],y[other], 'k.')


for ant,(xa,ya,za) in enumerate(zip(x,y,z)):
    if not opts.nonumbers:
	p.text(xa,ya, str(ant))
    	if ant in mark: p.text(xa,ya, str(ant), color='red')
p.grid()

#print n.sqrt((x[mark[0]] - x[mark[1]])**2 + (y[mark[0]] - y[mark[1]])**2)
#print x[mark[0]] - x[mark[1]]

#p.xlim(-100,100)
p.xlabel("East-West Antenna Position (m)")
p.ylabel("North-South Antenna Position (m)")
#p.ylim(-100,100)
a = p.gca()
if not opts.aspect_neq: a.set_aspect('equal')
p.show()
