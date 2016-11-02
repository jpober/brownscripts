import numpy as n

#Print an ant_layout for a given Nants as a hexagon.

Nants = 320

#def hex(n):
#    for x in [(n-abs(x-int(n/2))) for x in range(n)]:
#        for y in range(n-x):
#            print(' '),
#    for y in range(x):
#        print(' * '),
#    print('')
#Code to draw a hexagon.



def hex(m):
    na = 0
    rownum = 0
    ar = []
    for x in [(m-abs(x-n.ceil(m/2))) for x in range(0,m,1)]:
        ar.append([])
        for y in n.arange((m-x)):
            print('.'),
            ar[rownum].append(-1)
        for y in range(int(x)):
              print(' Y '),
              ar[rownum].append(na)
              na += 1
        for y in n.arange((m-x)):
            print('.'),
            ar[rownum].append(-1)
        print('')
        rownum += 1
    return ar

#t = n.floor(n.sqrt(2*Nants))
t = n.floor(2*n.sqrt(2*Nants/(3*n.sqrt(3))))
print hex(int(20))

