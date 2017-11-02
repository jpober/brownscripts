import pyephem as ephem

day = '1998/8/10 23:10:00'
longitude = ephem.degrees('-1.91667')
latitude = ephem.degrees('52.5')

star = ephem.FixedBody()
star._ra = '16:41:42.0'
star._dec = '36:28:00.0'

observer = ephem.Observer()
observer.date = day
observer.lon = longitude
observer.lat = latitude

star.compute(observer)

print 'Observer', observer
print 'RA', star.ra, 'DEC', star.dec
print 'AZ', star.az, 'ALT', star.alt
