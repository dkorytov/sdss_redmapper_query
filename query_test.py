#!/usr/bin/env python
import numpy as np
import sqlcl
from StringIO import StringIO
ra = 239.583329
dec = 0#27.233413
rad = 20.0
radians = (rad/60.0)*np.pi/180.0
solid_angle = 2.0*np.pi*(1.0-np.cos(radians))
solid_angle2 = 2.0*np.pi*(radians**2/2.0)
area = np.pi*rad**2
print "radians: %f solid angle: %e solid angle2: %e square deg: %f "%(radians,solid_angle,solid_angle2,solid_angle2/(np.pi*4)*41253)
print "sq arcmin: %f sq degree: %f" %(np.pi*rad**2, np.pi*(rad/60.0)**2)
result = sqlcl.query("select p.ra, p.dec from PhotoObjAll p,  dbo.fGetNearbyObjEq(%f,%f,%f) as r where p.ObjID = r.ObjID"%(ra,dec,rad)).read()

datagal = np.genfromtxt(StringIO(result),names=True,delimiter=",")
print datagal['dec'].size
print "number of elements per sq arcmin: %f" % (datagal['dec'].size/area)
ra1=datagal['ra'].min()
ra2=datagal['ra'].max()
dec1=datagal['dec'].min()
dec2=datagal['dec'].max()

print ra
print "min ra: %f max ra: %f diff: %f diff arcmin: %f"%(ra1,ra2,(ra2-ra1),(ra2-ra1)*60.0)
print "min dec: %f max dec: %f diff: %f diff arcmin: %f"%(dec1,dec2,(dec2-dec1),(dec2-dec1)*60.0)
exit()
result = sqlcl.query("select p.ra, p.dec from PhotoObjAll p join  dbo.fGetNearbyObjEq(%f,%f,%f) r on p.ObjID = r.ObjID "%(ra,dec,rad)).read()
datagal = np.genfromtxt(StringIO(result),names=True,delimiter=",")
exit()
print datagal['dec'].size
print "number of elements per sq arcmin: %f" % (datagal['dec'].size/area)
print sqlcl.query("select count(ra) from PhotoObjAll").read()
