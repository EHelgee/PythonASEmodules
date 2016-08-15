"""Try to find a good general way of constructing graphene grain boundaires"""

def graphenegb(Boundaryangle,xlength,zlength):
	"""graphenegb(angle=rotation angle, xlength=aimed for length,zlength=aimed for width) creates a grain boundary in (armchair oriented) graphene"""
	
	from ase import Atom,Atoms,visualize
	from bzo import axcut
	from sheetmake import strained_ribbon
	from numpy import transpose,argmin,argmax
	from math import pi
	
	#Grain characteristics
	xlength_grain=xlength*1.25/2.0
	zlength_grain=zlength*1.25
	rotangle=Boundaryangle*pi/2.0/180.0
	periodicity=GB_period(rotangle) # armchair 28.7 boundary: 8.983. 15.18 boundary: 16.338, 21.4: 11.417

	#Make grains
	grain1=strained_ribbon(xlength_grain,zlength_grain,1,'x')
	grain2=strained_ribbon(xlength_grain,zlength_grain,1,'x')

	#Rotate
	grain1.rotate('y',rotangle)
	grain2.rotate('y',-rotangle)

	#Cut perp. to new x direction
	#Find corner that is min. in z direction
	pos1=grain1.get_positions()
	pos1_transpose=pos1.transpose()
	z1min=argmin(pos1_transpose[2])
	z1max=argmax(pos1_transpose[2])
	minpos1=pos1[z1min]
	maxpos1=pos1[z1max]

	pos2=grain2.get_positions()
	pos2_transpose=pos2.transpose()
	z2min=argmin(pos2_transpose[2])
	z2max=argmax(pos2_transpose[2])
	minpos2=pos2[z2min]
	maxpos2=pos2[z2max]

	#Make x direction cuts
	axcut(grain1,'x',maxpos1[0],minpos1[0])
	axcut(grain2,'x',minpos2[0]-0.1,0)

	#Displace grain 2
	xdiff=minpos1[0]-minpos2[0]+0.05
	zdiff=minpos1[2]-minpos2[2]+periodicity/2

	grain2.translate([xdiff,0,zdiff])
	
	#Join grains, adjust cell
	grain1.extend(grain2)
	
	pos1=grain1.get_positions()
	pos1_transpose=pos1.transpose()
	xminind=argmin(pos1_transpose[0])
	xmaxind=argmax(pos1_transpose[0])
	zminind=argmin(pos1_transpose[2])
	zmaxind=argmax(pos1_transpose[2])
	
	xtop=pos1[xmaxind]
	xbot=pos1[xminind]
	ztop=pos1[zmaxind]
	zbot=pos1[zminind]

	axcut(grain1,'z',xtop[2],0)
	axcut(grain1,'x',0,ztop[0])

	pos1=grain1.get_positions()
	pos1_transpose=pos1.transpose()
	xminind=argmin(pos1_transpose[0])
	xmaxind=argmax(pos1_transpose[0])
	zminind=argmin(pos1_transpose[2])
	zmaxind=argmax(pos1_transpose[2])
	
	xtop=pos1[xmaxind]
	xbot=pos1[xminind]
	ztop=pos1[zmaxind]
	zbot=pos1[zminind]

	grain1.translate([-xbot[0],0,-zbot[2]])
	axcut(grain1,'z',0,periodicity-0.001)
	grain1.set_cell([xtop[0]-xbot[0],50,periodicity])
	grain1.center()
	
	visualize.view(grain1)
	return grain1

def zigzaggb(Boundaryangle,xlength,zlength):
	"""Grain boundaries from zigzag orientation"""

	from ase import Atom,Atoms,visualize
	from bzo import axcut
	from sheetmake import strained_ribbon_arm
	from numpy import transpose, argmax,argmin
	from math import pi
	
	#Grain characteristics
	xlength_grain=xlength*1.25/2.0
	zlength_grain=zlength*1.25
	rotangle=Boundaryangle*pi/2.0/180.0
	periodicity=GB_period(rotangle)

	#Make grains
	grain1=strained_ribbon_arm(xlength_grain,zlength_grain,1,'x')
	grain2=strained_ribbon_arm(xlength_grain,zlength_grain,1,'x')

	#Rotate
	grain1.rotate('y',rotangle)
	grain2.rotate('y',-rotangle)

	#Cut perp. to new x direction
	#Find corner that is min. in z direction
	pos1=grain1.get_positions()
	pos1_transpose=pos1.transpose()
	z1min=argmin(pos1_transpose[2])
	z1max=argmax(pos1_transpose[2])
	minpos1=pos1[z1min]
	maxpos1=pos1[z1max]

	pos2=grain2.get_positions()
	pos2_transpose=pos2.transpose()
	z2min=argmin(pos2_transpose[2])
	z2max=argmax(pos2_transpose[2])
	minpos2=pos2[z2min]
	maxpos2=pos2[z2max]

	#Make x direction cuts
	axcut(grain1,'x',maxpos1[0],minpos1[0])
	axcut(grain2,'x',minpos2[0]-0.1,0)

	#Displace grain 2
	xdiff=minpos1[0]-minpos2[0]
	zdiff=minpos1[2]-minpos2[2]

	grain2.translate([xdiff,0,zdiff])
	
	#Join grains, adjust cell
	grain1.extend(grain2)
	
	pos1=grain1.get_positions()
	pos1_transpose=pos1.transpose()
	xminind=argmin(pos1_transpose[0])
	xmaxind=argmax(pos1_transpose[0])
	zminind=argmin(pos1_transpose[2])
	zmaxind=argmax(pos1_transpose[2])
	
	xtop=pos1[xmaxind]
	xbot=pos1[xminind]
	ztop=pos1[zmaxind]
	zbot=pos1[zminind]

	axcut(grain1,'z',xtop[2],0)
	axcut(grain1,'x',0,ztop[0])

	pos1=grain1.get_positions()
	pos1_transpose=pos1.transpose()
	xminind=argmin(pos1_transpose[0])
	xmaxind=argmax(pos1_transpose[0])
	zminind=argmin(pos1_transpose[2])
	zmaxind=argmax(pos1_transpose[2])
	
	xtop=pos1[xmaxind]
	xbot=pos1[xminind]
	ztop=pos1[zmaxind]
	zbot=pos1[zminind]

	grain1.translate([-xbot[0],0,-zbot[2]])
	axcut(grain1,'z',0,periodicity-0.004)
	grain1.set_cell([xtop[0]-xbot[0],50,periodicity])
	grain1.center()
	
	visualize.view(grain1)
	return grain1	

def additional(sheet,x,y,zlist=[]):
	""" Add atoms to the boundaries"""

	from ase import Atoms,Atom
	
	natoms=len(zlist)
	symbols=str(natoms)+'C'
	
	if natoms==0:
		raise ValueError("Not adding any atoms!")

	poslist=[]
	for i in range(natoms):
		poslist.append([x,y,zlist[i]])

	addatoms=Atoms(symbols,poslist)
	
	sheet.extend(addatoms)

def makezlist(z0,N,start):
	""" makezlist(z0,N,start) Making the zlist for additional. if start==1, begin with long jump"""

	zl=[]
	zl.append(z0)
	zc=z0
	if start==1:
		a=1.44
		b=2.88
	else:
		a=2.88
		b=1.44
	
	for i in range(N):
		if (i/2)*2==i:
			z=zc+b
			zl.append(z)
		else:
			z=zc+a
			zl.append(z)

		zc=z

	return zl

def GB_period(angle):
	"""Find the GB period from ROTATION (not misorientation) angle in radians"""

	from math import atan,sqrt,pi
	import numpy as np

	angle=pi/6-angle
	a0=sqrt(3)*1.438
	N=50
	difflist=[]
	mlist=[]
	nlist=[]
	Llist=[]

	for n in range(N):
		for m in range(1,N):
			theta=atan(sqrt(3)*m/(2*n+m))
			diff=abs(theta-angle)
			if diff<=0.001:
				L=a0*sqrt(n**2+n*m+m**2)
				difflist.append(diff)
				mlist.append(m)
				nlist.append(n)
				Llist.append(L)
			

	indarray=np.array(mlist)+np.array(nlist)
	
	firstind=np.argmin(indarray)
	print 'm='+str(mlist[firstind])+' , n='+str(nlist[firstind])+' , L='+str(Llist[firstind])+' , diff='+str(difflist[firstind])

	return Llist[firstind]

def lengthen(inatoms,angle,gbpos,endlength):
	"""lengthen(atoms,angle,gbpos,endlength) makes a longer cell from gb cell atoms with misorientation angle angle and gb at gbpos. endlength is the aimed-for length. Assumes y perp. to plane and x perp to gb."""

	from ase import Atoms,Atom
	from ase.visualize import view
	from math import pi,floor
	from bzo import axcut
	
	rotangle=angle*pi/180/2-pi/6
	Xperiod=GB_period(rotangle)
	
	atoms=inatoms.copy()
	
	c=atoms.get_cell()
	
	addition=int(floor((endlength/Xperiod-1)/2))	#No of cells to add on each side

	rightside=atoms.copy()
	leftside=atoms.copy()
	rightedge=atoms.copy()
	leftedge=atoms.copy()
	
	#Taking on the right side
	axcut(rightside,'x',gbpos+Xperiod,gbpos+2*Xperiod)
	rc=rightside.get_cell()
	rightside.set_cell([Xperiod,rc[1][1],rc[2][2]])
	rightside.center()
	rightside=rightside.repeat([addition,1,1])
	rightside.translate([gbpos+Xperiod,0,0])

	#Left side
	axcut(leftside,'x',gbpos-2*Xperiod,gbpos-Xperiod)
	lc=leftside.get_cell()
	leftside.set_cell([Xperiod,lc[1][1],lc[2][2]])
	leftside.center()	
	leftside=leftside.repeat([addition,1,1])
	leftside.translate([-addition*Xperiod+(gbpos-Xperiod),0,0])


	axcut(leftedge,'x',0,gbpos-Xperiod)
	leftedge.set_cell([(gbpos-Xperiod),lc[1][1],lc[2][2]])
	leftedge.translate([(-addition*Xperiod),0,0])
	axcut(rightedge,'x',gbpos+Xperiod,c[0][0]*2)
	rightedge.set_cell([(rc[0][0]-gbpos-Xperiod),rc[1][1],rc[2][2]])
	rightedge.center()
	rightedge.translate([(gbpos+Xperiod+addition*Xperiod-0.3),0,0])

	leftside.extend(leftedge)
	rightside.extend(rightedge)

	axcut(atoms,'x',gbpos-Xperiod,gbpos+Xperiod)

	atoms.extend(rightside)
	atoms.extend(leftside)

	newxlength=2*Xperiod*(addition+1)+(gbpos-Xperiod)+rc[0][0]-gbpos-Xperiod
	atoms.set_cell([newxlength,c[1][1],c[2][2]])
	atoms.center()

	return atoms
		
def profile(atoms):
	""" plot the profile along x"""

	from ase import Atoms,Atom
	import pylab as p
	import numpy as np
	from math import ceil,floor
	import string

	pos=atoms.get_positions()
	N=atoms.get_number_of_atoms()
	
	xinterval_desired=1.4
	xlength=np.max(pos[:,0])-np.min(pos[:,0])
	xbins=floor(xlength/xinterval_desired)
	xinterval=xlength/xbins

	yvec=np.zeros(xbins)
	ycount=np.zeros(xbins)

	ymax=np.max(pos[:,1])
	ymin=np.min(pos[:,1])
	
	xvec=(np.array(range(int(xbins))))*xinterval

	for i in range(N):
		bin=floor((pos[i,0]-np.min(pos[:,0]))/xinterval)	
		yvec[bin-1]=yvec[bin-1]+pos[i,1]
		ycount[bin-1]=ycount[bin-1]+1

	yvec=yvec/ycount

	lowmean=floor(len(yvec)/4)-floor(len(yvec)/8)
	highmean=floor(len(yvec)/4)+floor(len(yvec)/8)
	ymean=np.mean(yvec[lowmean:highmean])
	if ymax-ymean > ymean-ymin:
		yheight=ymax-ymean
	else:
		yheight=ymean-ymin

	print xvec[lowmean],xvec[highmean],yheight

	strings=[]
	for j in range(len(yvec)):
		line=str(xvec[j])+' , '+str(yvec[j])
		strings.append(line)

	text=string.join(strings,'\n')
	f=open('profile.dat','w')
	f.write(text)
	f.close()

	p.figure()
	p.plot(xvec,yvec)
	p.xlabel('X (AA)')
	p.ylabel('Y (AA)')
	p.show()

def MAXprofile(atoms):
	""" plot the profile / height in y coord. as function of x coord."""

	from ase import Atoms,Atom
	import pylab as p
	import numpy as np
	from math import ceil,floor,sqrt
	import string

	pos=atoms.get_positions()
	N=atoms.get_number_of_atoms()
	
	xinterval_desired=2
	xlength=np.max(pos[:,0])-np.min(pos[:,0])
	xbins=floor(xlength/xinterval_desired)
	xinterval=xlength/xbins

	yvec=np.zeros(xbins)
	ycount=np.zeros(xbins)

	ymax=np.max(pos[:,1])
	ymin=np.min(pos[:,1])
	
	xvec=(np.array(range(int(xbins))))*xinterval

	for i in range(N):
		bin=floor((pos[i,0]-np.min(pos[:,0]))/xinterval)
		if yvec[bin-1]<	sqrt(pos[i,1]**2):
			yvec[bin-1]=pos[i,1]
		
	lowmean=floor(len(yvec)/4)-floor(len(yvec)/8)
	highmean=floor(len(yvec)/4)+floor(len(yvec)/8)
	ymean=np.mean(yvec[lowmean:highmean])
	if ymax-ymean > ymean-ymin:
		yheight=ymax-ymean
	else:
		yheight=ymean-ymin

	print xvec[lowmean],xvec[highmean],yheight

	strings=[]
	for j in range(len(yvec)):
		line=str(xvec[j])+' , '+str(yvec[j])
		strings.append(line)

	text=string.join(strings,'\n')
	f=open('profile.dat','w')
	f.write(text)
	f.close()

	p.figure()
	p.plot(xvec,yvec)
	p.xlabel('X (AA)')
	p.ylabel('Y (AA)')
	p.show()

def y_profile(atoms):
	""" The profile / height, counting max. positions, in y """

	from ase import Atoms,Atom
	import pylab as p
	import numpy as np
	from math import ceil,floor,sqrt
	import string

	pos=atoms.get_positions()
	N=atoms.get_number_of_atoms()

	zpos=pos[:,1]-np.mean(pos[:,1])
	peakind = np.where(zpos == max(zpos))[0][0]
	xmiddle = np.max(pos[:,0])/2.0
	
	yinterval_desired=1.5
	ylength=np.max(pos[:,2])-np.min(pos[:,2])
	ybins=floor(ylength/yinterval_desired)
	yinterval=ylength/ybins

	zvec=np.zeros(ybins)
	
	yvec=(np.array(range(int(ybins))))*yinterval

	for i in range(N):
		bin=floor((pos[i,2]-np.min(pos[:,2]))/yinterval)
		if pos[i,0] > xmiddle -200 and pos[i,0]<xmiddle+200:
			if sqrt(zvec[bin-1]**2) <	sqrt(zpos[i]**2):
				zvec[bin-1]=zpos[i]

	strings=[]
	for j in range(len(zvec)):
		line=str(yvec[j])+' , '+str(zvec[j])
		strings.append(line)

	print max(zvec)-min(zvec)

	text=string.join(strings,'\n')
	f=open('profile.dat','w')
	f.write(text)
	f.close()

	p.figure()
	p.plot(yvec,zvec)
	p.xlabel('Y (AA)')
	p.ylabel('Z (AA)')
	p.show()
	

def cf_profile(atoms1,atoms2):
	""" Calculate and plot change in profile along x"""
	from ase import Atoms,Atom
	import pylab as p
	import numpy as np
	from math import ceil,floor,sqrt
	import string

	pos1=atoms1.get_positions()
	pos2=atoms2.get_positions()
	N=atoms1.get_number_of_atoms()

	ydiff=pos2[:,1]-pos1[:,1]	

	xinterval_desired=1.4
	xlength=np.max(pos2[:,0])-np.min(pos2[:,0])
	xbins=floor(xlength/xinterval_desired)
	xinterval=xlength/xbins

	yvec=np.zeros(xbins)
	ycount=np.zeros(xbins)

	ymax=np.max(ydiff)
	ymin=np.min(ydiff)
	
	xvec=(np.array(range(int(xbins))))*xinterval

	for i in range(N):
		bin=floor((pos2[i,0]-np.min(pos2[:,0]))/xinterval)
		if yvec[bin-1]<	sqrt(ydiff[i]**2):
			yvec[bin-1]=ydiff[i]
		
	lowmean=floor(len(yvec)/4)-floor(len(yvec)/8)
	highmean=floor(len(yvec)/4)+floor(len(yvec)/8)
	ymean=np.mean(yvec[lowmean:highmean])
	if ymax-ymean > ymean-ymin:
		yheight=ymax-ymean
	else:
		yheight=ymean-ymin

	print xvec[lowmean],xvec[highmean],yheight

	strings=[]
	for j in range(len(yvec)):
		line=str(xvec[j])+' , '+str(yvec[j])
		strings.append(line)

	text=string.join(strings,'\n')
	f=open('profile.dat','w')
	f.write(text)
	f.close()

	p.figure()
	p.plot(xvec,yvec)
	p.xlabel('X')
	p.ylabel('Y')
	p.show()
