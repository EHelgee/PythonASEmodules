"""File contains attempt to generate gbs in cubic perovskites"""

def CutToGB(atoms,a0,xaxis,yaxis,zaxis):
	""" cut to gb *unit*cell, for large system"""
	
	from ase import Atom,Atoms
	import numpy as np
	import GBtools
	
	Natoms = atoms.get_number_of_atoms()

	# periodicities
	Xper = a0*np.sqrt(xaxis[0]**2+xaxis[1]**2+xaxis[2]**2)
	Yper = a0*np.sqrt(yaxis[0]**2+yaxis[1]**2+yaxis[2]**2)
	Zper = a0*np.sqrt(zaxis[0]**2+zaxis[1]**2+zaxis[2]**2)

	# Find atom close to middle
	pos = atoms.get_positions()
	xmid = (np.max(pos[:,0])+np.min(pos[:,0]))/2.0
	ymid = (np.max(pos[:,1])+np.min(pos[:,1]))/2.0
	zmid = (np.max(pos[:,2])+np.min(pos[:,2]))/2.0
	
	mInd = []

	for g in range(Natoms):
		if np.abs(pos[g,0]-xmid)<a0/2 and np.abs(pos[g,1]-ymid)<a0/2 and np.abs(pos[g,2]-zmid)<a0/2:
			mInd.append(g)

	MidAtom = mInd[0]
	MidPos = pos[MidAtom]
	
	Xpos_plus = MidPos[0]+Xper-0.01
	Ypos_plus = MidPos[1]+Yper-0.01
	Zpos_plus = MidPos[2]+Zper-0.01

	GBtools.axcut(atoms,'x',MidPos[0]-0.01,Xpos_plus)
	GBtools.axcut(atoms,'y',MidPos[1]-0.01,Ypos_plus)
	GBtools.axcut(atoms,'z',MidPos[2]-0.01,Zpos_plus)	

	atoms.set_cell([Xper,Yper,Zper])
	atoms.center()
	

def GBmake(oxide,plane,rotdir,xsize,midplane=['Ba','O']):
	""" atoms = GBmake(oxide,plane,rotdir,xsize) makes a symmetric tilt gb with tilt dir rotdir(list) and gb plane plane(list) in oxide, with each grain xsize periods long"""

	from ase import Atom,Atoms
	from ase.visualize import view
	import GBtools
	import OxideFactory
	import numpy as np

	# Find the crystal directions that should coincide with x,y,z directions
	NewX = np.array(plane)
	NewZ = np.array(rotdir)
	NewY = np.cross(NewZ,NewX)	
	
	# See if NewY can be divided by integer
	if 1 not in NewY:
		for i in range(7,1,-1):
			NewY_b = NewY/i*i
			trueval = np.sum(NewY == NewY_b)
			if trueval == len(NewY):
				NewY = NewY/i
				break


	# Find grain boundary periodicity
	PeriodX = np.sqrt(NewX[0]**2 + NewX[1]**2 + NewX[2]**2)
	PeriodY = np.sqrt(NewY[0]**2 + NewY[1]**2 + NewY[2]**2)
        PeriodZ = np.sqrt(NewZ[0]**2 + NewZ[1]**2 + NewZ[2]**2)

	# Make slabs of the right oxide and rotate
	if oxide == 'BaZrO3' or oxide == 'bzo' or oxide == 'BZO':
		a = 4.19134
		slab1 = OxideFactory.bzocube(10)

	elif oxide == 'BaCeO3' or oxide == 'bco' or oxide == 'BCO':
		a = 4.4778
		slab1 = OxideFactory.bcocube(10)
		
		
	GBtools.rotations(slab1,NewX,NewY,NewZ)

	CutToGB(slab1,a,NewX,NewY,NewZ)

	slab1 = slab1.repeat([xsize,1,1])
	
	# The following bit is added because I want a specific set of structures
	slab2 = slab1.repeat([1,1,2])
	GBtools.rotations(slab2,[1,0,0],[0,-1,0],[0,0,1])
	
	# Find the top layer in Z
	pos2 = slab2.get_positions()
	zmax2 = np.max(pos2[:,2])
	
	zlist = []
	for z in pos2:
		if z[2] > zmax2-2 and z[2] < zmax2-0.5:
			zlist.append(z[2])

	layerdiff_z = zmax2-np.max(zlist)
	
	# Cut back to 1 period in z
	lowlim = zmax2-layerdiff_z/2.0-PeriodZ*a
	highlim = zmax2-layerdiff_z/2.0

	GBtools.axcut(slab2,'z',lowlim,highlim)

	# Find middle plane
	pos1 = slab1.get_positions()
	xmax1 = np.max(pos1[:,0])
	ymin1 = np.min(pos1[:,1])
	ymax1 = np.max(pos2[:,1])

	pos2 = slab2.get_positions()
	xmin2 = np.min(pos2[:,0])
	
	symblist = []
	for atm in slab1:	
		if atm.position[0] > xmax1 - 0.5:
			symblist.append(atm.symbol)

	symblist2 = []
	for atm in slab2:
		if atm.position[0] < xmin2 + 0.5:
			symblist2.append(atm.symbol)
	
	xmax_is_mid = 1
	xmin_is_mid = 1
	for m in midplane:	
		if m not in symblist:
			xmax_is_mid = 0
		if m not in symblist2:
			xmin_is_mid = 0

	#print xmax_is_mid,xmin_is_mid
	
	# find all atom positions in the neighbour planes
	if xmax_is_mid:
		ypos_1_plane = 0
		zmax_1_plane = 0
		for p in pos1:
			if p[0] > xmax1 - a/2.0 and p[0] < xmax1 - 0.5:
				if p[2] > zmax_1_plane:
					zmax_1_plane = p[2]
					ypos_1_plane = p[1]

		ypos_2_plane = 0
		zmax_2_plane = 0
		for p in pos2:
			if p[0] < xmin2+0.5:
				if zmax_2_plane < p[2]:
					zmax_2_plane = p[2]
					ypos_2_plane = p[1]

	if xmin_is_mid:
		ypos_1_plane = 0
		zmax_1_plane = 0
		for p in pos1:
			if p[0] > xmax1 - 0.5:
				if p[2] > zmax_1_plane:
					zmax_1_plane = p[2]
					ypos_1_plane = p[1]
		
		ypos_2_plane = 0
		zmax_2_plane = 0
		for p in pos2:
			if p[0] > xmin2 + 0.5 and p[0] < xmin2 + a/2.0:
				if zmax_2_plane < p[2]:
					zmax_2_plane = p[2]
					ypos_2_plane = p[1]
			
	translate_y = ypos_1_plane - ypos_2_plane
	translate_z = zmax_1_plane - zmax_2_plane

	# align top z, top y, x

	plist = []
		
	for p in pos1:
		if p[0] > xmax1-2 and p[0] < xmax1-0.5:
			plist.append(p[0])	

	layerdiff_x = xmax1-np.max(plist)	

	slab2.translate([xmax1+layerdiff_x-xmin2,translate_y,translate_z])
	pos2 = slab2.get_positions()

	ymin2 = np.min(pos2[:,1])
	ymax2 = np.max(pos2[:,1])


	if ymin1 > ymin2:
		slab2 = slab2.repeat([1,3,1])
		GBtools.axcut(slab2,'y',ymin1-0.5,ymax1+0.5)
	else:
		slab1 = slab1.repeat([1,3,1])
		GBtools.axcut(slab1,'y',ymin2-0.5,ymax2+0.5)


	slab1.extend(slab2)

	slab1.set_cell([2*xsize*PeriodX*a,PeriodY*a,PeriodZ*a])
	slab1.center()

	return slab1
