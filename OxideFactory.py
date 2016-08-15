"""File for making oxides"""
def bzocube(x=0,y=0,z=0):
    	"""name=bzocube(x,y,z) makes bzo slab, xa by ya by za"""
    	from ase import Atom,Atoms
    	from math import sqrt

    	a=4.19134 #lattice parameter
    	d=a/2
    	l=(sqrt(2)*a)/2

	#create object
    	name=Atoms('BaZrOOO',[[0,0,0],[d,d,d],[d,d,0],[0,d,d],[d,0,d]])
    	name.set_cell([a,a,a])
	#Set mass and charge later
    	if y and z:
        	name=name.repeat([x,y,z]) #create larger cube
    	else:
        	name=name.repeat(x)
        

    	return name

def bcocube(x=0,y=0,z=0):
	"""name = bcocube(x,y,z) makes cubic bco slab, xa by ya by za"""
	
	from ase import Atoms,Atom
	
	# Lattice parameter from PBE calcs.
	a = 4.4778
	d = a/2
	
	# create unit cell
	name = Atoms('BaCeOOO',[[0,0,0],[d,d,d],[d,d,0],[0,d,d],[d,0,d]])
	name.set_cell([a,a,a])
  	
	# repeat acc. to input
	if y and z:
		name = name.repeat([x,y,z])
	else:
		name = name.repeat(x)

	return name

def BaO(x=1,y=0,z=0):
    """name=BaO(x,y,z) makes BaO slab, xa by ya by za"""
    from ase import Atoms,Atom

    a=5.5  #Lattice parameter
    d=a/2.

    #Create object
    name=Atoms('BaBaBaBaOOOO',[[0,0,0],[d,d,0],[0,d,d],[d,0,d],[d,d,d],[0,0,d],[d,0,0],[0,d,0]])
    name.set_cell([a,a,a])
    if y and z:
        name=name.repeat([x,y,z])
    else:
        name=name.repeat(x)

    return name

def ZrO2(x=1,y=0,z=0):
    """name=ZrO2(x,y,z) makes ZrO2 slab, xa by ya by za, cubic fluorite structure"""
    from ase import Atoms,Atom

    a=5.1328
    d=a/2
    l=a/4

    tile=Atoms('ZrOO',[[0,0,0],[l,l,l],[l,l,-l]])
    name=tile.copy()
    perm=[[d,d,0],[0,d,d],[d,0,d]]
    for p in perm:
        tile2=tile.copy()
        tile2.translate(p)
        name.extend(tile2)
        
    name.set_cell([a,a,a])
    if y and z:
        name=name.repeat([x,y,z])
    else:
        name=name.repeat(x)

    return name
