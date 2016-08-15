"""Functions to make graphene sheets with certain size and C-C distance (or strain in 1 dir). The sheets will be periodically repeatable"""

def rescale_ribbon(xmax,zmax,scaling):
    """rescale_ribbon(xmax,zmax,scaling) creates a graphene nanoribbon of size xmax times zmax or smaller, where the lattice is scaled by the parameter scaling"""

    from ase import Atoms,Atom
    from ase.structure import graphene_nanoribbon
    import lammpsio
    import math
    from numpy import dot

    #Unit cell param in x 2.13
    #Unit cell param in z 2.4595

    ax=2.13 # Unit cell param in x acc. to graphene_nanoribbon
    az=2.4595 # Unit cell param in z acc to graphene_nanoribbon
    X=math.floor(xmax/ax)
    Z=int(math.floor(zmax/az))

    #No of cells in x must be even for repeatability
    xhalf=X/2.0
    if xhalf!=math.floor(xhalf):
        X=X+1

    X=int(X)
    
    #Remove vacuum
    sheet=graphene_nanoribbon(X,Z)
    cell=sheet.get_cell()
    cell[0][0]=cell[0][0]-5
    sheet.set_cell(cell)

    #scale
    scalepos=sheet.get_scaled_positions()
    cell[0]=cell[0]*scaling
    cell[2]=cell[2]*scaling
    positions=dot(scalepos,cell)
    sheet.set_positions(positions)
    sheet.set_cell(cell)
    sheet.center()
    cell[1]=cell[1]
    sheet.set_cell(cell)

    return sheet
    
def strained_ribbon(xmax,zmax,scaling,axis):
	"""strained_ribbon(xmax,zmax,scaling,axis) creates a graphene nanoribbon of size xmax times zmax or smaller, where the distances in direction axis are scaled by parameter scaling RELATIVE TO THE C-C DIST 1.438785 GIVEN BY MODIFIED TERSOFF!!!"""

	from ase import Atoms,Atom
    	from ase.structure import graphene_nanoribbon
    	import lammpsio
    	import math
    	from numpy import dot
	
	sheet=rescale_ribbon(xmax,zmax,1.013)

    	#rescale
    	scalepos=sheet.get_scaled_positions()
	cell=sheet.get_cell()
	if axis=='x':    	
		cell[0]=cell[0]*scaling
	elif axis=='z':
    		cell[2]=cell[2]*scaling
	
    	positions=dot(scalepos,cell)
    	sheet.set_positions(positions)
    	sheet.set_cell(cell)
    	sheet.center()
    	cell[1]=cell[1]*6
    	sheet.set_cell(cell)

    	return sheet

def rescale_ribbon_arm(xmax,zmax,scaling):
    """rescale_ribbon(xmax,zmax,scaling) creates a graphene nanoribbon of size xmax times zmax or smaller, where the lattice is scaled by the parameter scaling"""

    from ase import Atoms,Atom
    from ase.structure import graphene_nanoribbon
    import lammpsio
    import math
    from numpy import dot

    #Unit cell param in x 2.13
    #Unit cell param in z 2.4595

    ax=2.4595 # Unit cell param in x acc. to graphene_nanoribbon
    az=4.26 # Unit cell param in z acc to graphene_nanoribbon
    X=math.floor(xmax/ax)
    Z=int(math.floor(zmax/az))

    #No of cells in x must be even for repeatability
    xhalf=X/2.0
    if xhalf!=math.floor(xhalf):
        X=X+1

    X=int(X)
    
    #Remove vacuum
    sheet=graphene_nanoribbon(X,Z,type='armchair')
    cell=sheet.get_cell()
    cell[0][0]=cell[0][0]-5
    sheet.set_cell(cell)

    #scale
    scalepos=sheet.get_scaled_positions()
    cell[0]=cell[0]*scaling
    cell[2]=cell[2]*scaling
    positions=dot(scalepos,cell)
    sheet.set_positions(positions)
    sheet.set_cell(cell)
    sheet.center()
    cell[1]=cell[1]
    sheet.set_cell(cell)

    return sheet
    
def strained_ribbon_arm(xmax,zmax,scaling,axis):
	"""strained_ribbon(xmax,zmax,scaling,axis) creates a graphene nanoribbon of size xmax times zmax or smaller, where the distances in direction axis are scaled by parameter scaling RELATIVE TO THE C-C DIST 1.438785 GIVEN BY MODIFIED TERSOFF!!!"""

	from ase import Atoms,Atom
    	from ase.structure import graphene_nanoribbon
    	import lammpsio
    	import math
    	from numpy import dot
	
	sheet=rescale_ribbon_arm(xmax,zmax,1.013)

    	#rescale
    	scalepos=sheet.get_scaled_positions()
	cell=sheet.get_cell()
	if axis=='x':    	
		cell[0]=cell[0]*scaling
	elif axis=='z':
    		cell[2]=cell[2]*scaling
	
    	positions=dot(scalepos,cell)
    	sheet.set_positions(positions)
    	sheet.set_cell(cell)
    	sheet.center()
    	cell[1]=cell[1]*6
    	sheet.set_cell(cell)

    	return sheet
