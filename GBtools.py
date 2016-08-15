"""File contains tools to handle gb production"""

#cutting cuts a bzo cube "name" perp. to 'vector' direction, leaving atoms that are closer than "k" or "m" times the lattice parameter to the plane
def cutting(name, xvector,yvector,zvector,k = 1, m = 1):
    """cutting(name,xvector,yvector,zvector,k,m) cuts a bzo cube name perp. to vector direction, leaving atoms that are closer than k or m times the lattice parameter to the plane"""
    
    from ase import Atoms,Atom
    import pylab as plt
    from math import sqrt
    g=name.get_center_of_mass() #moving COM to origo
    name.translate(-g)
    latticeparam = 4.19134
    l1 = latticeparam * k
    l2 = latticeparam * m
    absolute = sqrt(zvector[0]**2 + zvector[1]**2 + zvector[2]**2)
    normx = zvector[0] / absolute
    normy = zvector[1] / absolute
    normz = zvector[2] / absolute
    n = name.get_number_of_atoms()
    pos = name.get_positions()
    tagged = 0
    for i in range(0, n, 1):
        x = pos[i][0]
        y = pos[i][1]
        z = pos[i][2]
        dotprod = x * normx + y * normy + z * normz
        if abs(dotprod)>l2 and dotprod > 0:
            name[i].set_tag(1)
            tagged = tagged + 1
        elif abs(dotprod)>l1 and dotprod < 0:
            name[i].set_tag(1)
            tagged = tagged + 1
        else:
            name[i].set_tag(0)

    removed=0
    j=0
    k=0
    while removed<tagged and k<n:
        k=k+1
        s=name[j].get_tag()
        if s==1:
            name.pop(j)
            removed=removed+1
        else:
            j=j+1

    name.center()

        #Getting cartesian coordinates for vectors
    cell = name.cell
    vector0 = plt.dot(cell,xvector)
    vector1 = plt.dot(cell,yvector)
    vector2 = plt.dot(cell,zvector)

    #Introduce help atoms (yes this code is ugly)
    zer = Atom('Al',(0, 0, 0),tag = 10)
    xatom = Atom('Al',vector0,tag = 11)
    yatom = Atom('Al',vector1,tag = 12)
    zatom = Atom('Al',vector2,tag = 13)

    name.append(zer)
    name.append(xatom)
    name.append(yatom)
    name.append(zatom)

    #first rotation:xvector to x
    name.rotate(vector0,'x')

    num = name.get_number_of_atoms()

    #get new yvector
    for i in range(0, num):
        tag = name[i].get_tag()
        if tag == 12:
            ypos = name[i].get_position()

    #determine angle and rotate
    absolute=sqrt(ypos[0]**2+ypos[1]**2+ypos[2]**2)
    sine=ypos[1]/absolute
    angle=plt.arccos(sine)
    name.rotate('x',angle)

    removed = 0
    j = 0
    while removed < 4 and j < num:
        tag = name[j].get_tag()
        if tag > 9:
            name.pop(j)
            removed = removed + 1
        else:
            j = j + 1

    name.center()

    
#This function rotates a slab of material so that "xvector" coincides with
#the ASE's x direction, the y vector coincides with the ASE's y direction
#and so on. The vectors are assumed to be scaled with the cell
# rather than cartesian. 
def rotations(name,xvector,yvector,zvector):
    """rotations(name,xvector,yvector,zvector) rotates a slab of material so that xvector conincides with the ASE x direction and corresp. fro y and z. The vectors are assumed to be scaled with the cell rather than cartesian"""
    
    from ase import Atoms,Atom
    import pylab as plt
    from math import sqrt
    #Getting cartesian coordinates for vectors
    cell = name.cell
    vector0 = plt.dot(cell,xvector)
    vector1 = plt.dot(cell,yvector)
    vector2 = plt.dot(cell,zvector)

    #Introduce help atoms (yes this code is ugly)
    zer = Atom('Al',(0, 0, 0),tag = 10)
    xatom = Atom('Al',vector0,tag = 11)
    yatom = Atom('Al',vector1,tag = 12)
    zatom = Atom('Al',vector2,tag = 13)

    name.append(zer)
    name.append(xatom)
    name.append(yatom)
    name.append(zatom)

    #first rotation:xvector to x
    name.rotate(vector0,'x')

    g = name.get_number_of_atoms()

    #get new yvector
    for i in range(0, g):
        tag = name[i].tag
        if tag == 12:
            ypos = name[i].position

    #determine angle and rotate
    absolute=sqrt(ypos[0]**2+ypos[1]**2+ypos[2]**2)
    sine=ypos[1]/absolute
    angle=plt.arccos(sine)
    name.rotate('x',angle)

    removed = 0
    j = 0
    while removed < 4 and j < g:
        tag = name[j].tag
        if tag > 9:
            name.pop(j)
            removed = removed + 1
        else:
            j = j + 1
        

    name.center()    
    


def axcut(name,axis,min=0,max=0):
    """axcut(name,axis,min,max) cuts atoms from object name above and below specified max and min on axis axis"""
    
    from ase import Atoms,Atom
    import pylab as plt
    from math import sqrt
    n=name.get_number_of_atoms()
    if axis == 'x':
        coordinate=0
    elif axis == 'y':
        coordinate=1
    else:
        coordinate=2
        
    pos=name.get_positions()
    tagged=0
    ra=range(0,n,1)
    for i in ra:
        if min:
            if pos[i][coordinate]<=min:
                name[i].tag=1
                tagged=tagged+1
        

        if max:
            if pos[i][coordinate]>=max:
                name[i].tag=1
                tagged=tagged+1


    removed=0
    j=0
    k=0
    while removed<tagged and k<n:
        k=k+1
        s=name[j].tag
        if s==1:
            name.pop(j)
            removed=removed+1
        else:
            j=j+1
