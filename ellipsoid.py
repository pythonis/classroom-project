import os, sys, re, math
import sys
from math import *
from numpy import *
from pymol import cmd
from pymol.cgo import BEGIN, COLOR, TRIANGLES, VERTEX, NORMAL, END
##########################################################################
#
#
# Example: 2n8z, 5a2h
##########################################################################

def ellipsoid(sel='all'):
    
    model=cmd.get_model(sel+' and name ca')
    obj = cmd.get_object_list()[0]
    coord=[]
    for at in model.atom:
        coord.append(at.coord)
    cm=get_cm(coord)
       
    abc, g, v = get_moi(coord,cm)  # axis-length, eigenvalues, eigenvectors
    
    arrow_obj = []
    color = 'red red'
    pos2 = [x+y*abc[0] for x, y in zip(cm,v[0])]
    arrow_obj += cgo_arrow(cm,pos2,color=color,radius=0.15)

    color = 'green green'
    pos2 = [x+y*abc[1] for x, y in zip(cm,v[1])]
    arrow_obj += cgo_arrow(cm,pos2,color=color,radius=0.15)

    color = 'blue blue'
    pos2 = [x+y*abc[2] for x, y in zip(cm,v[2])]
    arrow_obj += cgo_arrow(cm,pos2,color=color,radius=0.15)

    name = cmd.get_unused_name('axis')
    cmd.load_cgo(arrow_obj, name)
    
    cgo = makeEllipsoid(cm[0], cm[1], cm[2], abc[0], abc[1], abc[2], m = v)
    
    
    obj_ellipsoid = cmd.get_unused_name('ellipsoid')
    cmd.load_cgo(cgo, obj_ellipsoid)

    # Do some cosmetic work
    cmd.set('cgo_transparency', 0.4, obj_ellipsoid)
    cmd.color('green',obj_ellipsoid)
    cmd.show_as('cartoon',obj)
    cmd.color('orange',obj)
    cmd.zoom(obj,animate=1)

    print('The volume of the ellipsoid (4/3 Pi*a*b*c):  %6.1f A' % \
              (4.0/3.0*math.pi*abc[0]*abc[1]*abc[2]))
    print('The volume of the spheroid (4/3 Pi*r**3):  %6.1f A' % \
              (4.0/3.0*math.pi*((abc[0]**2+abc[1]**2+abc[2]**2))**1.5))


cmd.extend('ellipsoid',ellipsoid)

def get_moi(atom,center):
    moment=array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
    abc=[0.0,0.0,0.0]
    n = len(atom) 
    for x,y,z in atom:
        x -= center[0]
        y -= center[1]
        z -= center[2]
        x2 = x**2; y2 = y**2; z2 = z**2
	moment[0][0] += y2+z2
        moment[1][1] += x2+z2
        moment[2][2] += x2+y2
        moment[1][0] -= x*y
        moment[2][0] -= z*x
        moment[2][1] -= y*z
        moment[0][1]=moment[1][0]
        moment[0][2]=moment[2][0]
        moment[1][2]=moment[2][1]

    rg=(moment[0][0]+moment[1][1]+moment[2][2])/n
    rg=sqrt(rg)
    print('Radius of Gyration %6.2f' % rg)

    v=linalg.eigh(moment) # get eigenvalues and eigenvectors

    abc = [sqrt(x/n) for x in v[0]]
    print('principal moi %6.2f %6.2f %6.2f' % (abc[0],abc[1],abc[2]))

    abc = [1.0/x for x in abc] # reciprocal moment of inertia
    r = sum([x*x for x in abc])**0.5 # r = sqrt(a*a+b*b+c*c)
    abc = [rg*x/r for x in abc] # normalize it to the length of rg
    return abc,v[0],v[1].transpose()


def get_cm(atom):
    center=[0,0,0]
    n = len(atom)
    for x,y,z in atom:
        center[0]+=x
        center[1]+=y
        center[2]+=z

    center=[x/n for x in center]

    print('Center of mass: %6.2f %6.2f %6.2f' % \
              (center[0],center[1],center[2]))
    return center

def cal_radius_of_gyration(particles,center):
    r=0.0
    n=0.0
    for x,y,z in particles:
        n=pow(x-center[0],2)+pow(y-center[1],2)+pow(z-center[2],2)
	r=r+n
    r=r/len(particles)
    r=sqrt(r)
    return r
"""
def cal_tensor(particles):
    tensor=array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
    Txx=0.0
    Tyy=0.0
    Tzz=0.0
    Txy=0.0
    Txz=0.0
    Tyz=0.0
    for xi,yi,zi in particles:
        for xj,yj,zj in particles:
            Txx+=pow(xi-xj,2)
            Tyy+=pow(yi-yj,2)
            Tzz+=pow(zi-zj,2)
            Txy+=(xi-xj)*(yi-yj)
            Txz+=(xi-xj)*(zi-zj)
            Tyz+=(yi-yj)*(zi-zj)
    #print Txx
    tensor[0][0]=Txx
    tensor[1][1]=Tyy
    tensor[2][2]=Tzz
    tensor[0][1]=Txy
    tensor[0][2]=Txz
    tensor[1][0]=Txy
    tensor[1][2]=Tyz
    tensor[2][0]=Txz
    tensor[2][1]=Tyz
    tensor=[elem/2/len(particles)/len(particles) for elem in tensor]
    print tensor
    #print ((Txx+Tyy+Tzz)/(2*len(particles)*len(particles)))
    #print (tensor[0][0]+tensor[1][1]+tensor[2][2])/(2*len(particles)*len(particles))
    #print tensor[0][0]+tensor[1][1]+tensor[2][2]
    #print sqrt(tensor[0][0]+tensor[1][1]+tensor[2][2])
    w=linalg.eigh(tensor)
    #print w[0][0]+w[0][1]+w[0][2]
    print sqrt(w[0][0]),sqrt(w[0][1]),sqrt(w[0][2])
    pa=[sqrt(w[0][0]),sqrt(w[0][1]),sqrt(w[0][2])]
    return pa
"""


##########################################################################
#
#  Graphics part: making 3D objects
#
#########################################################################
def signOfFloat(f):
        if f < 0: return -1
        if f > 0: return 1
        return 0

def sqC(v, n):
        return signOfFloat(math.cos(v)) * math.pow(math.fabs(math.cos(v)), n)

def sqS(v, n):
        return signOfFloat(math.sin(v)) * math.pow(math.fabs(math.sin(v)), n)

def sqEllipsoid(x, y, z, a1, a2, a3, u, v, n, e, m = array([[1,0,0],[0,1,0],[0,0,1]])):
    
    x0 = x; y0 = y; z0 = z #center

    x = a1 * sqC(u, n) * sqC(v, e)
    y = a2 * sqC(u, n) * sqS(v, e)
    z = a3 * sqS(u, n)

    nx = sqC(u, 2 - n) * sqC(v, 2 - e) / a1
    ny = sqC(u, 2 - n) * sqS(v, 2 - e) / a2
    nz = sqS(u, 2 - n) / a3
    
    '''
    x1 = m[0][0]*x+m[0][1]*y+m[0][2]*z
    y1 = m[1][0]*x+m[1][1]*y+m[1][2]*z
    z1 = m[2][0]*x+m[2][1]*y+m[2][2]*z
    nx1 = m[0][0]*nx+m[0][1]*ny+m[0][2]*nz
    ny1 = m[1][0]*nx+m[1][1]*ny+m[1][2]*nz
    nz1 = m[2][0]*nx+m[2][1]*ny+m[2][2]*nz
    '''

    # X_old = m*x_new (i.e., its principal axes are already aligned with XYZ)
    # Note the matrix is already transpose, so we need to transpose it back to the original matrix
    #
    x1 = m[0][0]*x+m[1][0]*y+m[2][0]*z
    y1 = m[0][1]*x+m[1][1]*y+m[2][1]*z
    z1 = m[0][2]*x+m[1][2]*y+m[2][2]*z
    nx1 = m[0][0]*nx+m[1][0]*ny+m[2][0]*nz
    ny1 = m[0][1]*nx+m[1][1]*ny+m[2][1]*nz
    nz1 = m[0][2]*nx+m[1][2]*ny+m[2][2]*nz
    
    x1 += x0; y1 += y0; z1 += z0

    return x1, y1, z1, nx1, ny1, nz1

def makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, \
                                  n, e, u1, u2, v1, v2, \
                                  u_segs, v_segs, \
                                  color=[1, 1, 1], \
                                  m = array([[1,0,0],[0,1,0],[0,0,1]])):

        r, g, b = color

        # Calculate delta variables */
        dU = (u2 - u1) / u_segs
        dV = (v2 - v1) / v_segs

        o = [ BEGIN, TRIANGLES, ALPHA, 0.6]

        U = u1
        for Y in range(0, u_segs):
                # Initialize variables for loop */
                V = v1
                for X in range(0, v_segs):
                        # VERTEX #1 */
                        x1, y1, z1, n1x, n1y, n1z = sqEllipsoid(x, y, z, a1, a2, a3, U, V, n, e, m)
                        x2, y2, z2, n2x, n2y, n2z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V, n, e, m)
                        x3, y3, z3, n3x, n3y, n3z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e, m)
                        x4, y4, z4, n4x, n4y, n4z = sqEllipsoid(x, y, z, a1, a2, a3, U, V + dV, n, e, m)

                        o.extend([COLOR, r, g, b, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                        o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                        o.extend([COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([COLOR, r, g, b, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                        o.extend([COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])

                        # Update variables for next loop */
                        V += dV
                # Update variables for next loop */
                U += dU
        o.append(END)
        return o

def makeEllipsoid(x, y, z, a1, a2, a3, m = array([[1,0,0],[0,1,0],[0,0,1]])):
	return makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, \
                                             1.0, 1.0, \
                                             -math.pi / 2, math.pi / 2, -math.pi, math.pi, \
                                             200, 200, m = m)

###########################################################################
#
#  cgo_arrow
#
###########################################################################
def cgo_arrow(atom1='pk1', atom2='pk2', radius=0.5, gap=0.0, hlength=-1, hradius=-1,
              color='blue red', name=''):
    '''
DESCRIPTION

    Create a CGO arrow between two picked atoms.

ARGUMENTS

    atom1 = string: single atom selection or list of 3 floats {default: pk1}

    atom2 = string: single atom selection or list of 3 floats {default: pk2}

    radius = float: arrow radius {default: 0.5}

    gap = float: gap between arrow tips and the two atoms {default: 0.0}

    hlength = float: length of head

    hradius = float: radius of head

    color = string: one or two color names {default: blue red}

    name = string: name of CGO object
    '''
    from chempy import cpv

    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)

    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color



    color1 = [x for x in cmd.get_color_tuple(color1)]
    color2 = [x for x in cmd.get_color_tuple(color2)]


    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
          [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
          [1.0, 0.0]


    return(obj)

    # Adding the above return, we will create one cgo object for all the arrows instead
    # of each cgo object for each arrow

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, name)

cmd.extend('cgo_arrow', cgo_arrow)


