
import sympy as smp
from scipy import integrate
import matplotlib.pyplot as plt
import math
import numpy as np

def cylinder(mass,length,radius,offx,offy,offz):
    I = np.zeros((3,3),dtype=np.float64)
    density = mass/(math.pi*(pow(radius[1],2)-pow(radius[0],2))*length)

    # Ixx
    theta = smp.symbols('theta')
    z = smp.symbols('z')
    r = smp.symbols('r')
    fsym = z**2*r+r**3*(smp.sin(theta))**2
    fint = smp.lambdify([z,r,(smp.sin(theta))],fsym)
    fint = lambda z,r,theta:z**2*r+r**3*(smp.sin(theta))**2
    resxx,errxx = integrate.tplquad(fint,0,2*math.pi,lambda theta:radius[0],lambda theta:radius[1],lambda r,theta:-length/2,lambda r,theta:length/2)
    Ixx = (resxx*density)+(pow(offy,2)+pow(offz,2))*mass
    I[0][0] = Ixx

    # Iyy
    theta = smp.symbols('theta')
    z = smp.symbols('z')
    r = smp.symbols('r')
    fsym = z**2*r+r**3*(smp.cos(theta))**2
    fint = smp.lambdify([z,r,(smp.cos(theta))],fsym)
    fint = lambda z,r,theta:z**2*r+r**3*(smp.cos(theta))**2
    resyy,erryy = integrate.tplquad(fint,0,2*math.pi,lambda theta:radius[0],lambda theta:radius[1],lambda r,theta:-length/2,lambda r,theta:length/2)
    Iyy = (resyy*density)+(pow(offz,2)+pow(offx,2))*mass
    I[1][1] = Iyy

    # Izz
    prod = 2*math.pi*length*density
    f = lambda r:r**3
    rlower = 0
    rupper = radius
    reszz, errzz = integrate.quad(f,radius[0],radius[1])
    Izz = (prod*reszz)+(pow(offx,2)+pow(offy,2))*mass
    I[2][2] = Izz  

    print(I)
    return I



mass = float(input("Enter the mass: "))
length = float(input("Enter the length: "))
radius1 = float(input("Enter the inner radius: "))
radius2 = float(input("Enter the outer radius: "))
radius = [radius1,radius2]
offsetx = float(input("Enter x offset distance: "))
offsety = float(input("Enter y offset distance: "))
offsetz = float(input("Enter z offset distance: "))

I = cylinder(mass,length,radius,offsetx,offsety,offsetz)




