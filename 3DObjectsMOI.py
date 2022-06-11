
import sympy as smp
from scipy import integrate
import matplotlib.pyplot as plt
import math
import numpy as np

class cylinder:
    def __init__(self,mass,length,radius,offx,offy,offz):
        self.mass = mass
        self.length = length
        self.radius = radius
        self.offx = offx
        self.offy = offy
        self.offz = offz

    def cylinderMOI(self,mass,length,radius,offx,offy,offz):
        I = np.zeros((3,3),dtype=np.float64)
        self.density = mass/(math.pi*(pow(radius[1],2)-pow(radius[0],2))*length)

        # Ixx
        theta = smp.symbols('theta')
        z = smp.symbols('z')
        r = smp.symbols('r')
        fsym = z**2*r+r**3*(smp.sin(theta))**2
        fint = smp.lambdify([z,r,(smp.sin(theta))],fsym)
        fint = lambda z,r,theta:z**2*r+r**3*(smp.sin(theta))**2
        resxx,errxx = integrate.tplquad(fint,0,2*math.pi,lambda theta:radius[0],lambda theta:radius[1],lambda r,theta:-length/2,lambda r,theta:length/2)
        Ixx = (resxx*self.density)+(pow(offy,2)+pow(offz,2))*mass
        I[0][0] = Ixx

        # Iyy
        theta = smp.symbols('theta')
        z = smp.symbols('z')
        r = smp.symbols('r')
        fsym = z**2*r+r**3*(smp.cos(theta))**2
        fint = smp.lambdify([z,r,(smp.cos(theta))],fsym)
        fint = lambda z,r,theta:z**2*r+r**3*(smp.cos(theta))**2
        resyy,erryy = integrate.tplquad(fint,0,2*math.pi,lambda theta:radius[0],lambda theta:radius[1],lambda r,theta:-length/2,lambda r,theta:length/2)
        Iyy = (resyy*self.density)+(pow(offz,2)+pow(offx,2))*mass
        I[1][1] = Iyy

        # Izz
        prod = 2*math.pi*length*self.density
        f = lambda r:r**3
        rlower = 0
        rupper = radius
        reszz, errzz = integrate.quad(f,radius[0],radius[1])
        Izz = (prod*reszz)+(pow(offx,2)+pow(offy,2))*mass
        I[2][2] = Izz  

        print(I)
        return I



class cube:
    def __init__(self,mass,length,offx,offy,offz):
        self.mass = mass
        self.length = length
        self.offx = offx
        self.offy = offy
        self.offz = offz
    def cubeMOI(self,mass,length,offx,offy,offz):
        I = np.zeros((3,3),dtype=np.float64)
        self.density = mass/(length[0]*length[1]*length[2])

        # Ixx
        fint = lambda x,y,z: 1*(y**2+z**2)
        resxx, errxx = integrate.tplquad(fint,-length[2]/2,length[2]/2,lambda z:-length[1]/2,length[1]/2,lambda y,z:-length[0]/2,lambda y,z:length[0]/2)
        Ixx = (resxx*self.density)+(pow(offy,2)+pow(offz,2))*mass
        I[0][0]=Ixx

        # Iyy
        fint = lambda x,y,z: 1*(x**2+z**2)
        resyy, erryy = integrate.tplquad(fint,-length[2]/2,length[2]/2,lambda z:-length[1]/2,length[1]/2,lambda y,z:-length[0]/2,lambda y,z:length[0]/2)
        Iyy = (resyy*self.density)+(pow(offx,2)+pow(offz,2))*mass
        I[1][1]=Iyy

        # Izz
        fint = lambda x,y,z: 1*(x**2+y**2)
        reszz, errzz = integrate.tplquad(fint,-length[2]/2,length[2]/2,lambda z:-length[1]/2,length[1]/2,lambda y,z:-length[0]/2,lambda y,z:length[0]/2)
        Izz = (reszz*self.density)+(pow(offx,2)+pow(offy,2))*mass
        I[2][2]=Izz

        print(I)
        return I


print("Select the operation")
print("1) MOI of single part")
print("2) MOI of a composite body")
operation = int(input("Enter the operation number: "))

while (operation==1):
    print("Which object do you want: ")
    print("1) Cylinder")
    print("2) Cube/cuboid")
    object = int(input("Enter your choice: "))
    if (object==1):
        mass = float(input("Enter the mass of the cylinder: "))
        length = float(input("Enter the length of the cylinder: "))
        radius1 = float(input("Enter the inner radius of the cylinder: "))
        radius2 = float(input("Enter the outer radius of the cylinder: "))
        radius = [radius1,radius2]
        offsetx = float(input("Enter x offset distance: "))
        offsety = float(input("Enter y offset distance: "))
        offsetz = float(input("Enter z offset distance: "))
        cyl = cylinder(mass,length,radius,offsetx,offsety,offsetz)
        I = cyl.cylinderMOI(mass,length,radius,offsetx,offsety,offsetz)
        
    elif (object==2):
        mass = float(input("Enter the mass of the cube: "))
        length1 = float(input("Enter the length 1(Length 1 is along the x axis): "))
        length2 = float(input("Enter the length 2(Length 2 is along the y axis): "))
        length3 = float(input("Enter the length 3(Length 3 is along the z axis): "))
        length = [length1,length2,length3]
        offsetx = float(input("Enter x offset distance: "))
        offsety = float(input("Enter y offset distance: "))
        offsetz = float(input("Enter z offset distance: "))
        cube = cube(mass,length,offsetx,offsety,offsetz)
        I = cube.cubeMOI(mass,length,offsetx,offsety,offsetz)

    print("Another object?")
    print("1) Yes")
    print("2) No")
    con = int(input("Enter your choice"))
    if (con==1):
        operation = 1
    elif (con==2):
        operation = 0


while (operation==2):
    print("Which body do you want first: ")
    print("1) Cylinder")
    print("2) Cube/cuboid")
    comp = int(input("Enter the first body: ")) 
    





# mass = float(input("Enter the mass of the cube: "))
# length1 = float(input("Enter the length 1(Length 1 is along the x axis): "))
# length2 = float(input("Enter the length 2(Length 2 is along the y axis): "))
# length3 = float(input("Enter the length 3(Length 3 is along the z axis): "))

# length = [length1,length2,length3]
# offsetx = float(input("Enter x offset distance: "))
# offsety = float(input("Enter y offset distance: "))
# offsetz = float(input("Enter z offset distance: "))

# I = cylinder(mass,length,offsetx,offsety,offsetz)




