
import sympy as smp
from scipy import integrate
import matplotlib.pyplot as plt
import math
import numpy as np


    
# def calcCentroid(KPA):




def cylinderMOI(mass,length,radius,offx,offy,offz):
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




def compCube(CGcube,start,lengthcube,con_point_cube,vector,lengthcyl):
    # For x coordinate of the centroid of cube/cuboid
    if vector[0][0]==0:
        CGcubex = 0
        # return self.CGcubex
    elif vector[0][0]==-1:
        CGcubex = CGcube[0][0]+start[0][0]+(lengthcube[0]/2-con_point_cube[0])
        # return self.CGcubex
    elif vector[0][0]==1:
        CGcubex = CGcube[0][0]+start[0][0]-(lengthcube[0]/2+con_point_cube[0])
        # return self.CGcubex

    # For y coordinate of the centroid of cube/cuboid
    if vector[1][0]==0:
        CGcubey = 0
        # return self.CGcubey
    elif vector[1][0]==-1:
        CGcubey = CGcube[1][0]+start[1][0]+(lengthcube[1]/2-con_point_cube[1])
        # return self.CGcubey
    elif vector[1][0]==1:
        CGcubey = CGcube[1][0]+start[1][0]-(lengthcube[1]/2+con_point_cube[1])
        # return self.CGcubey

    # For z coordinates of the centroid of cube/cuboid
    CGcubez = (lengthcyl+start[2][0])-CGcube[2][0]*con_point_cube[2]
    return CGcubex, CGcubey, CGcubez

def cubeMOI(mass,length,offx,offy,offz):
    I = np.zeros((3,3),dtype=np.float64)
    density = mass/(length[0]*length[1]*length[2])

    # Ixx
    fint = lambda x,y,z: 1*(y**2+z**2)
    resxx, errxx = integrate.tplquad(fint,-length[2]/2,length[2]/2,lambda z:-length[1]/2,length[1]/2,lambda y,z:-length[0]/2,lambda y,z:length[0]/2)
    Ixx = (resxx*density)+(pow(offy,2)+pow(offz,2))*mass
    I[0][0]=Ixx

    # Iyy
    fint = lambda x,y,z: 1*(x**2+z**2)
    resyy, erryy = integrate.tplquad(fint,-length[2]/2,length[2]/2,lambda z:-length[1]/2,length[1]/2,lambda y,z:-length[0]/2,lambda y,z:length[0]/2)
    Iyy = (resyy*density)+(pow(offx,2)+pow(offz,2))*mass
    I[1][1]=Iyy

    # Izz
    fint = lambda x,y,z: 1*(x**2+y**2)
    reszz, errzz = integrate.tplquad(fint,-length[2]/2,length[2]/2,lambda z:-length[1]/2,length[1]/2,lambda y,z:-length[0]/2,lambda y,z:length[0]/2)
    Izz = (reszz*density)+(pow(offx,2)+pow(offy,2))*mass
    I[2][2]=Izz

    print(I)
    return I

print("-------------------------------------------------------------")
print("Select the operation")
print("1) MOI of single part")
print("2) MOI of a composite body\n")
operation = int(input("Enter the operation number: "))
print("-------------------------------------------------------------")

while (operation==1):
    print("Which object do you want: ")
    print("1) Cylinder")
    print("2) Cube/cuboid\n")
    object = int(input("Enter your choice: "))
    print("-------------------------------------------------------------")
    if (object==1):
        mass = float(input("Enter the mass of the cylinder: "))
        length = float(input("Enter the length of the cylinder: "))
        radius1 = float(input("Enter the inner radius of the cylinder: "))
        radius2 = float(input("Enter the outer radius of the cylinder: "))
        radius = [radius1,radius2]
        offsetx = float(input("Enter x offset distance: "))
        offsety = float(input("Enter y offset distance: "))
        offsetz = float(input("Enter z offset distance: \n"))
        I = cylinderMOI(mass,length,radius,offsetx,offsety,offsetz)
        print("-------------------------------------------------------------")
        
    elif (object==2):
        mass = float(input("Enter the mass of the cube: "))
        length1 = float(input("Enter the length 1(Length 1 is along the x axis): "))
        length2 = float(input("Enter the length 2(Length 2 is along the y axis): "))
        length3 = float(input("Enter the length 3(Length 3 is along the z axis): "))
        length = [length1,length2,length3]
        offsetx = float(input("Enter x offset distance: "))
        offsety = float(input("Enter y offset distance: "))
        offsetz = float(input("Enter z offset distance: \n"))
        I = cubeMOI(mass,length,offsetx,offsety,offsetz)
        print("-------------------------------------------------------------")

    print("Another object?")
    print("1) Yes")
    print("2) No\n")
    con = int(input("Enter your choice: "))
    print("-------------------------------------------------------------")
    if (con==1):
        operation = 1
    elif (con==2):
        operation = 0


while (operation==2):
    print("Enter the number of objects: ")
    numObj = int(input("Number of objects: "))
    print("-------------------------------------------------------------")
    #  declaring the parameters here
    KPA = np.zeros((numObj,8),np.float64)
    # par2 = np.zeros((numObj,1))
    # par3 = np.zeros((numObj,1))
    # par4 = np.zeros((numObj,1))
    # par5 = np.zeros((numObj,1))
    # par6 = np.zeros((numObj,1))
    # par7 = np.zeros((numObj,1))
    # par8 = np.zeros((numObj,1))
    print("Enter the base")
    start = np.zeros((3,1),np.float64)
    for i in range(0,3):
        start[i][0] = float(input("Enter X" + str(i+1) + ":")) 

    print("-------------------------------------------------------------")
    print("Which body do you want first: ")
    print("1) Cylinder")
    print("2) Cube/cuboid")
    comp = int(input("Enter the first body: \n")) 
    print("-------------------------------------------------------------") 
    if (comp==1):
        masscyl = float(input("Enter the mass of the cylinder: "))
        lengthcyl = float(input("Enter the length of the cylinder: "))
        radius1 = float(input("Enter the inner radius of the cylinder: "))
        radius2 = float(input("Enter the outer radius of the cylinder: "))
        radius = [radius1,radius2]
        CGcyl = np.zeros((3,1),np.float64)
        CGcyl[0][0] = start[0][0]
        CGcyl[1][0] = start[1][0]
        CGcyl[2][0] = (lengthcyl/2)+start[2][0]
        KPA[0][0] = lengthcyl
        KPA[0][1] = radius[0]
        KPA[0][2] = radius[1]
        KPA[0][3] = masscyl
        KPA[0][4] = 1
        KPA[0][5] = CGcyl[0][0]
        KPA[0][6] = CGcyl[1][0]
        KPA[0][7] = CGcyl[2][0]
        
        
        
        print("-------------------------------------------------------------")
        print("Enter the point of contact on the cylinder:")
        con_point_cyl = np.zeros((3,1),np.float64)
        for i in range(0,3):
            con_point_cyl[i][0] = float(input("Enter X" + str(i+1) + ":")) 

        print("-------------------------------------------------------------")

        print("Which body do you want next?: ")
        print("1) Cylinder")
        print("2) Cube/cuboid")
        comp = int(input("Enter the next body: \n"))  
        print("-------------------------------------------------------------")
        if (comp==2):
            masscube = float(input("Enter the mass of the cube/cuboid: "))
            length1 = float(input("Enter the length 1(Length 1 is along the x axis): "))
            length2 = float(input("Enter the length 2(Length 2 is along the y axis): "))
            length3 = float(input("Enter the length 3(Length 3 is along the z axis): \n"))
            lengthcube = [length1,length2,length3]
            print("-------------------------------------------------------------")
            print("Enter the vector for the cube/cuboid:")
            vector = np.zeros((3,1))
            for i in range(0,3):
                vector[i][0] = int(input("Enter X" + str(i+1) + ":")) 
                # print(type(vector[i][0]))
                # while vector[i][0]!=-1 or vector[i][0]!=1 or vector[i][0]!=0:
                #     print("The vector should be either 1,-1 or 0")
                #     vector[i][0] = int(input("Enter X" + str(i+1) + ":")) 
            print("-------------------------------------------------------------")
            CGcube = np.zeros((3,1),np.float64)
            CGcube[0][0] = (lengthcube[0]/2)*vector[0][0]
            CGcube[1][0] = (lengthcube[1]/2)*vector[1][0]
            CGcube[2][0] = (lengthcube[2]/2)*vector[2][0]
            KPA[1][0] = lengthcube[0]
            KPA[1][1] = lengthcube[1]
            KPA[1][2] = lengthcube[2]
            KPA[1][3] = masscube
            KPA[1][4] = 2
            # KPA = keepaccount(param)
            print("Enter the point of contact on the cube/cuboid:")
            con_point_cube = np.zeros((3,1),np.float64)
            for i in range(0,3):
                con_point_cube[i][0] = float(input("Enter X" + str(i+1) + ":")) 
            # print(CGcube)
            print("-------------------------------------------------------------")
            CGcubex,CGcubey,CGcubez = compCube(CGcube,start,lengthcube,con_point_cube,vector,lengthcyl)
            CGcube = np.zeros((3,1),np.float64)
            CGcube[0][0] = CGcubex
            CGcube[1][0] = CGcubey
            CGcube[2][0] = CGcubez
            KPA[1][5] = CGcube[0][0]
            KPA[1][6] = CGcube[1][0]
            KPA[1][7] = CGcube[2][0]
            # print(CGcube)
            # print(KPA)





    





# mass = float(input("Enter the mass of the cube: "))
# length1 = float(input("Enter the length 1(Length 1 is along the x axis): "))
# length2 = float(input("Enter the length 2(Length 2 is along the y axis): "))
# length3 = float(input("Enter the length 3(Length 3 is along the z axis): "))

# length = [length1,length2,length3]
# offsetx = float(input("Enter x offset distance: "))
# offsety = float(input("Enter y offset distance: "))
# offsetz = float(input("Enter z offset distance: "))

# I = cylinder(mass,length,offsetx,offsety,offsetz)




