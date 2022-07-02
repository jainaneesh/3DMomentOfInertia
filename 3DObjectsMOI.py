
import sympy as smp
from scipy import integrate
import matplotlib.pyplot as plt
import math
import numpy as np

def display(con_point_cyl1,vector):
    print("--------------------------------------------------------------------------") 
    print("Contact point on the cylinder\t\t\t:" + str(con_point_cyl1[0][0]) + "\t\t" + str(con_point_cyl1[1][0]) + "\t\t" + str(con_point_cyl1[2][0]) )
    print("Vector of the cylinder \t\t\t\t:" + str(vector[0][0]) + "\t\t" + str(vector[1][0]) + "\t\t" + str(vector[2][0]) )
    print("--------------------------------------------------------------------------") 
    
def calcCentroid(KPA,numObj):
    centroidVal = np.zeros((1,4))
    i=0
    volumeMatrix = np.zeros((numObj,1))
    prodX = np.zeros((2,1))
    prodY = np.zeros((2,1))
    prodZ = np.zeros((2,1))
    while(i<numObj):
        if (KPA[i][4]==1):
            volume = math.pi*(pow(KPA[i][2],2)-pow(KPA[i][1],2))*KPA[i][0]
            prodx = volume*KPA[i][5]
            prody = volume*KPA[i][6]
            prodz = volume*KPA[i][7]
            
        elif (KPA[i][4]==2):
            volume = KPA[i][0]*KPA[i][1]*KPA[i][2]
            prodx = volume*KPA[i][5]
            prody = volume*KPA[i][6]
            prodz = volume*KPA[i][7]
        
        prodX[i][0] = prodx
        prodY[i][0] = prody
        prodZ[i][0] = prodz
        volumeMatrix[i][0] = volume
        
        i=i+1
    sumx=0
    sumy=0
    sumz=0
    volsum=0
    for j in range(0,numObj):
        sumx = sumx + prodX[j][0]
        sumy = sumy + prodY[j][0]
        sumz = sumz + prodZ[j][0]
        volsum = volsum + volumeMatrix[j][0]

    centroidX = sumx/volsum
    centroidY = sumy/volsum
    centroidZ = sumz/volsum
    # print("-------------------------------------------------------------")
    # print("Do you want the centroid?")
    # print("1) Yes")
    # print("2) No")
    # YN = int(input("Enter your choice: "))
    # print("-------------------------------------------------------------")
    # print("Centroid:")
    # print("X : " + str(centroidX) + "\tY : " + str(centroidY) + "\tZ : " + str(centroidZ))
    # centroidVal[0][0] = YN
    # centroidVal[0][1] = centroidX
    # centroidVal[0][2] = centroidY
    # centroidVal[0][3] = centroidZ
    print(centroidX,centroidY,centroidZ)
    return centroidX,centroidY,centroidZ


def compMOI(coords,KPA,numObj):
    if (numObj==2):
        if (coords[0][0]==coords[1][0] and coords[0][1]==coords[1][1]):
            print("DO you want MOI about start?")
            print("1) Yes")
            print("2) No")
            aMOI = int(input("Enter your choice: "))
            print("-------------------------------------------------------------")
            objCount = 1
            if (aMOI==1 and objCount==1):
                # calculating offsets
                for i in range(0,numObj):
                    if (KPA[i][4]==1):
                        # Calculating offset for cylinder
                        offsetx = KPA[i][5]-coords[0][0]
                        offsety = KPA[i][6]-coords[0][1]
                        offsetz = KPA[i][7]-coords[0][2]
                        radius = [KPA[i][1],KPA[i][2]]
                        Icyl = cylinderMOI(KPA[i][3],KPA[i][0],radius,offsetx,offsety,offsetz)
                        objCount = objCount+1
                    elif (KPA[i][4]==2):
                        # Calculating offset for cube/cuboid
                        offsetx = KPA[i][5]-coords[0][0]
                        offsety = KPA[i][6]-coords[0][1]
                        offsetz = KPA[i][7]-coords[0][2]
                        length = [KPA[i][0],KPA[i][1],KPA[i][2]]
                        Icube = cubeMOI(KPA[i][3],length,offsetx,offsety,offsetz)
                        objCount = objCount+1
                Itotal = Icyl+Icube
    print(Itotal)
    return Itotal

    






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


def compCyl(CGcyl,start,con_point_cyl1,vector,lengthcyl):
    
    if (vector[2][1]==1):
        # For x coordinate of the centroid of the cylinder
        if (vector[0][0]==0):
            CGxcyl = 0
        elif (vector[0][0]==-1):
            CGxcyl = CGcyl[0][0]+(start[0][0]-con_point_cyl1[0][0])
        elif (vector[0][0]==1):
            CGxcyl = CGcyl[0][0]+(start[0][0]-con_point_cyl1[0][0])
        
        # For y coordinate of the centroid of the cylinder
        if (vector[1][0]==0):
            CGycyl = 0
        elif (vector[1][0]==-1):
            CGycyl = CGcyl[1][0]+(start[1][0]-con_point_cyl1[1][0])
        elif (vector[0][0]==1):
            CGycyl = CGcyl[1][0]+(start[1][0]-con_point_cyl1[1][0])
        
        # For z coordinate of the centroid of the cylinder
        if (vector[2][0]==0):
            CGzcyl = 0
        elif (vector[2][0]==-1):
            if (start[2][0]==0):
                CGzcyl = CGcyl[2][0]+(lengthcyl/2-con_point_cyl1[2][0])
            elif (start[2][0]<0):
                CGzcyl = (CGcyl[2][0]+start[2][0])+(lengthcyl/2-con_point_cyl1[2][0])
        elif (vector[2][0]==1):
            if (start[2][0]==0):
                CGzcyl = CGcyl[2][0]-(lengthcyl/2+con_point_cyl1[2][0])
            elif (start[2][0]>0):
                CGzcyl = (CGcyl[2][0]+start[2][0])-(lengthcyl/2+con_point_cyl1[2][0])
    return CGxcyl,CGycyl,CGzcyl
        


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
    keepCount=0
    print("Enter the number of objects: ")
    numObj = int(input("Number of objects: "))
    
    coords = np.zeros((2*numObj,3),np.float64)
    
    print("------------------------------------------------------------------------------------------------")
    KPA = np.zeros((numObj,8))    
    print("------------------------------------------------------------------------------------------------")
    print("Which body do you want first: ")
    print("1) Cylinder")
    print("2) Cube/cuboid")
    comp = int(input("Enter the first body: \n")) 
    print("------------------------------------------------------------------------------------------------")
    for ii in range(0,numObj):
       
        if (comp==1 and ii==0):
            masscyl = float(input("Enter the mass of the cylinder: "))
            lengthcyl = float(input("Enter the length of the cylinder: "))
            radius1 = float(input("Enter the inner radius of the cylinder: "))
            radius2 = float(input("Enter the outer radius of the cylinder: "))
            radius = [radius1,radius2]
            CGcyl = np.zeros((3,1),np.float64)
            print("------------------------------------------------------------------------------------------------") 
            print("Enter the first point of contact on the cylinder")
            con_point_cyl1 = np.zeros((3,1),np.float64)
            for i in range(0,3):
                con_point_cyl1[i][0] = float(input("Enter  first point of contact X" + str(i+1) + ":")) 
                coords[1][i] = con_point_cyl1[i][0]
            print("------------------------------------------------------------------------------------------------")
            print("Enter the vectors for the cylinder:")
            vector = np.zeros((3*numObj,2),np.int8)
            for iii in range(0,3):
                vector[iii][0] = int(input("Enter X" + str(iii+1) + ":"))                 
                if vector[iii][0]==-1:
                    arg1=-1
                    arg2=1
                    arg3=0
                elif vector[iii][0]==1:
                    arg1=-1
                    arg2=1
                    arg3=0
                elif vector[iii][0]==0:
                    arg1=-1
                    arg2=1
                    arg3=0
                else:
                    arg1 = vector[iii][0]
                    arg2 = vector[iii][0]
                    arg3 = vector[iii][0]
                while (arg1!=-1 or arg2!=1 or arg3!=0):
                    print("The vector should be either 1,-1 or 0")
                    vector[iii] = int(input("Enter X" + str(iii+1) + ":")) 
                    if vector[iii][0]==-1:
                        arg1=-1
                        arg2=1
                        arg3=0
                    elif vector[iii][0]==1:
                        arg1=-1
                        arg2=1
                        arg3=0
                    elif vector[iii][0]==0:
                        arg1=-1
                        arg2=1
                        arg3=0
            print("------------------------------------------------------------------------------------------------")
            print("Vertical in?")
            print("1) X")
            print("2) Y")
            print("3) Z")
            vert = int(input("Enter the axis value: "))  
            print("------------------------------------------------------------------------------------------------")
            if (vert ==1):
                vector[0][1] = 1 
            elif (vert == 2):
                vector[1][1] = 1
            elif (vert ==3):
                vector[2][1] =1
            print("------------------------------------------------------------------------------------------------")
            print("Enter the base")
            start = np.zeros((3,1),np.float64)
            for i in range(0,3):
                if (vector[2][1]==1):
                    start[i][0] = float(input("Enter  base X" + str(i+1) + ":")) 
                    coords[0][i] = start[i][0]
                    if (i==0):
                        if (vector[i][0]==0):
                        # For x coordinate of the base
                            while (start[i][0]!=0 or con_point_cyl1[i][0]!=0):
                                while (start[i][0]!=0):
                                    display(con_point_cyl1,vector)
                                    print("For X vector 0| X base should be 0")
                                    print("------------------------------------------------------------------------------------------------") 
                                    print("Enter the base")
                                    start[i][0] = float(input("Enter start point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                    coords[0][i] = start[i][0]
                                while (con_point_cyl1[i][0]!=0):
                                    display(con_point_cyl1,vector)
                                    print("For X vector 0| X base 0| X contact point 1 should be 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                    coords[1][i] = con_point_cyl1[i][0]
                        elif (vector[i][0]==1):
                            if (start[i][0]==0):
                                while (con_point_cyl1[i][0]>=0 or con_point_cyl1[i][0]<-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For X vector 1| X base 0| X contact point 1 should be less than 0 and greater than " + str(-radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                coords[1][i] = con_point_cyl1[i][0] 
                            elif (start[i][0]< -radius[1] or start[i][0]< 0): 
                                while (start[i][0]>=0 or start[i][0]<=-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For X vector 1| X base should be greater than " + str(-radius[1]) + " and smaller than 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    print("Enter the base")
                                    start[i][0] = float(input("Enter start point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                coords[0][i] = start[i][0]                            
                                while (con_point_cyl1[i][0]<-radius[1] or con_point_cyl1[i][0]>start[i][0]):
                                    display(con_point_cyl1,vector)
                                    print("For X vector 1| X base less than 0| X start should be less than " + str(start[i][0]) + " and greater than " + str(-radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":"))    
                                    print("================================================================================================")                             
                                coords[1][i] = con_point_cyl1[i][0]
                            elif (start[i][0]>0 and start[i][0]<radius[1]):
                                while(con_point_cyl1[i][0]>0 or con_point_cyl1[i][0]<-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For X vector 1| X base greater than 0| X contact point on the cylinder should be greater than " + str(-radius[1]) + " and less than 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================")                                  
                                coords[1][i] = con_point_cyl1[i][0]
                            elif (start[i][0]>radius[1]):
                                while(con_point_cyl1[i][0]>radius[1] or con_point_cyl1[i][0]<-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For X vector 1| X base greater than radius| X contact point on the cylinder should be greater than " + str(-radius[1]) + " and less than " + str(radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":"))    
                                    print("================================================================================================")                              
                                coords[1][i] = con_point_cyl1[i][0]
                        elif (vector[i][0]==-1):
                            if (start[i][0]==0):
                                while (con_point_cyl1[i][0]<=0 or con_point_cyl1[i][0]>radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For X vector -1| X base 0| X contact point 1 should be greater than 0 and less than" + str(radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================")    
                                coords[1][i] = con_point_cyl1[i][0] 
                            if (start[i][0]>radius[1] or start[i][0]> 0): 
                                while (start[i][0]<=0 or start[i][0]>=radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For X vector -1| X base should be greater than 0 and smaller than" + str(radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    print("Enter the base")
                                    start[i][0] = float(input("Enter start point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                coords[0][i] = start[i][0]                            
                                while (con_point_cyl1[i][0]>=radius[1] or con_point_cyl1[i][0]<start[i][0]):
                                    display(con_point_cyl1,vector)
                                    print("For x vector -1 x base greater than 1 x start should be less than" + str(radius[1]) +  "and greater than" + start[i][0])
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":"))     
                                    print("================================================================================================")                             
                                coords[1][i] = con_point_cyl1[i][0]
                            elif (start[i][0]<0 and start[i][0]>-radius[1]):
                                while(con_point_cyl1[i][0]<0 or con_point_cyl1[i][0]>radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For X vector -1| X base less than 0| X contact point on the cylinder should be greater than" + str(-radius[1]) + "and less than 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":"))    
                                    print("================================================================================================")                               
                                coords[1][i]= con_point_cyl1[i][0]
                            elif (start[i][0]<-radius[1]):
                                while(con_point_cyl1[i][0]>radius[1] or con_point_cyl1[i][0]<-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For X vector 1| X base greater than radius| X contact point on the cylinder should be greater than" + str(-radius[1]) + "and less than" + str(radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================")                                 
                                coords[1][i] = con_point_cyl1[i][0]
                        
                    elif (i==1):
                        if (vector[i][0]==0):
                        # For y coordinate of the base
                            while (start[i][0]!=0 and con_point_cyl1[i][0]!=0):
                                while (start[i][0]!=0):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector 0| Y base should be 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    print("Enter the base")
                                    start[i][0] = float(input("Enter start point X" + str(i+1) + ":"))
                                    print("================================================================================================")  
                                coords[0][i] = start[i][0]
                                while (con_point_cyl1[i][0]!=0):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector 0| Y base 0| Y contact point 1 should be 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                coords[1][i] = con_point_cyl1[i][0]
                        elif (vector[i][0]==1):
                            if (start[i][0]==0):
                                while (con_point_cyl1[i][0]>=0 or con_point_cyl1[i][0]<-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector 1| Y base 0| Y contact point 1 should be less than 0 and greater than" + str(-radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                coords[1][i] = con_point_cyl1[i][0] 
                            elif (start[i][0]< -radius[1] or start[i][0]< 0): 
                                while (start[i][0]>=0 or start[i][0]<=-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector 1| Y base should be greater than" + str(-radius[1]) + "and smaller than 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    print("Enter the base")
                                    start[i][0] = float(input("Enter start point X" + str(i+1) + ":")) 
                                    print("================================================================================================")  
                                coords[0][i] = start[i][0]                            
                                while (con_point_cyl1[i][0]<=-radius[1] or con_point_cyl1[i][0]>start[i][0]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector 1| Y base less than 1| Y start should be less than" + start[i][0] + "and greater than" + str(-radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":"))   
                                    print("================================================================================================")                               
                                coords[1][i] = con_point_cyl1[i][0]
                            elif (start[i][0]>0 and start[i][0]<radius[1]):
                                while(con_point_cyl1[i][0]>0 or con_point_cyl1[i][0]<-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector 1| Y base greater than 0| Y contact point on the cylinder should be greater than" + str(-radius[1]) + "and less than 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":"))      
                                    print("================================================================================================")                             
                                coords[1][i] = con_point_cyl1[i][0]
                            elif (start[i][0]>radius[1]):
                                while(con_point_cyl1[i][0]>radius[1] or con_point_cyl1[i][0]<-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector 1| Y base greater than radius| Y contact point on the cylinder should be greater than" + str(-radius[1]) + "and less than" + str(radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":"))         
                                    print("================================================================================================")                         
                                coords[1][i] = con_point_cyl1[i][0]
                        elif (vector[i][0]==-1):
                            if (start[i][0]==0):
                                while (con_point_cyl1<=0 or con_point_cyl1>radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector -1| Y base 0| Y contact point 1 should be greater than 0 and less than" + str(radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                coords[1][i] = con_point_cyl1[i][0] 
                            if (start[i][0]>radius[1] or start[i][0]> 0): 
                                while (start[i][0]<=0 or start[i][0]>=radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector -1| Y base should be greater than 0 and smaller than" + str(radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    print("Enter the base")
                                    start[i][0] = float(input("Enter start point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                coords[0][i] = start[i][0]                            
                                while (con_point_cyl1[i][0]>=radius[1] or con_point_cyl1[i][0]<start[i][0]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector -1| Y base greater than 1| Y start should be less than" + str(radius[1]) +  "and greater than" + start[i][0])
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":"))   
                                    print("================================================================================================")                                
                                coords[1][i] = con_point_cyl1[i][0]
                            elif (start[i][0]<0 and start[i][0]>-radius[1]):
                                while(con_point_cyl1[i][0]<0 or con_point_cyl1[i][0]>radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector -1| Y base less than 0| Y contact point on the cylinder should be greater than" + str(-radius[1]) + "and less than 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================")                                  
                                coords[1][i]= con_point_cyl1[i][0]
                            elif (start[i][0]<-radius[1]):
                                while(con_point_cyl1[i][0]>radius[1] or con_point_cyl1[i][0]<-radius[1]):
                                    display(con_point_cyl1,vector)
                                    print("For Y vector 1| Y base greater than radius| Y contact point on the cylinder should be greater than" + str(-radius[1]) + "and less than" + str(radius[1]))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":"))        
                                    print("================================================================================================")                         
                                coords[1][i] = con_point_cyl1[i][0]
                    elif (i==2):
                        # For z coordinate of the base
                        if (vector[2][0]==0):
                            while (start[i][0]!=0 or con_point_cyl1[i][0]!=0):
                                while (start[i][0]!=0):
                                    display(con_point_cyl1,vector)
                                    print("For Z vector 0| Z base should be 0")
                                    print("------------------------------------------------------------------------------------------------")
                                    print("Enter the base")
                                    start[i][0] = float(input("Enter start point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                    coords[0][i] = start[i][0]
                                while (con_point_cyl1[i][0]!=0):
                                    display(con_point_cyl1,vector)
                                    print("For Z vector 0| Z base 0| Z contact point 1 should be 0")
                                    print("------------------------------------------------------------------------------------------------") 
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    coords[1][i] = con_point_cyl1[i][0]
                        elif (vector[2][0]==1):
                            if (start[i][0]==0):
                                while (con_point_cyl1[i][0]>=0 or con_point_cyl1[i][0]<-lengthcyl/2):
                                   display(con_point_cyl1,vector)
                                   print("For Z vector 1| Z base 0| Z contact point should be less than 0 and greater than" + str(-lengthcyl/2)) 
                                   print("------------------------------------------------------------------------------------------------")
                                   con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                   print("================================================================================================") 
                                   coords[1][i] = con_point_cyl1[i][0] 
                            elif (start[i][0]<0 or start[i][0]<-lengthcyl/4):
                                while (start[i][0]>=0 or start[i][0]<-lengthcyl/4):
                                   display(con_point_cyl1,vector)
                                   print("For Z vector 1| Z base should be less than 0 and greater than" + str(-lengthcyl/4))
                                   print("------------------------------------------------------------------------------------------------")
                                   print("Enter the base")
                                   start[i][0] = float(input("Enter start point X" + str(i+1) + ":")) 
                                   print("================================================================================================") 
                                   coords[0][i] = start[i][0]
                                while (con_point_cyl1[i][0]>start[i][0] or con_point_cyl1[i][0]<-lengthcyl/2):
                                   display(con_point_cyl1,vector)
                                   print("For Z vector 1| Z base less than 0| Z contact point should be less than" + str(start[i][0] + " and greater than " + str(-lengthcyl/2)))
                                   print("------------------------------------------------------------------------------------------------")
                                   con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                   print("================================================================================================") 
                                   coords[1][i] = con_point_cyl1[i][0]  
                            elif (start[i][0]>0 and start[i][0]<=lengthcyl/4):
                                while (con_point_cyl1[i][0]>0 or con_point_cyl1<-lengthcyl/2):
                                    display(con_point_cyl1,vector)
                                    print("For Z vector 1| Z base greater than 0 and less than " + str(lengthcyl/4) + "| Z contact point should be less than 0 and greater than " + str(-lengthcyl/2))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                            elif (start[i][0]>lengthcyl/4):
                                while (con_point_cyl1[i][0]>lengthcyl/2 or con_point_cyl1<-lengthcyl/2):
                                    display(con_point_cyl1,vector)
                                    print("For Z vector 1| Z base greater than " + str(lengthcyl/4) + "| Z contact point should be less than " + str(lengthcyl/2) + " and greater than " + str(-lengthcyl/2))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                        elif (vector[i][0]==-1):
                            if (start[i][0]==0):
                                while (con_point_cyl1[i][0]<0 or con_point_cyl1[i][0]>lengthcyl/2):
                                    display(con_point_cyl1,vector)
                                    print("For Z vector -1| Z base 0| Z contact point should be greater than 0 and less than" + str(lengthcyl/2)) 
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                    coords[1][i] = con_point_cyl1[i][0] 
                            elif (start[i][0]>0 or start[i][0]>lengthcyl/2):
                                while (start[i][0]<=0 or start[i][0]>lengthcyl/4):
                                    display(con_point_cyl1,vector)
                                    print("For Z vector -1| Z base should be greater than 0 and less than " + str(lengthcyl/4))
                                    print("------------------------------------------------------------------------------------------------")
                                    print("Enter the base")
                                    start[i][0] = float(input("Enter start point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                    coords[0][i] = start[i][0]
                                while (con_point_cyl1[i][0]<start[i][0] or con_point_cyl1[i][0]>lengthcyl/2):
                                    display(con_point_cyl1,vector)
                                    print("For Z vector -1| Z base greater than 0| Z contact point should be greater than" + str(start[i][0]) + " and less than " + str(lengthcyl/2))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                    coords[1][i] = con_point_cyl1[i][0]  
                            elif (start[i][0]<0 and start[i][0]>=lengthcyl/4):
                                while (con_point_cyl1[i][0]<0 or con_point_cyl1>lengthcyl/2):
                                    display(con_point_cyl1,vector)
                                    print("For Z vector -1| Z base less than 0 and greater than " + str(-lengthcyl/4) + "| Z contact point should be greater than 0 and less than " + str(lengthcyl/2))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================") 
                                    coords[1][i] = con_point_cyl1[i][0] 
                            elif (start[i][0]<-lengthcyl/4):
                                while (con_point_cyl1[i][0]>lengthcyl/2 or con_point_cyl1[i][0]<-lengthcyl/2):
                                    display(con_point_cyl1,vector)
                                    print("For Z vector -1| Z base less than " + str(-lengthcyl/4) + "| Z contact point should be less than " + str(lengthcyl/2) + " and greater than " + str(-lengthcyl/2))
                                    print("------------------------------------------------------------------------------------------------")
                                    con_point_cyl1[i][0] = float(input("Enter contact point X" + str(i+1) + ":")) 
                                    print("================================================================================================")    
                                    coords[1][i] = con_point_cyl1[i][0]        
                                
            
            if (vector[2][1]==1):
                origCGxcyl = 0
                origCGycyl = 0
                origCGzcyl = (lengthcyl/2)*vector[2][0]
            CGcyl = [[origCGxcyl],[origCGycyl],[origCGzcyl]]
            CGxcyl,CGycyl,CGzcyl = compCyl(CGcyl,start,con_point_cyl1,vector,lengthcyl)
            print(CGxcyl,CGycyl,CGzcyl)
            KPA[0][0] = lengthcyl
            KPA[0][1] = radius[0]
            KPA[0][2] = radius[1]
            KPA[0][3] = masscyl
            KPA[0][4] = 1
            KPA[0][5] = CGxcyl
            KPA[0][6] = CGycyl
            KPA[0][7] = CGzcyl           
            print("------------------------------------------------------------------------------------------------")
            con_point_cyl2 = np.zeros((3,1),np.float64)
            print("Enter the contact point 2 on cylinder: ")
            for j in range(0,3):
                con_point_cyl2[j][0] = float(input("Enter the second point of contact on the cylinder X" + str(j+1) + ": "))
                if (j==0):
                    while (con_point_cyl2[j][0]<KPA[0][5]-radius[1] or con_point_cyl2[j][0]>KPA[0][5]+radius[1]):
                        print("------------------------------------------------------------------------------------------------")
                        print("X Contact point 2 should lie between " + str(KPA[0][5]-radius[1]) + " and " + str(KPA[0][5]+radius[1]))
                        print("------------------------------------------------------------------------------------------------")
                        con_point_cyl2[j][0] = float(input("Enter the second point of contact on the cylinder X" + str(j+1) + ": "))
                        flag = 0
                elif (j==1):
                    while (con_point_cyl2[j][0]<KPA[0][6]-radius[1] or con_point_cyl2[j][0]>KPA[0][6]+radius[1]):
                        if (flag==0):
                            print("Y Contact point 2 should lie between " + str(KPA[0][6]-radius[1]) + " and " + str(KPA[0][6]+radius[1]))
                            print("------------------------------------------------------------------------------------------------")
                            flag = 1
                        else: 
                            print("------------------------------------------------------------------------------------------------")
                            print("Y Contact point 2 should lie between " + str(KPA[0][6]-radius[1]) + " and " + str(KPA[0][6]+radius[1]))
                            print("------------------------------------------------------------------------------------------------")
                            con_point_cyl2[j][0] = float(input("Enter the second point of contact on the cylinder X" + str(j+1) + ": "))
                            flag = 1
                elif (j==2):
                    while (con_point_cyl2[j][0]!=KPA[0][7]+lengthcyl/2):
                        if (flag==1):
                            print("Z Contact point 2 should be " + str(KPA[0][7]+lengthcyl/2))
                            print("------------------------------------------------------------------------------------------------")
                        else:
                            print("------------------------------------------------------------------------------------------------")
                            print("Z Contact point 2 should be " + str(KPA[0][7]+lengthcyl/2))
                            print("------------------------------------------------------------------------------------------------")
                        con_point_cyl2[j][0] = float(input("Enter the second point of contact on the cylinder X" + str(j+1) + ": "))
                coords[2][j] = con_point_cyl2[j][0]


            
                
        elif (comp==2 and ii==1):
            masscube = float(input("Enter the mass of the cube/cuboid: "))
            length1 = float(input("Enter the length 1(Length 1 is along the x axis): "))
            length2 = float(input("Enter the length 2(Length 2 is along the y axis): "))
            length3 = float(input("Enter the length 3(Length 3 is along the z axis): "))
            lengthcube = [length1,length2,length3]
            print("-------------------------------------------------------------")
            print("Enter the vector for the cube/cuboid:")
            # vectorcube = np.zeros((3,1)) 
            for jj in range(3,6):
                vector[jj][0] = int(input("Enter X" + str(jj+1) + ":"))                
                if vector[jj][0]==-1:
                    arg1=-1
                    arg2=1
                    arg3=0
                elif vector[jj][0]==1:
                    arg1=-1
                    arg2=1
                    arg3=0
                elif vector[jj][0]==0:
                    arg1=-1
                    arg2=1
                    arg3=0
                else:
                    arg1 = vector[jj][0]
                    arg2 = vector[jj][0]
                    arg3 = vector[jj][0]
                while (arg1!=-1 or arg2!=1 or arg3!=0):
                    print("The vector should be either 1,-1 or 0")
                    vector[jj][0] = int(input("Enter X" + str(jj+1) + ":"))   
                    if vector[jj][0]==-1:
                        arg1=-1
                        arg2=1
                        arg3=0
                    elif vector[jj][0]==1:
                        arg1=-1
                        arg2=1
                        arg3=0
                    elif vector[jj][0]==0:
                        arg1=-1
                        arg2=1
                        arg3=0 
            con_point_cube = np.zeros((3,1),np.float64)
            for jjj in range(0,3):
                if vector[jjj+3][0]==1:
                    if KPA[0][5]==0:
                        if con_point_cyl2[jjj][0]==0:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>=0 or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than 0" )
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]>0 and con_point_cyl2[jjj][0]<=KPA[0][5]+radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>=con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than or equal to 0" )
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]<0 and con_point_cyl2[jjj][0]>=KPA[0][5]-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>=con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]>0 and KPA[0][5]<=lengthcube[0]/2:
                        if con_point_cyl2[jjj][0]==0:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>=KPA[0][5] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]>0 and con_point_cyl2[jjj][0]<=radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>=KPA[0][5]+con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]+con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]<0 and con_point_cyl2[jjj][0]>=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>=KPA[0][5]+con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]+con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]>lengthcube[0]/2 and KPA[0][5]<=lengthcube[0]/2+radius1:
                        if con_point_cyl2[jjj][0]<=-(KPA[0][5]-lengthcube[0]/2) and con_point_cyl2[jjj][0]>=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>=KPA[0][5]+con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]+con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        else:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>lengthcube[0]/2 or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than or equal to " + str(lengthcube[0]/2))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]>lengthcube[0]/2+radius1:
                        print("Enter the contact point on the cube")
                        con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                        while(con_point_cube[jjj][0]>lengthcube[0]/2 or con_point_cube[jjj][0]<-lengthcube[0]/2):
                            print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than or equal to " + str(lengthcube[0]/2))
                            con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]<0 and KPA[0][5]>=-(lengthcube[0]/2-radius1):
                        if con_point_cyl2[jjj][0]==0:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:")) 
                            while(con_point_cube[jjj][0]>=KPA[0][5] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]>0 and con_point_cyl2[jjj][0]<=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:")) 
                            while(con_point_cube[jjj][0]>=KPA[0][5]-con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]-con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]<0 and con_point_cyl2[jjj][0]>=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:")) 
                            while(con_point_cube[jjj][0]>=KPA[0][5]+con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]+con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]==-lengthcube[0]/2:
                        if con_point_cyl2[jjj][0]>0 and con_point_cyl2[jjj][0]<=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:")) 
                            while(con_point_cube[jjj][0]>KPA[0][5]-con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]-con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        else:
                            print("Incorrect or invalid contact point interface entered on cylinder!")
                    elif KPA[0][5]<-lengthcube[0]/2:
                        print("The entered values are not valid")

                elif vector[jjj+3][0]==-1:
                    if KPA[0][5]==0:
                        if con_point_cyl2[jjj][0]==0:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]<=0 or con_point_cube[jjj][0]>lengthcube[0]/2):
                                print("The contact point X on the cube should be less than or equal to " + str(lengthcube[0]/2) + " and greater than 0" )
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]>0 and con_point_cyl2[jjj][0]<=KPA[0][5]+radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]<=con_point_cyl2[jjj][0] or con_point_cube[jjj][0]>lengthcube[0]/2):
                                print("The contact point X on the cube should be less than or equal to " + str(lengthcube[0]/2) + " and greater than 0" )
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]<0 and con_point_cyl2[jjj][0]>=KPA[0][5]-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]<=con_point_cyl2[jjj][0] or con_point_cube[jjj][0]>lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than " + str(con_point_cyl2[jjj][0]) + " and less than " + str(lengthcube[0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]<0 and KPA[0][5]>=-lengthcube[0]/2:
                        if con_point_cyl2[jjj][0]==0:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]<=KPA[0][5] or con_point_cube[jjj][0]>lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than " + str(KPA[0][5]) + " and less than or equal to " + str(lengthcube[0]/2))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]>0 and con_point_cyl2[jjj][0]<=radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]<=KPA[0][5]-con_point_cyl2[jjj][0] or con_point_cube[jjj][0]>lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than " + str(KPA[0][5]-con_point_cyl2[jjj][0]) + " and less than or equal to " + str(lengthcube[0]/2))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]<0 and con_point_cyl2[jjj][0]>=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>=KPA[0][5]+con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]+con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]>lengthcube[0]/2 and KPA[0][5]<=lengthcube[0]/2+radius1:
                        if con_point_cyl2[jjj][0]<=-(KPA[0][5]-lengthcube[0]/2) and con_point_cyl2[jjj][0]>=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>=KPA[0][5]+con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]+con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        else:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                            while(con_point_cube[jjj][0]>lengthcube[0]/2 or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than or equal to " + str(lengthcube[0]/2))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]>lengthcube[0]/2+radius1:
                        print("Enter the contact point on the cube")
                        con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:"))
                        while(con_point_cube[jjj][0]>lengthcube[0]/2 or con_point_cube[jjj][0]<-lengthcube[0]/2):
                            print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than or equal to " + str(lengthcube[0]/2))
                            con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]<0 and KPA[0][5]>=-(lengthcube[0]/2-radius1):
                        if con_point_cyl2[jjj][0]==0:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:")) 
                            while(con_point_cube[jjj][0]>=KPA[0][5] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]>0 and con_point_cyl2[jjj][0]<=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:")) 
                            while(con_point_cube[jjj][0]>=KPA[0][5]-con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]-con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        elif con_point_cyl2[jjj][0]<0 and con_point_cyl2[jjj][0]>=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:")) 
                            while(con_point_cube[jjj][0]>=KPA[0][5]+con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]+con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                    elif KPA[0][5]==-lengthcube[0]/2:
                        if con_point_cyl2[jjj][0]>0 and con_point_cyl2[jjj][0]<=-radius1:
                            print("Enter the contact point on the cube")
                            con_point_cube[jjj][0] = float(input("Enter the contact point X" + str(jjj+1) + " on the cube:")) 
                            while(con_point_cube[jjj][0]>KPA[0][5]-con_point_cyl2[jjj][0] or con_point_cube[jjj][0]<-lengthcube[0]/2):
                                print("The contact point X on the cube should be greater than or equal to " + str(-lengthcube[0]/2) + " and less than " + str(KPA[0][5]-con_point_cyl2[jjj][0]))
                                con_point_cube[jjj][0] = float(input("Enter the correct contact point X" + str(jjj+1) + " on the cube:"))
                        else:
                            print("Incorrect or invalid contact point interface entered on cylinder!")
                    elif KPA[0][5]<-lengthcube[0]/2:
                        print("The entered values are not valid")






                              














                 

                    

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
            print("Enter the point of contact on the cube/cuboid:")
            con_point_cube = np.zeros((3,1),np.float64)
            for i in range(0,3):
                con_point_cube[i][0] = float(input("Enter X" + str(i+1) + ":")) 

            print("-------------------------------------------------------------")
            CGcubex,CGcubey,CGcubez = compCube(CGcube,start,lengthcube,con_point_cube,vector,lengthcyl)
            CGcube = np.zeros((3,1),np.float64)
            CGcube[0][0] = CGcubex
            CGcube[1][0] = CGcubey
            CGcube[2][0] = CGcubez
            KPA[1][5] = CGcube[0][0]
            KPA[1][6] = CGcube[1][0]
            KPA[1][7] = CGcube[2][0] 
            
        
        if (ii<numObj-1):
            print("------------------------------------------------------------------------------------------------")
            print("Which body do you want next?: ")
            print("1) Cylinder")
            print("2) Cube/cuboid")
            comp = int(input("Enter the next body: "))  
            print("------------------------------------------------------------------------------------------------")

    centroidX,centroidY,centroidZ = calcCentroid(KPA,numObj)
    Itotal = compMOI(coords,KPA,numObj)








    





# mass = float(input("Enter the mass of the cube: "))
# length1 = float(input("Enter the length 1(Length 1 is along the x axis): "))
# length2 = float(input("Enter the length 2(Length 2 is along the y axis): "))
# length3 = float(input("Enter the length 3(Length 3 is along the z axis): "))

# length = [length1,length2,length3]
# offsetx = float(input("Enter x offset distance: "))
# offsety = float(input("Enter y offset distance: "))
# offsetz = float(input("Enter z offset distance: "))

# I = cylinder(mass,length,offsetx,offsety,offsetz)




