import numpy as np
import matplotlib.pyplot as plt
import math 
import pandas as pd
from PIL import Image
import random
from matplotlib.patches import Circle
def plot_hex(x,y):
    plt.plot(x,y)
def circlebased(rmax,L_RVE,Vf):
    #In-put for the cricel packing in the RVE. #Can be provided in the main file and input to the fucntion

    IFD=0.01*rmax #Important paramter to achieve packing fraction

    fn=1 #fn indicates the fiber number #randomly choose first fiber
    X=random.randint(0,L_RVE) #initate array to store the fiber centers (X,Y) & radius (R)
    Y=random.randint(0,L_RVE) 
    R=rmax
    Acc_fibers = np.empty((0, 3))
    Acc_fibers = np.append(Acc_fibers, np.array([[X,Y,R]]), axis=0)
    while fn<=100:
        fiber_found=0
        if fn==1:
            Xc=X; Yc=Y; Rc=R; #Subscripit C represnts the center fiber to pack next fibers
            (Xn,Yn,Rn,fiber_found)=NNA(Xc,Yc,Rc,IFD,rmax,Acc_fibers,L_RVE) #new fibers provided by the NNA
            Acc_fibers = np.append(Acc_fibers, np.array([[Xn,Yn,Rn]]), axis=0)
            fn=fn+1
        else:
            Xc=Acc_fibers[-1][0];Yc=Acc_fibers[-1][1];Rc=Acc_fibers[-1][2]
            (Xn,Yn,Rn,fiber_found)=NNA(Xc,Yc,Rc,IFD,rmax,Acc_fibers,L_RVE)
            if fiber_found==1:
                Acc_fibers = np.append(Acc_fibers, np.array([[Xn,Yn,Rn]]), axis=0)
                fn=fn+1
            else:
                while fiber_found==0:
                    for j in range(0,len(Acc_fibers)):
                        Xc=Acc_fibers[j][0];Yc=Acc_fibers[j][1];Rc=Acc_fibers[j][2]
                        (Xn,Yn,Rn,fiber_found)=NNA(Xc,Yc,Rc,IFD,rmax,Acc_fibers,L_RVE)
                        if fiber_found==1:
                            Acc_fibers = np.append(Acc_fibers, np.array([[Xn,Yn,Rn]]), axis=0)
                            fn=fn+1
                            break
    return(Acc_fibers)
def NNA(Xc,Yc,Rc,IFD,rmax,Acc_fibers,L_RVE):
    fiber_found=0
    attempt=1; max_attempts=10; #limit for NNA 
    temp_R=rmax; temp_IFD=IFD; #temporary fiber radius and IFD
    theta=0; #first search angle
    while fiber_found==0:
        theta=theta+random.randint(0,20);
        d=Rc+temp_R+1.01*temp_IFD #distance betweent the fiber centers
        (temp_X,temp_Y)=Find_center(Xc,Yc,d,theta)
        fiber_found=compatability_check(temp_X,temp_Y,temp_R,L_RVE,Acc_fibers,temp_IFD,rmax,fiber_found)
        #print(fiber_found)
        #do the check for the
            #1)Overlap of the accepted fibers
            #2)Fiber center in the domain or not
        if fiber_found==1:
            Xn=temp_X;Yn=temp_Y;Rn=temp_R;
            attempt=attempt+1
            break
        else:
            attempt=attempt+1
        if attempt>max_attempts:
            Xn=temp_X;Yn=temp_Y;Rn=temp_R;
            fiber_found=0
            break
    return (Xn,Yn,Rn,fiber_found)
def compatability_check(temp_X,temp_Y,temp_R,L_RVE,Acc_fibers,temp_IFD,rmax,fiber_found):
    Tol=0.75*rmax
    Pass=0 
    if temp_X<-Tol or temp_Y<-Tol or temp_X>=L_RVE+Tol or temp_Y>=L_RVE+Tol : #the temp_X and temp_Y are in bounds (defines the shape of the packing)
        fiber_found=0
    else:
        for i in range(0,len(Acc_fibers)):
            Xo=Acc_fibers[i][0];Yo=Acc_fibers[i][1];Ro=Acc_fibers[i][2] 
            Distance=math.sqrt((Xo-temp_X)**2+(Yo-temp_Y)**2)
            if Distance>Ro+temp_R+temp_IFD:
                Pass=Pass+1
            else:
                fiber_found=0
        if Pass==len(Acc_fibers):
            fiber_found=1
        else:
            fiber_found=0
    return (fiber_found)
def Find_center(Xc,Yc,d,theta):
    temp_X=Xc-d*np.cos(theta) #theta input in degrees
    temp_Y=Yc-d*np.sin(theta)
    return (temp_X,temp_Y) 
def plot_circles(adsorbed_x, adsorbed_y, radius, width, reference_indices=[]):
    # Plot each run
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for p in range(len(adsorbed_x)):
        if len(np.where(reference_indices == p)[0]) > 0:
            ax.add_patch(Circle((adsorbed_x[p], adsorbed_y[p]), radius[p], 
                edgecolor='black', facecolor='black'))
        else:
            ax.add_patch(Circle((adsorbed_x[p], adsorbed_y[p]), radius[p],
                edgecolor='black', facecolor='black'))

    ax.set_aspect(1.0)
    
    plt.axhline(y=0, color='k')
    plt.axhline(y=width, color='k')
    plt.axvline(x=0, color='k')
    plt.axvline(x=width, color='k')

    plt.xlabel("non-dimensional x")
    plt.ylabel("non-dimensional y")

    return (ax)
def fit_polygon(x,y,width,x_max,x_min,RR):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect(1)    
    
    for ii in range(0,len(x)):
        Poly_area=[]
        points_x=[]
        points_y=[]
        cents=[x[ii],y[ii]]
        angs= np.random.random(n_p).cumsum()
        angs = (angs - angs.min()) / angs.ptp()
        angs = (x_max - x_min)*angs + x_min
        
        for i in range(0,n_p-1):

            r_s=random.choice(RR)
            cc_x= cents[0]+r_s*np.cos(angs[i])
            cc_y=cents[1]+r_s*np.sin(angs[i])
            
            points_x.append(cc_x)
            points_y.append(cc_y)
            
        points_x.append(points_x[0])
        points_y.append(points_y[0])
        
        #ax.add_patch(plt.Circle((cents[0], cents[1]), RR[0], edgecolor='black',facecolor='none'))
        #ax.add_patch(plt.Circle((cents[0], cents[1]), RR[1], edgecolor='red',facecolor='none'))
        
        plt.plot(points_x,points_y,color='k')
        plt.fill(points_x,points_y, "k")  
        
    return(fig)
#Input data
rmax=10
rmin=8
n_p=6
# packing the circles
delta=20
L_RVE=rmax*delta
Vf=60
#
RR=[rmax, rmin]
x_min, x_max= 0, 2*math.pi
#Main
#Packing cricles
Acc_fibers=circlebased(rmax,L_RVE,Vf)
#Plotting cricles
plot_circles([i[0] for i in Acc_fibers], [i[1] for i in Acc_fibers], [i [2] for i in Acc_fibers], L_RVE, reference_indices=[])
#Fit the polygons
fit_polygon([i[0] for i in Acc_fibers],[i[1] for i in Acc_fibers],L_RVE,x_max,x_min,RR)
#
