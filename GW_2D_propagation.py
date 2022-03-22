import numpy as np
import matplotlib.pyplot as plt



# import appropriate simulation configuration file(s)
from config_baseline import *
from VarCoeffImplicitCNMatrix import VarCoeffImplicitCNMatrix
## Additional configuration from inputs of the config file

# note on indexing: X(row,col) --> X(z_ind, x_ind)

# managing z axis indices
LastDomainZindex = np.nonzero(np.logical_and(z_c > ZDomainEnd,z_c < ZDomainEnd+dz));         # last Z index for the 'physically valid' domain
FirstSpongeZindex = LastDomainZindex + np.array([1]);
if IsTopSpongeLayer == 0:  # if no sponge layer is implemented, just take last 2 indices out for top BCs
    LastDomainZindex = LastDomainZindex - np.array([2]);


# setting viscosity coefficient constant in sponge layer to prevent diffusion timestep being too low
# may not be a concern with implicit method
if IsTopSpongeLayer == 1:
     kinvisc[FirstSpongeZindex:,:] = kinvisc[LastDomainZindex,0];   


# model gravity to maintain hydrostatic equilibrium initially (g dimension is Z-1 x X-1)
g = (P0[1:,0]-P0[0:-1,0])/(-0.5*dz*(rho0[1:,0]+rho0[0:-1,0]));
g=np.transpose(np.tile(g,(np.size(X, axis=1)-1,1)))                                          

# ---- initial timestep calculation
dt = dCFL*np.minimum(dx,dz)/np.amax(C);   #limited by speed of sound          

# generate coefficient matrix for implicit diffusion                                     
if IsDiffusionImplicit:
    A_visc = VarCoeffImplicitCNMatrix(kinvisc[:,0],dz,I,J);         
    A_therm = VarCoeffImplicitCNMatrix(thermdiffus[:,0],dz,I,J);
else:
    A_visc = 0; # need to pass some value to the function later
    A_therm = 0;

## Simulation setup

# Initializing arrays:
    #our PDE system is: dQ/dt + dF/dx + dG/dz = S 
    #all these are zero arrays of size X in 4D 
F = np.dstack([0*X,0*X,0*X,0*X]);   #Fluxes for x-split(concantenate along 3rd dim to yield 3D array)    not very sure
G = F;                            #Fluxes for y-split
Q = F;                            #Solution Variables (rho, rho*u, rho*w, E)
S = F;                            #Source Terms
Q_save = np.zeros((J,I,4,167));     #4D array to hold Q results at certain cadence



# ---- Initial Conditions ----
P_pert=0*X;    # zero pressure perturbation
#P_pert=1000*exp(-.5*((X-0)**2)/(Xdomain/10)**2)*exp(-.5*((Z-10000/2)**2)/(Zdomain/10)**2);
Q[:,:,0] = rho0;    #rho
Q[:,:,1] = rho0*wind;    # rho*u
Q[:,:,2] = 0;             # rho*w (forcing is added in BCs)
Q[:,:,3] =(P_pert+P0*X**0)/(gamma-1)+0.5*rho0*wind**2; # E for ideal gas



# Set domain indices for evaluations (leave out 1st and last gridcenter)
iD = np.arange(2,I)-1;   # x indices
jD = np.arange(2,J)-1;   # z indices
[ID,JD] = np.meshgrid(iD,jD)

# Store initial state
Q_save[:,:,:,nframe] = Q;   #currently nframe is 0

#Define fuctions
## ---- Boundary Conditions ----
#function Q = bc[Q,t]
def bc(Q,t):
# This function applies boundary conditions to input Q at time t and
# returns Q after modifying it

# refer to Leveque's textbook for details on BCs

    global g, R, P0, rho0, gamma, C; 
    global wind, forcing;
    
    # ---- Bottom ----
    # Outflow condition: Zero order extrapolation + scaling to account for
    # stratification (refer Leveque)
    Q[0:2,:,0] = rho0[0:2,:]+(Q[2,:,0]-rho0[2,:])*(rho0[0:2,:]/rho0[2,:])**(0.5); # rho
    Q[0:2,:,1] = rho0[0:2,:]*wind[0:2,:]+(Q[2,:,1]-rho0[2,:]*wind[2,:])*(rho0[0:2,:]/rho0[2,:])**(0.5); #rho*u
    
    # Adding bottom forcing for rho*w
    if (forcing.verticalvelocity==False):   # i.e. if no vertical velocity forcing, use reflective BC for rho*w at domain bottom
        Q[0:2,:,2] = -Q[2,:,2]*(rho0[0:2,:]/rho0[2,:])**(0.5); 
    else: # enforce vertical velocity forcing
        w = forcing.amp*np.cos(forcing.omega*(t-forcing.t0)-forcing.kxx)*np.exp(-(t-forcing.t0)**2/(2*forcing.sigmat**2)); 
        #w = Tsunami_forcing(t); 
        Q[0:2,:,2] = w*rho0[0:2,:];
    
    # bottom for E
    Q[0:2,:,3] = P0[0:2,:]/(gamma[0:2,:]-1)+(Q[2,:,3]-P0[2,:]/(gamma[2,:]-1)-0.5*rho0[2,:]*wind[2,:]**2)*(rho0[0:2,:]/rho0[2,:])**(0.5)+0.5*rho0[0:2,:]*wind[0:2,:]**2;
    
    # ---- Top ----
    # open top (outflow) for all 4 quantities STOPPED HERE FOR INDEXING CHECK
    Q[-2:,:,0] = rho0[-2:,:]+(Q[-3,:,0]-rho0[-3,:])*(rho0[-2:,:]/rho0[-3,:])**(0.5);
    Q[-2:,:,1] = rho0[-2:,:]*wind[-2:,:]+(Q[-3,:,1]-rho0[-3,:]*wind[-3,:])*(rho0[-2:,:]/rho0[-3,:])**(0.5);
    Q[-2:,:,2] = Q[-3,:,2]*(rho0[-2:,:]/rho0[-3,:])**(0.5);
    Q[-2:,:,3] = P0[-2:,:]/(gamma[-2:,:]-1)+(Q[-3,:,3]-P0[-3,:]/(gamma[-3,:]-1)-0.5*rho0[-3,:]*wind[-3,:]**2)*(rho0[-2:,:]/rho0[-3,:])**(0.5)+0.5*rho0[-2:,:]*wind[-2:,:]**2;
        
        # ---- Sides ----
        # periodic for both sides (+x and -x)
        # note that indices 1,2,end-1,end are outside computational domain
    Q[:,1,:] = Q[:,-3,:];
    Q[:,0,:] = Q[:,-4,:];
    Q[:,-1,:]  = Q[:,3,:];
    Q[:,-2,:]= Q[:,2,:]; 

    return Q

# Set BCs for first time
Q = bc(Q,0);

## ---- Flux terms ----
#function F = Fflux[Q]
def Fflux(Q):
    global g, R, P0, rho0, gamma, C;  
    
    # ensuring that gamma is of same size as Q (when half step Qs is being passed)
    a,b,_ = np.shape(Q);
      
    gamm = gamma[0:a,0:b]; # gamm is just gamma of the correct size
    
    # compute pressure from ideal gas equation
    P=np.zeros(Q.shape)
    P = (gamm-1)*(Q[:,:,3]-(0.5*(Q[:,:,1]**2+Q[:,:,2]*2)/Q[:,:,0]));
    F=np.zeros(np.shape(Q))
    F[:,:,0] = Q[:,:,1];
    F[:,:,1] = Q[:,:,1]*Q[:,:,1]/Q[:,:,0] + P;
    F[:,:,2] = Q[:,:,1]*Q[:,:,2]/Q[:,:,0];
    F[:,:,3] = (Q[:,:,3]+P)*Q[:,:,1]/Q[:,:,0];

    return F

#function G = Gflux[Q]
def Gflux(Q):
    global g, R, P0, rho0, gamma, C; 
    
    # ensuring that gamma is of same size as Q (when half step Qs is being passed)
    a,b,_ = np.shape(Q);
    gamm = gamma[0:a,0:b]; # gamm is just gamma of the correct size
    
    P = (gamm-1)*(Q[:,:,3]-(0.5*(Q[:,:,1]**2+Q[:,:,2]**2)/Q[:,:,0]));
    G=np.zeros(Q.shape)
    G[:,:,0] = Q[:,:,2];
    G[:,:,1] = Q[:,:,1]*Q[:,:,2]/Q[:,:,0];
    G[:,:,2] = Q[:,:,2]*Q[:,:,2]/Q[:,:,0]+P;
    G[:,:,3] = (Q[:,:,3]+P)*Q[:,:,2]/Q[:,:,0];

    return G
## ---- Source term ----
#function S = Source[Q,g,X,Z,t]
def Source(Q,g,X,Z,t):
    global forcing
    S=np.zeros(Q.shape)
    S[:,:,0] = 0*Q[:,:,0];
    S[:,:,1] = 0*Q[:,:,0];
    S[:,:,2] = -Q[:,:,0]*g;
    if (forcing.thermal==False): #if thermal forcing is not applied 
        S[:,:,3] = -Q[:,:,2]*g;   # no thermal forcing (just -rho*g*w)
    else:
        a,b,_ = np.shape(Q);
        S[:,:,3] = -Q[:,:,2]*g + Q[:,:,0]*forcing.amp*np.exp(-((X[0:a,0:b]-forcing.x0)**2)/(2*forcing.sigmax**2))*np.exp(-((Z[0:a,0:b]-forcing.z0)**2)/(2*forcing.sigmaz**2))*np.exp(-((t-forcing.t0)**2)/(2*forcing.sigmat**2)); #matirx power 

    return S
## ---- Viscous terms ----

#function[Q] = MolecularViscosity[kinvisc,difCFL,dt,dx,dz,jD,iD,Q,t,IsDiffusionImplicit,A_visc,I,J]
def MolecularViscosity(kinvisc,difCFL,dt,dx,dz,jD,iD,Q,t,IsDiffusionImplicit,A_visc,I,J):
    global wind
    # This function solves the diffusion equation for molecular viscosity
    # inputs: kinvic -> array of kinematic viscosity values
             # difCFL -> CFL number for solving diffusive equations (must be <0.5)
             # dt     -> time step for the main (advective) loop
             # dz     -> spatial grid resolution
             # jD     -> indices for computation in the grid
             # Q      -> Array of prognostic vars at a timestep
             # t      -> current simulation epoch time
     # outputs: Q-> updated Q after solving for viscosity and updating u and E   
     
    
    if IsDiffusionImplicit == 1: # use Implicit method
        # solve for u
        u_diff = Q[:,:,1]/Q[:,:,0] - wind; # perturbation in u
        u_diff = SolveImplicitDiffusion(u_diff,A_visc,I,J);
        
        # repeat for w
        w_diff = Q[:,:,2]/Q[:,:,0];
        w_diff = SolveImplicitDiffusion(w_diff,A_visc,I,J);
        
        # form Q
        Q[:,:,1] = (u_diff + wind)*Q[:,:,0];
        Q[:,:,2] = w_diff*Q[:,:,0];
        
        # apply BCs to new Q
        Q = bc(Q,t);
        
    else:
        # First calculating the number of sub-steps for integration of diffusion equation 
        max_visc = np.amax(kinvisc);  #max value of viscosity in the domain
        N_substeps = np.ceil(dt*max_visc/(difCFL*min(dx,dz)**2));   #no of substeps required to solve diffusion equation... 
        #... based on Von Neumann Number (dt = N_substeps x dt_sub)

        # Main Substepping Loop
        for m in np.arange(1,N_substeps):
     
            # Substep Timestep
            dt_sub = dt/N_substeps;
        
            # x-split for diffusion equation ----
            # compute intermediate values for ease:
            Q[:,:,3] = Q[:,:,3]-0.5*(Q[:,:,1]**2+Q[:,:,2]**2)/Q[:,:,0]; # compute P/(gamma-1) 
            Q[:,:,1:3] = Q[:,:,1:3]/Q[:,:,0];  # compute velocity
            # Using an explicit scheme: forward Euler in time and centrered difference in space
            Q[JD,ID,1:3] = Q[JD,ID,0]*(Q[JD,ID,1:3]+kinvisc[JD,ID]*(dt_sub/dx**2)*(Q[JD,ID+1,1:3]-2*Q[JD,ID,1:3]+Q[JD,ID-1,1:3])); # get updated rho*u and rho*w.
            Q[JD,ID,3] = Q[JD,ID,3]+0.5*(Q[JD,ID,1]**2+Q[JD,ID,2]**2)/Q[JD,ID,0]; # get updated E
            # apply BCs
            Q = bc(Q,t);

            # z-split for diffusion equation ----
            Q[:,:,3] = Q[:,:,3]-0.5*(Q[:,:,1]**2+Q[:,:,2]**2)/Q[:,:,0];
            Q[:,:,1:3]=Q[:,:,1:3]/Q[:,:,0];
            # Explicit method
            Q[JD,ID,1:3] = Q[JD,ID,0]*(Q[JD,ID,1:3]+kinvisc[JD,ID]*(dt_sub/dz**2)*(Q[JD+1,ID,1:3]-2*Q[JD,ID,1:3]+Q[JD-1,ID,1:3]));
            Q[JD,ID,3] = Q[JD,ID,3]+0.5*(Q[JD,ID,1]**2+Q[JD,ID,2]**2)/Q[JD,ID,0];
            # apply BCs
            Q = bc(Q,t);
         
    return Q


#function[Q] = ThermalConduction[thermdiffus,difCFL,T_ref,dt,dx,dz,x_c,z_c,jD,iD,Q,t,IsDiffusionImplicit,A_therm,I,J]
def ThermalConduction(thermdiffus,difCFL,T_ref,dt,dx,dz,x_c,z_c,jD,iD,Q,t,IsDiffusionImplicit,A_therm,I,J):
    global R, gamma, P0
    # This function solves the diffusion equation for thermal conduction
    # inputs: thermdiffus -> array of thermal diffusivity values
             # T_ref -> backgrounf temperature field
             # difCFL -> CFL number for solving diffusive equations (must be <0.5)
             # dt     -> time step for the main (advective) loop
             # dz     -> spatial grid resolution
             # jD     -> indices for computation in the grid
             # Q      -> Array of prognostic vars at a timestep
             # t      -> current simulation epoch time
     # outputs: Q-> updated Q after solving for viscosity and updating u and E   
     

    
    if IsDiffusionImplicit == 1:
        Q[:,:,3] = Q[:,:,3]-0.5*(Q[:,:,1]**2+Q[:,:,2]**2)/Q[:,:,0]; # compute P/(gamma-1) 
        T = (Q[:,:,3]*(gamma-1))/(R*Q[:,:,0]);  # compute T
        T_diff = T - T_ref; # diffusion will be applied only to deviation from reference state
        T_diff = SolveImplicitDiffusion(T_diff,A_therm,I,J);
       
        # form Q
        P = Q[:,:,0]*R*(T_diff + T_ref); # get updated P
        Q[JD,ID,3] = P[JD,ID]/(gamma[JD,ID]-1) + 0.5*(Q[JD,ID,1]**2+Q[JD,ID,2]**2)/Q[JD,ID,0]; # get updated E
        
        # apply BCs to new Q
        Q = bc(Q,t);
    else:
        # First calculating the number of sub-steps for integration of diffusion equation 
        max_diffusivity = np.amax(thermdiffus);  #max value of viscosity in the domain
        N_substeps = np.ceil(dt*max_diffusivity/(difCFL*min(dx,dz)**2));   #no of substeps required to solve diffusion equation.
        #... based on Von Neumann Number (dt = N_substeps x dt_sub)

    # Main Substepping Loop
        for m in np.arange(1,N_substeps):
        
        # Substep Timestep
            dt_sub = dt/N_substeps;
        
            # x-split for diffusion equation ----
            # compute intermediate values for ease:
            Q[:,:,3] = Q[:,:,3]-0.5*(Q[:,:,1]**2+Q[:,:,2]**2)/Q[:,:,0]; # compute P/(gamma-1) 
            T = (Q[:,:,3]*(gamma-1))/(R*Q[:,:,0]);  # compute T
            T_diff = T - T_ref; # diffusion will be applied only to deviation from reference state
            # Using an explicit scheme: forward Euler in time and centrered difference in space
            T_diff[JD,ID] = T_diff[JD,ID] + thermdiffus[JD,ID]*(dt_sub/dx**2)*(T_diff[JD,ID+1]-2*T_diff[JD,ID]+T_diff[JD,ID-1]); # get updated T.
                       
            P = Q[:,:,0]*R*(T_diff + T_ref); # get updated P
            Q[JD,ID,3] = P[JD,ID]/(gamma[JD,ID]-1) + 0.5*(Q[JD,ID,1]**2+Q[JD,ID,2]**2)/Q[JD,ID,0]; # get updated E

            # apply BCs
            Q = bc(Q,t);

            # z-split for diffusion equation ----
            Q[:,:,3] = Q[:,:,3]-0.5*(Q[:,:,1]**2+Q[:,:,2]**2)/Q[:,:,0]; # compute P/(gamma-1) 
            T = (Q[:,:,3]*(gamma-1))/(R*Q[:,:,0]);  # compute T
            T_diff = T - T_ref; # diffusion will be applied only to deviation from reference state
            # Using an explicit scheme: forward Euler in time and centrered difference in space
            T_diff[JD,ID] = T_diff[JD,ID] + thermdiffus[JD,ID]*(dt_sub/dz**2)*(T_diff[jD+1,iD]-2*T_diff[JD,ID]+T_diff[JD-1,ID]); # get updated T.
            
            P = Q[:,:,0]*R*(T_diff + T_ref); # get updated P
            Q[JD,ID,3] = P[JD,ID]/(gamma[JD,ID]-1) + 0.5*(Q[JD,ID,1]**2+Q[JD,ID,2]**2)/Q[JD,ID,0]; # get updated E
        
   
        # apply BCs
        Q = bc(Q,t)

    return Q


#function[u_new] = SolveImplicitDiffusion[u_old,A,I,J]
def SolveImplicitDiffusion(u_old,A,I,J):
# This equation solves the matrix equation for solving the 2D diffusion equation discretized using Crank-Nicholson

# INPUTS:
# u_old -> previous timestep solution matrix (JxI)
# A -> coefficient matrix (see VarCoeffImplicitCNMatrix.m or ConstCoeffImplicitCNMatrix.m)
# I -> no. of x values
# J -> no. of z values
     b = np.reshape(np.transpose(u_old[:,:]),(-1,1)); #RHS of the matrix equation
     u_new = np.linalg.solve(A,b); # Direct method using LU factorization,matrix division. np.linalg
     u_new = np.transpose(np.reshape(u_new,(I,J))); # reshape the solution as matrix

     return u_new

##-----------------------------------------------------------------------------------------------##
## Computations (Lax-Wendroff 2-step)
[Qdim1,Qdim2,Qdim3]= np.shape(Q)
Qs=np.zeros((Qdim1-1,Qdim2-1,Qdim3))

while t < Tmax and nframe <= 126:       #last index of T_arr is 166
    
    # ---- x-split ---- (no source used in x split since our sources are height dependent)
    F=Fflux(Q); # compute flux F    
    # half-step
    Qs[JD,ID,:]=0.5*(Q[JD,ID,:]+Q[JD,ID+1,:])-(dt/(2*dx))*(F[JD,ID+1,:]-F[JD,ID,:]);
    Fs=Fflux(Qs); # update flux
    # full-step in x
    Q[JD,ID,:]=Q[JD,ID,:]-(dt/dx)*(Fs[JD,ID,:]-Fs[JD,ID-1,:]);
    # apply BCs
    Q=bc(Q,t);
    
    # z-split
    G=Gflux(Q); # compute flux G
    S=Source(0.5*(Q[JD,ID,:]+Q[JD+1,ID,:]),g[JD,ID],X,Z,t); # compute source
    # half step in z
    Qs[JD,ID,:]=0.5*(Q[JD,ID,:]+Q[JD+1,ID,:])-(dt/(2*dz))*(G[JD+1,ID,:]-G[JD,ID,:])+(dt/2)*S[:,:,:];
    Gs=Gflux(Qs);    # update flux
    S=Source(Qs,g,X,Z,t); # update source
    # full step in z
    Q[JD,ID,:]=Q[JD,ID,:]-(dt/dz)*(Gs[JD,ID,:]-Gs[JD-1,ID,:])+dt*0.5*(S[JD,ID,:]+S[JD-1,ID,:]);
    # apply BCs
    Q=bc(Q,t);
    
    # Solve for diffusion terms             
    # Molecular Diffusion
    if IsViscosity == 1:
        Q = MolecularViscosity(kinvisc,difCFL,dt,dx,dz,jD,iD,Q,t,IsDiffusionImplicit,A_visc,I,J);
    
    
    # Thermal conductivity
    if IsConduction == 1:
        Q = ThermalConduction(thermdiffus,difCFL,T0,dt,dx,dz,x_c,z_c,jD,iD,Q,t,IsDiffusionImplicit,A_therm,I,J);
    
    
    # Sponge layer implementation
    if IsTopSpongeLayer == 1:
        Q = MolecularViscosity(kinvisc,difCFL,dt,dx,dz,jD,iD,Q,t,IsDiffusionImplicit,A_visc,I,J);        
    
    
    # ---- Update time ----
    t=t+dt;
    n=n+1;
    
    # Update advective (main loop) timestep (adaptive)  
    dt = dCFL*np.minimum(dx,dz)/(np.amax(C)+np.amax((abs(Q[:,:,1:3]/np.dstack([Q[:,:,0]]))))); # when Q(:,:,1) -> 0, dt becomes 0 & program breaks.
    #This happens with ideal exponential density (hence, limit them to ~300
    #km). Should not be an issue with realistic density.
    if ((np.isnan(dt)) or (dt==0)):
        break;
    
    
    # Store results     
    if (round(t)%skipT==0) and (T_arr[nframe]!=round(t)):
        nframe=nframe+1
        Q_save[:,:,:,nframe]=Q
        T_arr[nframe]=t
        print(dt,n,t);
    
    


## Outputs (sim outputs are in all caps)
# compute fluid properties from saved Q data (rho, rho*u, rho*w and E)
# all outputs are only taken from indices (3:-2) in X and (3:LastDomainZindex) in Z since that is the
# computational domain, after excluding 2 ghost cells on either sides and sponge layer, if implemented.

# These values are 3D arrays (z-x-t)        
KE = np.squeeze(0.5*(Q_save[2:LastDomainZindex.item(0),2:-2,1,:]**2+Q_save[2:LastDomainZindex.item(0),2:-2,2,:]**2)/Q_save[2:LastDomainZindex.item(0),2:-2,0,:]);
P_PERT = (np.squeeze(Q_save[2:LastDomainZindex.item(0),2:-2,3,:])-KE)*np.dstack([gamma[2:LastDomainZindex.item(0),2:-2]-1])-np.dstack([P0[2:LastDomainZindex.item(0),2:-2]]);
T_PERT = P_PERT/(np.dstack([R[2:LastDomainZindex.item(0),2:-2]])*np.squeeze(Q_save[2:LastDomainZindex.item(0),2:-2,0,:]));
U = np.squeeze(Q_save[2:LastDomainZindex.item(0),2:-2,1,:]/Q_save[2:LastDomainZindex.item(0),2:-2,0,:]);    #horiz. wind (not perturbation)
W = np.squeeze(Q_save[2:LastDomainZindex.item(0),2:-2,2,:]/Q_save[2:LastDomainZindex.item(0),2:-2,0,:]);    # vertical wind (not perturbation) 

SCALING_FACTOR = np.sqrt(rho0[2:LastDomainZindex.item(0),2:-2]/rho0[2,2:-2]); # an 2d Z-X matrix
Z_KM = z_c[2:LastDomainZindex.item(0)]/1000; # grid center arrays for plotting the computational domain
X_KM = x_c[2:-2]/1000;

## Optional plots
plt.figure # wave travel plot for x in the middle of domain
plt.contourf(T_arr,Z_KM,W[:,np.size(X_KM)//2-1,:]*np.column_stack([SCALING_FACTOR[:,np.size(X_KM/2)-1]]))
plt.xlabel('time (s)')
plt.ylabel('z (km)')

# BVf = BV(g(:,1),C(:,1),gamma(:,1));
# Ri = RichardsonNumber(BVf,LastDomainZindex,U(:,length(X_KM)/2,500),500);
# figure
# plot(Ri,Z_KM,'LineWidth',2)
# xlim([-1 1])


