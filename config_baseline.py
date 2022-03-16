# Config file for baseline simulations (for testing purposes)
import numpy as np
import math

## ---- Time settings ----
Tmin  = 0;    # Initial time
Tmax  = 5000; # Final time in seconds
skipT = 30;  # Number of seconds to skip storing results
# computation happens after every dt but only limited data is stored
n = 0;       # First Step n = 0 (n counts the current stored frame)
t = Tmin;    # First time t = Tmin
nframe = 0;  # First frame = 0 (nframe is the total no of stored frames)
T_arr = np.arange(0,5001,30)    #T_arr(nframe) = 0, this translation may be incorrect; # T_arr is the time array of stored frames
#arange
## ---- Viscous phenomenon flags ----
IsTopSpongeLayer = 0; # flag to include a sponge layer on top (20 km gives good results)
IsViscosity = 1;# flag to solve for molecular viscosity
IsConduction = 1; # flag to solve for thermal conduction  
IsDiffusionImplicit = 0;
## ---- Domain settings ----
# note on indexing: X(row,col) --> X(z_ind, x_ind)
Xmin = -20000;
Xmax = 20000;
Zmin = 0;
Zmax = 160000;
dx = 500; # horizontal resolution
dz = 500; # vertical resolution
SpongeHeight = 20000; # sponge layer thickness in meters

ZDomainEnd = Zmax; # last value of Z for ' physically valid' domain (pre-sponge)

if IsTopSpongeLayer == 1:
    Zmax = Zmax + SpongeHeight;    # extend Z by 50 km for computations


Xdomain = Xmax-Xmin;
Zdomain = Zmax-Zmin;
# Define positions at cell centers
x_c = np.arange(Xmin-3*dx/2,Xmax+3*dx/2+1,dx);   #2 centre points outside on either boundary. x_c = Xmin-3*dx/2:dx:Xmax+3*dx/2. Translation may need correction
z_c = np.arange(Zmin-3*dz/2,Zmax+3*dz/2+1,dz);   #z_c = Zmin-3*dz/2:dz:Zmax+3*dz/2. Translation may need correction
[X,Z] = np.meshgrid(x_c,z_c);    # grid of cell centers
[J,I] = np.shape(X);  #J gives no of z levels and I gives no of x levels

## ---- CFL ----
dCFL = 0.8; # desired Courant-Friedrichs-Lewy number
difCFL = 0.2; # CFL for diffusion problem

## ---- Background atmosphere ----
g = None
R = None
P0 = None
rho0 = None
gamma = None
C = None
#global g, R, P0, rho0, gamma, C; 

# If using Earth isothermal model
from Earth_iso import Earth_isothermal
(T0,rho0,P0,R,gamma,kinvisc,thermdiffus,H,C) = Earth_isothermal(Z);
 
# If using Earth MSIS model
#[T0,rho0,P0,R,gamma,kinvisc,thermdiffus,H,C,_,_] = Earth_MSIS(Z,10,180,2020,1,0);

## ---- Background wind ----
# only horizontal wind is specified -> time invariant. Vertical wind is zero.
# will need recheck equations to subtract from rho*u and rho*w if this changes 

global wind
# Gaussian wind shear
u_max = 0;    # wind amplitude (m/s) 
u_zloc = 100000;    # z location of wind peak (m)
u_sig = 10000;    # stdev of wind profile (m)
wind = u_max*np.exp(-(Z-u_zloc)**2/(2*u_sig**2));    # wind profile, also a matrix of size X=Z. 

# linear wind shear
# wind = linspace(0,u_max,length(z_c));
# wind = repmat(wind',1,length(x_c));    # also a matrix of size X=Z

## ---- Wave forcing ----
# A lower boundary Source is simulated as Gaussian w perturbation

class forcing:
        thermal = False;     #if thermal forcing is applied
        verticalvelocity = True;    

        if thermal == True:         #forcing_thermal=True?
            amp = 100;      # amplitude (K/s), typically of form: A x Cp
            x0 = 0;  #  forcing center x location (m)
            sigmax = 500;     # forcing half width x (m)
            z0 = 20000; #  forcing center z location (m)
            sigmaz = 1000;     # forcing half width z (m)
            t0 = 1200;      # time at forcing maxima (s)
            sigmat = 600;     # forcing half width time (s)
        else:
            # parameters for vertical velocity type forcing
            amp = 0.001;      # amplitude (m/s)
            omega = 0.007;  # centered frequency
            kx = 2*math.pi / (Xmax-Xmin);    # One horizontal wavelength per domain is set (lambda_x = x domain length)
            kxx = x_c*kx;  # computing kx*x
            t0 = 1200;      # time at forcing maxima (s)
            sigmat=600;     # forcing half width time (s)

