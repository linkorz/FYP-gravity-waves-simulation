import numpy as np
import math

def Earth_isothermal(z_array):

# Earth isothermal model to provide atmospheric state
# input: z_array -> Z x X array 
# all outputs have the same dimension as z_array

    T_surface = 239;
    p_surface = 1E5; 
    rho_surface = 1.2;
    gamma = 1.4;    #constant gamma
    R = 287;
    Pr = 0.7;

    H = p_surface/(rho_surface*9.8); # scaled height

    p = p_surface*math.exp(-z_array/H);
    rho = rho_surface*math.exp(-z_array/H);
    T = T_surface*np.ones(np.shape(z_array));
    R = R*np.ones(np.shape(z_array));
    C = math.sqrt(gamma*R*T);
    C = C*np.ones(np.shape(z_array));
    H = H*np.ones(np.shape(z_array));
    gamma = gamma*np.ones(np.shape(z_array));
    #N = sqrt(gamma-1)*9.8/C;
    kvisc = (3.5e-7)*(T**0.69)/rho; # Banks & Kockarts, 1973
    thermdiffus = kvisc/Pr;

    return T,rho,p,R,gamma,kvisc,thermdiffus,H,C