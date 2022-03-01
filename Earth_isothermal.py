import numpy as np
import math
from config_baseline import *

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

    H_2D = np.zeros((X,Z))
    p_2D = np.zeros((X,Z))
    rho_2D = np.zeros((X,Z))
    T_2D = np.zeros((X,Z))
    R_2D = np.zeros((X,Z))
    C_2D = np.zeros((X,Z))
    gamma_2D = np.zeros((X,Z))
    kvisc_2D = np.zeros((X,Z))
    thermdiffus_2D = np.zeros((X,Z))

    H = p_surface/(rho_surface*9.8); # scaled height

    p_2D = p_surface*math.exp(-z_array/H);
    rho_2D = rho_surface*math.exp(-z_array/H);
    T_2D = T_surface*np.ones(np.shape(z_array));
    R_2D = R*np.ones(np.shape(z_array));
    C_2D = math.sqrt(gamma*R*T_2D);
    C_2D = C*np.ones(np.shape(z_array));
    H_2D = H*np.ones(np.shape(z_array));
    gamma_2D = gamma*np.ones(np.shape(z_array));
    #N = sqrt(gamma-1)*9.8/C;
    kvisc_2D = (3.5e-7)*(T_2D**0.69)/rho_2D; # Banks & Kockarts, 1973
    thermdiffus_2D = kvisc_2D/Pr;

    return T_2D,rho_2D,p_2D,R_2D,gamma_2D,kvisc_2D,thermdiffus_2D,H_2D,C_2D