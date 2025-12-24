#CONSTANTS
G_cgs = 6.67e-8
kB_cgs = 1.38e-16
mu = 0.6
mp_cgs = 1.67e-24
gamma = 5.0/3.0
s_Myr = 3.154e13       #s


#Code units
length_cgs  = 3.086e+21                     # 1 kpc  'L0'
mass_cgs    = 4.907996409352e+39            # 1 (kpc^3)*(0.1mp/cm^3)   'L0^3*rho0'
time_cgs    = 3.0857e+13                    # 0.97 Myr; v0 = 1000 km/s

rho_cgs     = mass_cgs/(length_cgs**3)
v_cgs       = length_cgs/time_cgs
pres_cgs    = rho_cgs*(v_cgs**2)

#Bondi Radius
R_bondi = 1.25

Temp_norm = 1.60
Velr_scale = 1000.0