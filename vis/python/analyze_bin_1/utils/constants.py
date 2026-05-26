#CONSTANTS
G_cgs                       = 6.67e-8
kB_cgs                      = 1.38e-16
mu                          = 0.6
mp_cgs                      = 1.67e-24
M_solar                     = 2e33

gamma                       = 5.0/3.0
s_Myr                       = 3.154e13      #s
Msolar_per_kpc_sq_per_yr    = 6.7e-18       #g cm^(-2)s^(-1)


#Code units
# length_cgs                  = 3.086e+21     # 1 kpc  'L0'
# mass_cgs                    = 4.907e+40     # 1 (kpc^3)*(1mp/cm^3)   'L0^3*rho0'
# time_cgs                    = 3.0857e+13    # 0.97 Myr; v0 = 1000 km/s

length_cgs  = 3.086e+21         # 1 kpc  'L0'
mass_cgs    = 4.907e+40         # 1 (kpc^3)*(1mp/cm^3)   'L0^3*rho0'
time_cgs    = 3.086e+13         # 0.97 Myr; v0 = 1000 km/s

rho_cgs                     = mass_cgs/(length_cgs**3)
v_cgs                       = length_cgs/time_cgs
pres_cgs                    = rho_cgs*(v_cgs**2)

#Parameters
M_bh_cgs                    = 2e7*M_solar
R_bondi                     = 1.25

#Normalizations
Temp_norm                   = 1.0
Velr_scale                  = 1000.0
dens_scale                  = 1.0
pres_scale                  = 1.0

#Thresholds
tcool_threshold_Myr         = 1.0e3
t_in_threshold_Myr          = 1.0e5