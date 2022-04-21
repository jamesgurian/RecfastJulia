

## Hubble constant, in s^-1
function H_z_loc(params, z)
    z1 = 1.0 + z
    Fnu = params.n_eff * 7/8 * (4/11)^(4/3)
    Zeq = (3.0 * (params.H0 * RF_cLight)^2 /
        (8 * pi * RF_G * RF_aRad * (1 + Fnu)) /
        params.T0^4 * params.Omega_M - 1.0)
    return params.H0 * sqrt(params.Omega_L + z1^2 *
    (params.Omega_K + z1 * params.Omega_M * (1 + z1/(1 + Zeq))))
end

#Hydrogen number density, in m^-3
function NH(params, z)
    mu_H = 1 / (1 - params.Yp)
    return 3 * (params.H0)^2 * params.Omega_B /
    (8 * pi * RF_G *RF_mHatom * mu_H) * (1 + z)^3
end

#CMB Temp, K
TCMB(params, z) = params.T0 * (1 + z)

#Calculate the relativistic contribution
#This is used to calculate the dark energy when not specified.
function calc_Orel(TCMB0, Nnu, h100)
    H0 = h100 * 100 * 1e5 /RF_Mpc
    a = RF_aRad * TCMB0^4 / RF_cLight^2
    b = 3 * H0^2 / (8 * pi * RF_G)
    return a / b * (1 + Nnu * (7/8) * (4/11)^(4/3))
end
