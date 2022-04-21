
"""
Struct to control the rescaling. 
"""
@with_kw mutable struct Variables_Fund_Consts
    RF_aScale::Float64 = 1.0 # alpha/alpha_r
    RF_mScale::Float64 = 1.0 # mS = me/me_r
    RF_mxScale::Float64 = 1.0 # mx/mp
    fT::Float64 = 1.0 # Temperature scale factor
    fl::Float64 = 1.0 # Wavelength scale factor
    fA2g::Float64 = 1.0 # Two photon scale factor
    fC::Float64 = 1.0 # Compton cooling scale factor
    fAp::Float64 = 1.0 # Normalization phot-coefficient scale factor
    fBp::Float64 = 1.0 # Normalization rec-coefficient scale factor
    fm::Float64 = 1.0 # Mass scale factor for saha
    delta_Evib::Float64 = 1.0 #Ratio of H2 vibrational energy levels
    delta_Erot::Float64 = 1.0
    langevin::Float64 = 1.0 #langevin rate rescaling

end

"""
Struct holding all the cosmology parameters. 
"""
@with_kw mutable struct Params
    Yp::Float64 = 0.24
    T0::Float64 = 2.725
    Omega_M::Float64 = 0.26
    Omega_B::Float64 = 0.044
    Omega_K::Float64 = 0.0


    h100::Float64 = 0.6774
    n_eff::Float64 = 3.04
    F::Float64 = 1.14

    rescaling::Variables_Fund_Consts = Variables_Fund_Consts()

    fHe = Yp / (4.0 * RF_fac_mHemH * (1 - Yp))
    H0 = h100 * 100 * 1e5 / RF_Mpc
    Omega_L::Float64 = 1.0-Omega_K - Omega_M -
    calc_Orel(T0, n_eff, h100);
    artificial::SVector{16, Float64} = @SVector ones(16)

end

Broadcast.broadcastable(m) = Ref(m)
