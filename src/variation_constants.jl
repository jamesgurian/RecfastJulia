#Fundamental constant scaling



"""

Function to turn off all rescaling
$(TYPEDSIGNATURES)
"""
function Switch_VFC_off!(r::Variables_Fund_Consts)
    r.RF_aScale = r.RF_mScale = 1.0;
    r.RF_aScale = r.RF_mScale = r.fT = r.fl = r.fA2g = r.fC = r.fAp =
        r.fBp = r.fm = r.RF_mxScale = r.delta_Evib = r.langevin = 1.0
    return
end

"""
$(TYPEDSIGNATURES)
"""
function Switch_VFC_off!(p)
    Switch_VFC_off!(p.rescaling)
    return
end

"""
$(TYPEDSIGNATURES)
Function to set the parameters fields.
"""
function Set_VFC_params!(r::Variables_Fund_Consts;
    aS = 1, mS = 1, mxS = 1)

    # Reset to avoid "double dipping"
    Switch_VFC_off!(r)
    r.RF_aScale = aS
    r.RF_mScale = mS
    r.RF_mxScale = mxS

    #Rescale temp
    (r.fT = aS^2 * mS)

    # Rescale Thomson cross section and cooling
    # Extra 1/me form rho_g/me in coefficient is neglected
    (r.fC = (aS/mS)^2)

    # Rescale two photon rate for lambda values
    (r.fA2g = aS^8 * mS)

    #Rescale photon-ionization, recombination coefficients & Saha

    r.fAp = (aS/mS)^2
    r.fBp = aS^5 * mS
    r.fm = mS

    #Rescaling of wavelength to change the Ly-a channel
    (r.fl = aS^2 * mS)

    #rescale everything we need for molecular hydrogen

    r.delta_Evib = aS^2 * mS^(3/2) * mxS^(-1/2)
    r.delta_Erot = aS^2 * mS^2 /mxS
    r.langevin = aS^(-1) * mS ^(-3/2) * mxS ^(-1/2)
    return
end

"""
$(TYPEDSIGNATURES)
"""
function Set_VFC_params!(p;
    aS = 1, mS = 1, mxS = 1)

    Set_VFC_params!(p.rescaling, aS = aS, mS = mS, mxS = mxS)
end

"""
Function to set artificial reaction rates
"""
function Set_artificial_rates!(p, rates)
        p.artificial = rates
    end
