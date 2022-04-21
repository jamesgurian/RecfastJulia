



#rates and notation from Galli and Palla (1998)

#H + e -> H^- + \gamma
function ch3(p, TM)
    if TM < 0
        return 0
    end
    rescaling_factor = p.rescaling.RF_aScale^2* p.rescaling.RF_mScale^(-2) * p.artificial[2]
    return 1.4e-18/1e6 * rescaling_factor * (TM/p.rescaling.fl)^.928 * exp(-TM/(p.rescaling.fl *16200))
end

#H^- + \gamma -> H + e
function ch4(p, TR)
    if TR <0
        return 0
    end
    rescaling_factor = p.rescaling.RF_aScale^5 * p.rescaling.RF_mScale * p.artificial[3]
    return 1.1e-1 * rescaling_factor * (TR/p.rescaling.fl)^2.13 * exp(-8823 * p.rescaling.fl / TR)
end

# H^- + H -> H_2 + e
function ch5(p, TM)

    if TM /p.rescaling.fl < 300
        return p.rescaling.langevin * 1.5e-9/1e6 * p.artificial[4]
    else
        return p.rescaling.langevin * 4.0e-9/1e6 * (TM/ p.rescaling.fl)^-.17 * p.artificial[4]
    end
end

#H^- + H^+ -> 2H
function ch7(p, TM)
    if TM < 0
        return 0
    end
    rescale_factor = p.rescaling.RF_aScale^-3 * p.rescaling.RF_mScale^(-3) * p.artificial[5]
    return rescale_factor * 4e-6 * (TM/p.rescaling.fl)^-0.5/1e6
end

# H + H^+ -> H_2^+ + \gamma
function ch8(p, TM)
    if TM < 0
        return 0
    end
    rescale_factor = p.rescaling.RF_mxScale^-1 * p.rescaling.RF_mScale^-1 * p.rescaling.RF_aScale^-2 *
    p.rescaling.fl^2 * p.artificial[6]
    return rescale_factor * 1.85e-23 * (TM/p.rescaling.fl)^1.8 /1e6
end

# H_2^+ + \gamma -> H + H^+
function ch9(p, TR)

    rescale_factor = p.rescaling.RF_aScale^5 * p.rescaling.RF_mScale^0.5 * p.rescaling.RF_mxScale^0.5 *
     p.artificial[7]

    return rescale_factor * 1.63e7 * exp(-32400 * p.rescaling.fl/TR)
end

# H_2^+ + H -> H_2 + H^+
function ch10(p)
   return p.rescaling.langevin * 6.4e-10/1e6 * p.artificial[8]
end

#H_2^+ + H_2 -> H_3^+ + H
function ch13(p, TM)
    return  p.rescaling.langevin * 2e-9/1e6 * p.artificial[9]
end


#H_2 + H^+ -> H_2^+ + H
function ch15(p, TM)
    if TM/p.rescaling.fl < 1e4
        return 3e-10/1e6 * p.rescaling.langevin * exp(-21050 * p.rescaling.fl /TM) * p.artificial[10]
    else
        return 1.5e-10/1e6 * p.rescaling.langevin * exp(-14000 * p.rescaling.fl /TM) * p.artificial[10]
    end
end


#H3+ + e -> H2 + H
function ch20(p, TM)
    if TM < 0
        return 0
    end
    return p.rescaling.RF_aScale^(-1) *p.rescaling.RF_mScale^(-2) * 4.6e-6/1e6 * (TM/p.rescaling.fl)^(-0.65)* p.artificial[11]
end

function h2hdis(p, NH, T)
    if NH <= 0 || T <= 0
        return 0
    end
    NH *= 1e-6
    Tv = T/p.rescaling.delta_Evib
    Ta = T/p.rescaling.fl
    if Ta > 45000
        Ta = 45000
    end
    if Tv > 45000
        Tv = 45000
    end
    cdscaling = p.rescaling.RF_aScale^8 * p.rescaling.RF_mScale^(19/4) * p.rescaling.RF_mxScale^(-7/4)
    overall_scaling = p.rescaling.RF_aScale^-1 * p.rescaling.RF_mScale^(-3/2) * p.rescaling.RF_mxScale^(-1/2)
    gcid = [-1.784239e2,
            -6.842243e1,
            4.320243e1,
            -4.633167e0,
            6.970086e1,
            4.087038e4,
            -2.370570e4,
            1.288953e2,
            -5.391334e1,
             5.315517,
            -1.973427e1,
            1.678095e4,
            -2.578611e4,
            1.482123e1,
            -4.890915e0,
            4.749030e-1,
            -1.338283e2,
            -1.164408e0,
            8.227443e-1,
            5.864073e-1,
            -2.056313e0]

    gdt = [-1.427664e2,
            4.270741e1,
            -2.027365e0,
            -2.582097e-1,
            2.136094e1,
            2.753531e4,
            -2.146779e4,
            6.034928e1,
            -2.743096e1,
            2.676150e0,
            -1.128215e1,
            1.425455e4,
            -2.312520e4,
            9.305564e0,
            -2.464009e0,
            1.985955e-1,
            7.430600e2,
            -1.174242e0,
            7.502286e-1,
            2.358848e-1,
            2.937507e0]
    rate = 0
    for cs in [gcid, gdt]
        log_gh1 = cs[1]  + cs[2] * log10(Ta) + cs[3] * log10(Ta)^2 + cs[4] * log10(Ta)^3 + cs[5] * log10(1 + cs[6]/Ta)
        log_gl1 = cs[8] + cs[9] * log10(Ta) + cs[10] * log10(Ta)^2 + cs[11] * log10(1 + cs[12]/Ta)
        log_nc1 = cs[14] + cs[15] * log10(Tv) + cs[16] * log10(Tv)^2 + cs[17] /Tv
        log_gh2 = cs[7]/Ta
        log_gl2 = cs[13]/Ta
        log_nc2 = cs[18] + log_nc1

        pv = cs[19] + cs[20] * exp(-Tv/1850) + cs[21] * exp(-Tv/440)

        rate += 10 ^ (log_gh1 - (log_gh1 - log_gl1)/(1 + (NH/ (cdscaling * 10 ^ log_nc1))^pv) + log_gh2 -
        (log_gh2 - log_gl2)/(1 + (NH/(cdscaling * 10 ^ log_nc2))^pv))
    end
    return overall_scaling * rate/1e6 * p.artificial[12]
end

## H + H + H -> H2 + H
function threeh_association(p, TM)
    if TM < 0
        return 0
    end
    return 1e-12 *p.rescaling.RF_mxScale^-1 * p.rescaling.RF_mScale^-4 * p.rescaling.RF_aScale^-4 * (6e-32 * (TM/p.rescaling.fl)^(-1/4) + 2e-31 *(TM/p.rescaling.fl)^(-1/2))* p.artificial[13]
end

# H + H + H2 -> 2H2
#From yoshida 2006 this is 1/8 the previous rate.
function h_h_h2_association(p, TM)
    if TM < 0
        return 0
    end
    return 1e-12 *p.rescaling.RF_mxScale^-1 * p.rescaling.RF_mScale^-4 * p.rescaling.RF_aScale^-4 * (1/8) *(6e-32 * (TM/p.rescaling.fl)^(-1/4) + 2e-31 *(TM/p.rescaling.fl)^(-1/2))* p.artificial[14]
end

## H + H + H+ -> H2 + H+
function h_h_hp_to_h2(p,TM)
    if TM < 0
        return 0
    end
   kda = 2.817 * (TM/p.rescaling.fl)^(-1.141)
   kcta =1.201 * (TM/p.rescaling.fl)^(-1.083)
   #statistical weight * (direct association + charge transfer) * (cm^6 to m^6)
   return 0.25* (kda + kcta)* 1e-12 * 1e-29 * p.rescaling.RF_mxScale^-1 * p.rescaling.RF_mScale^-4 * p.rescaling.RF_aScale^-4 * p.artificial[15]
end

# H + H + H+ -> H2+ + H
function h_h_hp_to_h2p(p,TM)
    if TM < 0
        return 0
    end
   kda = 3.397 * (TM/p.rescaling.fl)^(-1.051)
   kcta = 1.429 * (TM/p.rescaling.fl)^(-1.021)
   #statistical weight * (direct association + charge transfer) * (cm^6 to m^6)
   return 0.25* (kda + kcta)* 1e-12 * 1e-29 * p.rescaling.RF_mxScale^-1 * p.rescaling.RF_mScale^-4 * p.rescaling.RF_aScale^-4 * p.artificial[16]

end


function Evaluate_Recfast_System_h2(p, logz, y)
    #calc h2 here
    cosmoVals = calc_cosmoVals_h2(p, 10^logz, y)
    calc_dxe_h2(p, logz, cosmoVals)

end

function calc_cosmoVals_h2(p, z, y)
    nHTot = NH(p, z)
    fHe = p.fHe
    xHep = y[1]
    xp = y[2]
    TM = y[3]
    x_Hm = y[4]
    x_H2p = y[5]
    x_H2 = y[6]
    x_H3p = y[7]
    TR = TCMB(p, z)
    xe = max(0, xp + xHep -  x_Hm + x_H2p + x_H3p)
    xH1 = max(0,(1 - xp - x_Hm - 2 * (x_H2p + x_H2) - 3 * x_H3p))
    nH1 =   xH1 * nHTot
    np = xp * nHTot
    nH2 = x_H2 * nHTot
    nH3p = x_H3p * nHTot
    nHe1 = nHTot * fHe - nHTot * xHep
    Hz = H_z_loc(p, z)
    #The total number of particles
    nTot = (xH1 + xp + xe + x_H2p + x_Hm + x_H2 + x_H3p + fHe) * ( nHTot)

    return @SVector [nTot, nHTot, xHep, xp, TM, x_Hm, x_H2p, x_H2, TR,
    xe, xH1, nH1, np, nH2, nHe1, Hz,  x_H3p,  nH3p]
end
function calc_dxe_h2(p, logz, cosmoVals)
    nTot, nHTot, xHep, xp, TM, x_Hm, x_H2p, x_H2,  TR,
    xe, xH1, nH1, np, nH2, nHe1, Hz, x_H3p, nH3p = cosmoVals
    z = 10^logz
            #compute the adjustment to the recombination rate from these reactions


            dHm = -(z * log(10)) * (ch3(p, TM) * xe * nH1 -  ch4(p, TR) * x_Hm - ch7(p, TM) * np * x_Hm
                        -  ch5(p, TM) * x_Hm * nH1) / (Hz * (1 + z))
            dH2p = - (z * log(10)) * (ch8(p, TM) * xp * nH1 - ch9(p, TR) * x_H2p -  ch10(p) * x_H2p * nH1 -ch13(p, TM) * nH2 * x_H2p
                        + ch15(p, TM) * x_H2 * np *x_H2 + h_h_hp_to_h2p(p, TM) * nH1^2 * xp ) / (Hz * (1 + z))
            dH2 = -  (z * log(10)) * ( ch5(p, TM) * x_Hm * nH1  + ch10(p) * x_H2p *nH1 -ch13(p, TM) * nH2 * x_H2p
                    -  ch15(p, TM) * x_H2 * np+ ch20(p,TM)*nH3p *xe - h2hdis(p, nH1, TM) *x_H2 * nH1
                    + threeh_association(p, TM) * nH1^2 * xH1 + h_h_h2_association(p, TM) * nH1^2 * x_H2
                    + h_h_hp_to_h2(p,TM) * nH1^2 * xp) / ((1 + z) * Hz)
            dH3p = -(z * log(10)) * ( ch13(p, TM) * nH2 * x_H2p - ch20(p,TM)*nH3p *xe ) / ((1 + z) * Hz)

            #compute the adjustment to the recombination rate from these reactions
            dHp =  -  (z * log(10)) * (-ch7(p,TM)* np * x_Hm -ch8(p, TM) * xp * nH1 + ch9(p, TR) * x_H2p
                        +  ch10(p) * x_H2p * nH1 -  ch15(p, TM) * x_H2 * np
                        - h_h_hp_to_h2(p,TM) * nH1^2 * xp)/((1 + z) * Hz)
            recombination = calc_dxe(p, logz, [nTot, nHTot, xHep, xp, TM, TR, xe, xH1, nH1, nHe1, Hz])

    return @SVector [recombination[1], recombination[2] + dHp, recombination[3], dHm, dH2p, dH2, dH3p]


end

"""
Evaluate the recombination history, including molecular hydrogen. Note that we assume everything is completely
ionized at `zstart`. 
"""
function Evaluate_recombination_h2(p; logzstart = 5., logzend = 1, dt = 1e-5,
    reltol = 1e-10, abstol = 1e-8, dtmax = 0.5, dtmin = 1e-20)
    fHe = p.fHe
    Xe_He0 = fHe
    Xe_H0 = 1
    TM0 = p.T0 * (1 + 10^logzstart)


    X_Hm = 0
    X_H2p = 0
    X_H2 = 0

    X_H3p = 0
    y0 = @SVector [Xe_He0, Xe_H0, TM0, X_Hm, X_H2p, X_H2, X_H3p]
    function evalODE(y, k, logz)
        dy = Evaluate_Recfast_System_h2(p, logz, y)

    end
    zrange = (logzstart, logzend)
    prob = ODEProblem(evalODE, y0, zrange)
    cf(integrator) = choice_function(integrator, p)
    alg_switch = CompositeAlgorithm((Rodas4P2(), Rodas5(autodiff=false)), cf)
    sol = solve(prob, alg_switch, isoutofdomain=(u,p,t) -> any(x -> x < 0, u), dt = dt, dtmax = dtmax, reltol = reltol, abstol = abstol, dtmin = dtmin, maxiters=1e6)

    t =  10 .^sol.t
    u = sol.u
    if u[end][2] < u[end][7]
        @warn "More H3+ than H+"
    end

    function mySol(z)
        t
        u
        if z < t[1]
            return sol(log10(z))
        else
            return [u[1][1], u[1][2], u[1][3] * (1 +z)/(1 + t[1]), u[1][4], u[1][5], u[1][6], u[1][7]]
        end
    end
    

    return mySol

end
