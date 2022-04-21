"""
The implicit solver which we need at high temperature seems to have problems
at low T, so we switch.
"""
function choice_function(integrator, p)
    if (integrator.u[3]/p.rescaling.fl < 100)
        choice = 2
    else
        choice = 1
    end
    return choice
end

"""
Evaluate the recombination history. Note that we assume everything is completely
ionized at `zstart`. 
"""
function Evaluate_recombination(p; logzstart = 4., logzend = 1., dt = 1e-5,
    reltol = 1e-10, abstol = 1e-8, dtmax = 0.5, dtmin = 1e-10)
    fHe = p.fHe
    Xe_He0 = fHe
    Xe_H0 = 1
    TM0 = p.T0 * (1 + 10^logzstart)
    y0 = [Xe_He0, Xe_H0, TM0]
    function evalODE(y, k, logz)
        dy = Evaluate_Recfast_System(p, logz, y)
    end
    zrange = (logzstart, logzend)
    prob = ODEProblem(evalODE, y0, zrange)
    cf(integrator) = choice_function(integrator, p)
    alg_switch = CompositeAlgorithm((Rodas4P2(), Rodas5(autodiff=false)), cf)
    sol = solve(prob, alg_switch, dt = dt, dtmax = dtmax, reltol = reltol, abstol = abstol, dtmin = dtmin)
    t =  10 .^sol.t
    u = sol.u

    function mySol(z)
        t
        u
        if z < t[1]
            return sol(log10(z))
        else
            return [u[1][1], u[1][2], u[1][3] * (1 +z)/(1 + t[1])]
        end
    end


    return mySol
end
