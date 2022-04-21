module Recfast
    using StaticArrays
    using DocStringExtensions
    using Parameters
    using DifferentialEquations
    using QuadGK

    include("params.jl")

    include("variation_constants.jl")
    include("constants.jl")

    include("h2.jl")
    include("cosmology.jl")
    include("recombination.Recfast.jl")
    include("evalode.jl")

    export Evaluate_recombination
    export Evaluate_recombination_h2
    export Evaluate_abundances_nconst
    export Set_VFC_params!
    export Switch_VFC_off!
    export Set_artificial_rates!
    export Params
end
