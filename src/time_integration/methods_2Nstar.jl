# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
using JSON

@muladd begin
#! format: noindent

# Abstract base type for time integration schemes of storage class `2N*`
abstract type SimpleAlgorithm2Nstar end

"""
Low Storage 2N* SSP s-stage explicit Runge Kutta method of order 4 
with only two storage registers
"""
struct SSPs4_2Nstar <: SimpleAlgorithm2Nstar
    NumStages::Int
    lambda::SVector
    gamma::SVector
    c::SVector

    function SSPs4_2Nstar(NumStages::Int)

        data = JSON.parsefile("src/time_integration/SSP_2Nstar_methods.json")

        lambda = SVector{2*NumStages-1, Float64}(data[string(NumStages)]["Lambda"])
        gamma = SVector{2*NumStages-1, Float64}(data[string(NumStages)]["Gamma"])
        c = SVector{NumStages, Float64}(data[string(NumStages)]["c"])

        new(NumStages, lambda, gamma, c)
    end
end

"""
Low Storage 2N* SSP 6-stage explicit Runge Kutta method of order 4 
with only two storage registers
"""
struct SSP64_2Nstar <: SimpleAlgorithm2Nstar
    lambda::SVector{11, Float64}
    gamma::SVector{11, Float64}
    c::SVector{6, Float64}

    function SSP64_2Nstar()
        lambda = SVector(1.0,
                        0.9182743963949238,
                        0.0817256036050762,
                        0.9513953767623428,
                        0.0486046232376572,
                        0.7597661646688474,
                        0.2402338353311526,
                        6.99957e-11,
                        0.9999999999300043,
                        0.2613787016800961,
                        0.7386212983199039)

        gamma = SVector(0.6104727475803353,
                        0.0590270558926983,
                        0.4166581726008912,
                        0.2935587144522069,
                        0.1135070598044559,
                        2.6817e-12,
                        0.5610214523406655,
                        4.1357e-12,
                        0.4056356699993987,
                        0.1076836983827033,
                        0.1015594372977109)

        c = SVector(0.0,
                    0.6104727475803353,
                    0.5255764822740417,
                    0.4326112211601659,
                    0.6649493052099473,
                    1.070584975166938)

        new(lambda, gamma, c)
    end
end

"""
Low Storage 2N* SSP 10-stage explicit Runge Kutta method of order 4 
with only two storage registers
"""
struct SSP104_2Nstar <: SimpleAlgorithm2Nstar
    lambda::SVector{19, Float64}
    gamma::SVector{19, Float64}
    c::SVector{10, Float64}

    function SSP104_2Nstar()
        lambda = SVector(1.0,
                         0.7461586730434318, 0.2538413269565682,
                         0.3360004731162419, 0.6639995268837581,
                         0.3871068892184111, 0.6128931107815889,
                         0.9748974615674895, 0.0251025384325105,
                         0.4736841498257829, 0.5263158501742171,
                         0.3676801369895344, 0.6323198630104656,
                         0.0441436069847376, 0.9558563930152624,
                         0.0000000000368702, 0.9999999999631298,
                         0.1468854337190545, 0.8531145662809454)

        gamma = SVector(0.5801774807348781,
                        0.2774179223775330, 0.4246909440456424,
                        0.0000000000291864, 0.3852375727963142,
                        0.0000000000306165, 0.3555867810453154,
                        0.2144917481287411, 0.0145639275107956,
                        0.0000000000032062, 0.3053566040927783,
                        0.0000000000020217, 0.3668577452121914,
                        0.0000000000016537, 0.2596556653402283,
                        0.0000000000026306, 0.1285869610319747,
                        0.0615546206861932, 0.0820758149079128)

        c = SVector(0.0,
                    0.580177480734878,
                    0.849381888003236,
                    0.949226744603282,
                    0.937361313412918,
                    0.252585824034633,
                    0.438296526814728,
                    0.644001345007665,
                    0.875228468077886,
                    1.00381542908022)

        new(lambda, gamma, c)
    end
end

"""
Low Storage 2N* SSP 14-stage explicit Runge Kutta method of order 4 
with only two storage registers
"""
struct SSP144_2Nstar <: SimpleAlgorithm2Nstar
    lambda::SVector{27, Float64}
    gamma::SVector{27, Float64}
    c::SVector{14, Float64}

    function SSP144_2Nstar()
        lambda = SVector(1.0,
                        0.8334108980877719,
                        0.1665891019122281,
                        0.1032947056381604,
                        0.8967052943618397,
                        0.1830143341065496,
                        0.8169856658934505,
                        0.1930749676770172,
                        0.8069250323229828,
                        0.2076136990573719,
                        0.7923863009426282,
                        0.9805567938107509,
                        0.0194432061892491,
                        0.3584179884135753,
                        0.6415820115864247,
                        0.2402807846126498,
                        0.7597192153873502,
                        0.1749841662040309,
                        0.8250158337959691,
                        0.0256869322991015,
                        0.9743130677008985,
                        4.62799e-11,
                        0.9999999999537201,
                        3.78555e-11,
                        0.9999999999621445,
                        0.1011277304506511,
                        0.8988722695493488)

        gamma = SVector(0.3128028073986172,
                        0.2606932685872564,
                        0.312802807353276,
                        2.42634e-11,
                        0.2804919334933156,
                        1.91822e-11,
                        0.2555554099030914,
                        1.91364e-11,
                        0.2524084154775574,
                        2.69325e-11,
                        0.2478606594852272,
                        0.1786898804726819,
                        0.0060818894809712,
                        3.0452e-12,
                        0.2006886544057519,
                        1.7358e-12,
                        0.2376423034206306,
                        1.3275e-12,
                        0.2580672689526694,
                        1.1882e-12,
                        0.1705555524739845,
                        1.3961e-12,
                        0.1020159399787256,
                        2.3342e-12,
                        0.0794125638639779,
                        0.0428049958695966,
                        0.0606593646782706)

        c = SVector(0.0,
                    0.3128028073986172,
                    0.6256056147006916,
                    0.8414758004021823,
                    0.9430290770470746,
                    1.0133621839744174,
                    1.0508349719867913,
                    0.20520337098486602,
                    0.3323434459495827,
                    0.49013000541831186,
                    0.6624322840426085,
                    0.8159719832848393,
                    0.9179879232271981,
                    0.9974004870587593)

        new(lambda, gamma, c)
    end
end

"""
Low Storage 2N* SSP 15-stage explicit Runge Kutta method of order 4 
with only two storage registers
"""
struct SSP154_2Nstar <: SimpleAlgorithm2Nstar
    lambda::SVector{29, Float64}
    gamma::SVector{29, Float64}
    c::SVector{15, Float64}

    function SSP154_2Nstar()
        lambda = SVector(1.0,
                        0.8523352138377114,
                        0.1476647861622886,
                        0.1666104720490982,
                        0.8333895279509018,
                        0.1860951577164843,
                        0.8139048422835157,
                        0.1920781653557145,
                        0.8079218346442856,
                        0.2040624560612781,
                        0.795937543938722,
                        0.978973912989085,
                        0.021026087010915,
                        0.3670798745644486,
                        0.6329201254355514,
                        0.242963372653435,
                        0.757036627346565,
                        0.1474154587865577,
                        0.8525845412134423,
                        0.0125940747063368,
                        0.9874059252936632,
                        1.419882e-10,
                        0.9999999998580118,
                        7.1378e-11,
                        0.999999999928622,
                        8.22213e-11,
                        0.9999999999177787,
                        0.0860899606051551,
                        0.9139100393948448)

        gamma = SVector(0.3001460971382099,
                        0.255825087199926,
                        0.3001460964946184,
                        8.0014e-11,
                        0.2501386142597081,
                        6.01914e-11,
                        0.2442903619006913,
                        5.71371e-11,
                        0.2424945855060537,
                        7.65756e-11,
                        0.2388975474202729,
                        0.1925980689368956,
                        0.0063108979554327,
                        8.3303e-12,
                        0.1899685054708881,
                        8.7052e-12,
                        0.2272215889077047,
                        7.6354e-12,
                        0.2247096327018077,
                        7.8622e-12,
                        0.1230520444518983,
                        1.21192e-11,
                        0.096517927584308,
                        3.1347e-11,
                        0.080219240841554,
                        0.0022301820900746,
                        0.0697445633367235,
                        0.0335174862492144,
                        0.0575463517000687)

        c = SVector(0.0,
                    0.3001460971382099,
                    0.6002921929459037,
                    0.7504158416515205,
                    0.8550574492073152,
                    0.933314168653028,
                    0.9817573346177497,
                    0.21955148203360508,
                    0.3289270570274888,
                    0.4762314188115314,
                    0.6307371784282993,
                    0.7458456717428698,
                    0.8423635992333957,
                    0.9225828400461705,
                    0.9945575853971129)

        new(lambda, gamma, c)
    end
end

"""
Low Storage 2N* SSP 5-stage explicit Runge Kutta method of order 3 
with only two storage registers
(I. Higueras, T. Roldan (2019) New Third Order Low-Storage SSP Explicit Runge–Kutta Methods) 
"""
struct SSP53_2Nstar <: SimpleAlgorithm2Nstar
    lambda::SVector{9, Float64}
    gamma::SVector{9, Float64}
    c::SVector{5, Float64}
    A::SMatrix{5, 5, Float64}
    b::SVector{5, Float64}

    function SSP53_2Nstar()
        # SSP_2N*(5,3) Butcher and Shu Osher
        A = SMatrix{5, 5}(0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
                          0.0000000000000000, 0.0000000000000000,
                          0.3772689153349139, 0.0000000000000000, 0.0000000000000000,
                          0.0000000000000000, 0.0000000000000000,
                          0.3772689153210351, 0.3772689153210356, 0.0000000000000000,
                          0.0000000000000000, 0.0000000000000000,
                          0.2273646445334142, 0.2273646445166876, 0.2273646445250521,
                          0.0000000000000000, 0.0000000000000000,
                          0.2118911342119499, 0.1277410428919109, 0.1277410428966103,
                          0.2119622635173820, 0.0000000000000000)
        A = A'
        c = Vector{Float64}(undef, 5)
        for i in 1:5
            c[i] = sum(A[i, j] for j in 1:5)
        end
        b = SVector(0.2256709383174730, 0.1170972518383502, 0.1170972518426580,
                    0.1943008917840403, 0.3458336662174786)
        lambda = SVector(1.0000000000000000,
                         0.5834464195819975, 0.4165535804180025,
                         0.3973406361263226, 0.6026593638736774,
                         0.4381666368425402, 0.5618333631574598,
                         0.0833231889500624, 0.9166768110499376)

        gamma = SVector(0.3772689153349139,
                        0.2201161978578605, 0.3772689153210356,
                        0.0000000000167269, 0.2273646445250521,
                        0.0841500913106415, 0.2119622635173820,
                        0.0314352491183084, 0.3458336662174786)

        new(lambda, gamma, c)
    end
end

# This struct is needed to fake https://github.com/SciML/OrdinaryDiffEq.jl/blob/0c2048a502101647ac35faabd80da8a5645beac7/src/integrators/type.jl#L1
mutable struct SimpleIntegrator2NstarOptions{Callback}
    callback::Callback # callbacks; used in Trixi.jl
    adaptive::Bool # whether the algorithm is adaptive; ignored
    dtmax::Float64 # ignored
    maxiters::Int # maximal number of time steps
    tstops::Vector{Float64} # tstops from https://diffeq.sciml.ai/v6.8/basics/common_solver_opts/#Output-Control-1; ignored
end

function SimpleIntegrator2NstarOptions(callback, tspan; maxiters = typemax(Int),
                                       kwargs...)
    SimpleIntegrator2NstarOptions{typeof(callback)}(callback, false, Inf, maxiters,
                                                    [last(tspan)])
end

# This struct is needed to fake https://github.com/SciML/OrdinaryDiffEq.jl/blob/0c2048a502101647ac35faabd80da8a5645beac7/src/integrators/type.jl#L77
# This implements the interface components described at
# https://diffeq.sciml.ai/v6.8/basics/integrator/#Handing-Integrators-1
# which are used in Trixi.jl.
mutable struct SimpleIntegrator2Nstar{RealT <: Real, uType, Params, Sol, F, Alg,
                                      SimpleIntegrator2NstarOptions} <: AbstractTimeIntegrator
    u::uType #
    du::uType
    u_tmp::uType
    du_tmp::uType
    t::RealT
    dt::RealT # current time step
    dtcache::RealT # ignored
    iter::Int # current number of time steps (iteration)
    p::Params # will be the semidiscretization from Trixi.jl
    sol::Sol # faked
    f::F
    alg::Alg
    opts::SimpleIntegrator2NstarOptions
    finalstep::Bool # added for convenience
end

# Forward integrator.stats.naccept to integrator.iter (see GitHub PR#771)
function Base.getproperty(integrator::SimpleIntegrator2Nstar, field::Symbol)
    if field === :stats
        return (naccept = getfield(integrator, :iter),)
    end
    # general fallback
    return getfield(integrator, field)
end

function init(ode::ODEProblem, alg::SimpleAlgorithm2Nstar;
              dt, callback::Union{CallbackSet, Nothing} = nothing, kwargs...)
    u = copy(ode.u0)
    du = zero(u)
    u_tmp = zero(u)
    du_tmp = zero(du)

    t = first(ode.tspan)
    iter = 0
    integrator = SimpleIntegrator2Nstar(u, du, u_tmp, du_tmp, t, dt, zero(dt), iter,
                                        ode.p,
                                        (prob = ode,), ode.f, alg,
                                        SimpleIntegrator2NstarOptions(callback,
                                                                      ode.tspan;
                                                                      kwargs...), false)

    # initialize callbacks
    if callback isa CallbackSet
        foreach(callback.continuous_callbacks) do cb
            throw(ArgumentError("Continuous callbacks are unsupported with the 2N storage time integration methods."))
        end
        foreach(callback.discrete_callbacks) do cb
            cb.initialize(cb, integrator.u, integrator.t, integrator)
        end
    end

    return integrator
end

# Fakes `solve`: https://diffeq.sciml.ai/v6.8/basics/overview/#Solving-the-Problems-1
function solve(ode::ODEProblem, alg::SimpleAlgorithm2Nstar;
               dt, callback = nothing, kwargs...)
    integrator = init(ode, alg, dt = dt, callback = callback; kwargs...)

    # Start actual solve
    solve!(integrator)
end

function solve!(integrator::SimpleIntegrator2Nstar)
    @unpack prob = integrator.sol

    integrator.finalstep = false

    @trixi_timeit timer() "main loop" while !integrator.finalstep
        step!(integrator)
    end # "main loop" timer

    finalize_callbacks(integrator)
    return TimeIntegratorSolution((first(prob.tspan), integrator.t),
                                  (prob.u0, integrator.u),
                                  integrator.sol.prob)
end

function step!(integrator::SimpleIntegrator2Nstar)
    @unpack prob = integrator.sol
    @unpack alg = integrator
    t_end = last(prob.tspan)
    callbacks = integrator.opts.callback

    @assert !integrator.finalstep
    if isnan(integrator.dt)
        error("time step size `dt` is NaN")
    end

    # if the next iteration would push the simulation beyond the end time, set dt accordingly
    if integrator.t + integrator.dt > t_end ||
       isapprox(integrator.t + integrator.dt, t_end)
        integrator.dt = t_end - integrator.t
        terminate!(integrator)
    end

    @trixi_timeit timer() "SSP Low-Storage Runge-Kutta step" begin
        # q2 = q1
        @threaded for i in eachindex(integrator.u)
            integrator.u_tmp[i] = integrator.u[i]
        end

        t_stage = integrator.t
        integrator.f(integrator.du_tmp, integrator.u, prob.p, t_stage) # f(q2) = f(q1)

        # q1 = q1 + gamma_21 * dt * f(q1) = q1 + gamma_21 * dt * f(q2)
        @threaded for i in eachindex(integrator.u)
            integrator.u[i] += alg.gamma[1] * integrator.du_tmp[i] * integrator.dt
        end

        for stage in 2:length(alg.c)
            t_stage = integrator.t + integrator.dt * alg.c[stage]
            integrator.f(integrator.du, integrator.u, prob.p, t_stage) # f(q1)

            # First column of lambda matrix
            lambda_stage = alg.lambda[2 * (stage - 1)]

            # First column of gamma matrix
            gamma_stage_1 = alg.gamma[2 * (stage - 1)]
            # Subdiagonal of gamma matrix
            gamma_stage_2 = alg.gamma[2 * (stage - 1) + 1]

            @threaded for i in eachindex(integrator.u)
                integrator.u[i] = lambda_stage * integrator.u_tmp[i] +
                                  (1 - lambda_stage) * integrator.u[i] +
                                  integrator.dt *
                                  (gamma_stage_1 * integrator.du_tmp[i] +
                                   gamma_stage_2 * integrator.du[i])
            end
        end
    end

    integrator.iter += 1
    integrator.t += integrator.dt

    # handle callbacks
    if callbacks isa CallbackSet
        foreach(callbacks.discrete_callbacks) do cb
            if cb.condition(integrator.u, integrator.t, integrator)
                cb.affect!(integrator)
            end
            return nothing
        end
    end

    # respect maximum number of iterations
    if integrator.iter >= integrator.opts.maxiters && !integrator.finalstep
        @warn "Interrupted. Larger maxiters is needed."
        terminate!(integrator)
    end
end

# get a cache where the RHS can be stored
get_du(integrator::SimpleIntegrator2Nstar) = integrator.du
get_tmp_cache(integrator::SimpleIntegrator2Nstar) = (integrator.u_tmp,
                                                     integrator.du_tmp)

# some algorithms from DiffEq like FSAL-ones need to be informed when a callback has modified u
u_modified!(integrator::SimpleIntegrator2Nstar, ::Bool) = false

# used by adaptive timestepping algorithms in DiffEq
function set_proposed_dt!(integrator::SimpleIntegrator2Nstar, dt)
    integrator.dt = dt
end

# Required e.g. for `glm_speed_callback` 
function get_proposed_dt(integrator::SimpleIntegrator2Nstar)
    return integrator.dt
end

# stop the time integration
function terminate!(integrator::SimpleIntegrator2Nstar)
    integrator.finalstep = true
    empty!(integrator.opts.tstops)
end

# used for AMR
function Base.resize!(integrator::SimpleIntegrator2Nstar, new_size)
    resize!(integrator.u, new_size)
    resize!(integrator.du, new_size)
    resize!(integrator.u_tmp, new_size)
    resize!(integrator.du_tmp, new_size)
end
end # @muladd
