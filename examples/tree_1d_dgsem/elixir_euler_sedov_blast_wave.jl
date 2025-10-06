using OrdinaryDiffEqLowStorageRK
using Trixi

###############################################################################
# semidiscretization of the compressible Euler equations
equations = CompressibleEulerEquations1D(1.4)

"""
    initial_condition_sedov_blast_wave(x, t, equations::CompressibleEulerEquations1D)

The Sedov blast wave setup based on Flash
- https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node187.html#SECTION010114000000000000000
"""
function initial_condition_sedov_blast_wave(x, t, equations::CompressibleEulerEquations1D)
    # Set up polar coordinates
    RealT = eltype(x)
    inicenter = SVector(0)
    x_norm = x[1] - inicenter[1]
    r = abs(x_norm)

    # Setup based on https://flash.rochester.edu/site/flashcode/user_support/flash_ug_devel/node187.html#SECTION010114000000000000000
    r0 = 0.21875f0 # = 3.5 * smallest dx (for domain length=4 and max-ref=6)
    # r0 = 0.5 # = more reasonable setup
    E = 1
    p0_inner = 6 * (equations.gamma - 1) * E / (3 * convert(RealT, pi) * r0)
    p0_outer = convert(RealT, 1.0e-5) # = true Sedov setup
    # p0_outer = 1.0e-3 # = more reasonable setup

    # Calculate primitive variables
    rho = 1
    v1 = 0
    p = r > r0 ? p0_outer : p0_inner

    return prim2cons(SVector(rho, v1, p), equations)
end
initial_condition = initial_condition_sedov_blast_wave

# Up to version 0.13.0, `max_abs_speed_naive` was used as the default wave speed estimate of
# `const flux_lax_friedrichs = FluxLaxFriedrichs(), i.e., `FluxLaxFriedrichs(max_abs_speed = max_abs_speed_naive)`.
# In the `StepsizeCallback`, though, the less diffusive `max_abs_speeds` is employed which is consistent with `max_abs_speed`.
# Thus, we exchanged in PR#2458 the default wave speed used in the LLF flux to `max_abs_speed`.
# To ensure that every example still runs we specify explicitly `FluxLaxFriedrichs(max_abs_speed_naive)`.
# We remark, however, that the now default `max_abs_speed` is in general recommended due to compliance with the 
# `StepsizeCallback` (CFL-Condition) and less diffusion.
surface_flux = FluxLaxFriedrichs(max_abs_speed_naive)
volume_flux = flux_chandrashekar
basis = LobattoLegendreBasis(3)
shock_indicator_variable = density_pressure
indicator_sc = IndicatorHennemannGassner(equations, basis,
                                         alpha_max = 0.5,
                                         alpha_min = 0.001,
                                         alpha_smooth = true,
                                         variable = shock_indicator_variable)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

coordinates_min = (-2.0,)
coordinates_max = (2.0,)
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 6,
                n_cells_max = 10_000)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 12.5)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     analysis_errors = Symbol[],
                                     analysis_integrals = ())

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(interval = 100,
                                     save_initial_solution = true,
                                     save_final_solution = true,
                                     solution_variables = cons2prim)

amr_indicator = IndicatorHennemannGassner(semi,
                                          alpha_max = 0.5,
                                          alpha_min = 0.001,
                                          alpha_smooth = true,
                                          variable = density_pressure)
amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level = 4,
                                      max_level = 6, max_threshold = 0.01)
amr_callback = AMRCallback(semi, amr_controller,
                           interval = 5,
                           adapt_initial_condition = true,
                           adapt_initial_condition_only_refine = true)

stepsize_callback = StepsizeCallback(cfl = 1.5)

callbacks = CallbackSet(summary_callback,
                        analysis_callback, alive_callback,
                        # save_solution,
                        amr_callback, stepsize_callback)

###############################################################################
# run the simulation

sol = Trixi.solve(ode, Trixi.SSP144_2Nstar();
# sol = solve(ode, SSPRK104();
            dt = stepsize_callback(ode), # solve needs some value here but it will be overwritten by the stepsize_callback
            ode_default_options()..., callback=callbacks);

println(analysis_callback(sol))

# Plots
using Plots
pd = PlotData1D(sol)

plot(getmesh(pd), label = "", ylim = (-0.3, 1.4))

plot!(pd["rho"], xlabel = "\$x\$",
      label = "\$\\rho\$",
      linewidth = 3, color = RGB(0, 84 / 256, 159 / 256),
      guidefont = font("Computer Modern", 16), tickfont = font("Computer Modern", 14),
      yticks = -0.2:0.4:1.4,
      xtick = ([-2, -1, -0.5, 0, 0.5, 1, 2], ["-2", "-1", "-0.5", "0", "0.5", "1", "2"]),
      legend = true)

plot!(pd["v1"],
      label = "\$v\$",
      linewidth = 3, color = RGB(246 / 256, 169 / 256, 0), legend = true)

plot!(pd["p"], title = "Euler Sedov Blast Wave",
      label = "\$p\$",
      linewidth = 3, color = RGB(70 / 256, 171 / 256, 39 / 256),
      titlefont = font("Computer Modern", 18),
      legendfont = font("Computer Modern", 16),
      legend = :left)

savefig("examples/sedov_blast_wave.pdf")
