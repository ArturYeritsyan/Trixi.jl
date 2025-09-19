using OrdinaryDiffEq
using Trixi
using Printf

# Solving nonlinear Bernoulli equation: u' = t*u^2 - u
# Analytical solution: u(t) = 1/(t+1-0.5*exp(t))
# Calculating Order of convergence

analytical_solution(t) = 1 / (t + 1 - 0.5*exp(t))

# Setup ODE
function f(du, u, p, t)
  du[1] = t * u[1]^2 - u[1]
end

u0 = [2.0]
tspan = (0, 1.0)
ode = ODEProblem(f, u0, tspan)

# Solve ODE for different time steps dt
dts = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625]
errors = []
for dt in dts
  sol = Trixi.solve(ode, Trixi.SSP144_2Nstar(); dt = dt, callback=nothing)
  u_numerical = sol.u[end][1]
  u_exact = analytical_solution(tspan[2])
  push!(errors, abs(u_numerical - u_exact))
end

# Print errors and convergence order for every time step dt
println("      dt    |     Error     |   Order")
println("----------------------------------------------")
@printf(" %10.5f | %13.16e |    -    \n", dts[1], errors[1])

for i in 1:(length(errors) - 1)
  order = log(errors[i] / errors[i+1]) / log(dts[i] / dts[i+1])
  @printf(" %10.5f | %13.16e |  %.4f \n", dts[i+1], errors[i+1], order)
end
println("----------------------------------------------")

