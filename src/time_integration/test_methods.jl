### Finding the maximal stable CFL number for a chosen elixir ###

# elixir = "../../examples/tree_1d_dgsem/elixir_euler_sedov_blast_wave.jl"  
elixir = "../../examples/tree_1d_dgsem/elixir_mhd_torrilhon_shock_tube.jl"  
# elixir = "../../examples/tree_2d_dgsem/elixir_euler_kelvin_helmholtz_instability_amr.jl"
# elixir = "../../examples/tree_2d_dgsem/elixir_mhdmulti_rotor.jl"
# elixir = "../../examples/tree_2d_dgsem/elixir_mhd_orszag_tang.jl"

global cfl_number = 3.0
println("Current CFL Number: ", cfl_number)
include(elixir)
println("--- ", analysis_callback(sol)[1][1])
sol_ref = sol
sol_l2err_ref = analysis_callback(sol_ref)[1][1]

eps = 10
sol_l2err_vector = []
push!(sol_l2err_vector, sol_l2err_ref)

while abs(maximum(sol_l2err_vector))-sol_l2err_ref <= eps
   global cfl_number += 0.1
   println("Current CFL Number: ", cfl_number)
   include(elixir)
   println("--- ", analysis_callback(sol)[1][1])
   push!(sol_l2err_vector, analysis_callback(sol)[1][1])
end

println("Max CFL Number: ", (cfl_number-0.1))
println("L2-Error: ")
display(sol_l2err_vector)
