### Measuring run time of a elixir ### (obsolte, not been used)
# elixir = "../../examples/tree_1d_dgsem/elixir_euler_sedov_blast_wave.jl"  
elixir = "../../examples/tree_1d_dgsem/elixir_mhd_torrilhon_shock_tube.jl"  
# elixir = "../../examples/tree_2d_dgsem/elixir_euler_kelvin_helmholtz_instability_amr.jl"
# elixir = "../../examples/tree_2d_dgsem/elixir_mhdmulti_rotor.jl"
# elixir = "../../examples/tree_2d_dgsem/elixir_mhd_orszag_tang.jl"

using Statistics

global cfl_number = 1.3

elapsed_times = []
for i=1:11
   t1 = time()
   include(elixir)
   elapsed_time = time() - t1;
   if i!=1  # ignore first simulation
      push!(elapsed_times, elapsed_time)
   end
end

println("elapsed times: ", elapsed_times)
println("median time: ", median(elapsed_times))