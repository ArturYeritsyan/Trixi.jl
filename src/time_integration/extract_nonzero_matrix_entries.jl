### Writing the found low storage methods into a json file ### (not needed anymore)
using JSON
using DelimitedFiles

# read Butcher K, Lambda, Gamma from csv file
matrices = readdlm("src/time_integration/matrices.csv";)  # 3x (s+1)x(s+1) matrices (Butcher K, Lambda, Gamma in this order)
display(matrices)

NumStages = size(matrices)[2]-1
lambda = Vector{Float64}(undef, 2*NumStages-1)
gamma = Vector{Float64}(undef, 2*NumStages-1)
c = zeros(NumStages)

# read lambda, gamma
global k = 1
for i in 2:NumStages+1
   for j in 1:NumStages+1
      if j==1 || i==j+1
         lambda[k] = matrices[i+(NumStages+1),j]  # lines 8 to 14
         gamma[k] = matrices[i+2*(NumStages+1),j]  # lines 15 to 21
         global k += 1
      end
   end
end

display(lambda)
display(gamma)

# calculate c = sum(Aij) over j
for i in 1:NumStages
   for j in 1:NumStages
      c[i] += matrices[i,j]
   end
end
display(c)


# write lambda, gamma, c to json file
data_dict = Dict(
    string(NumStages) => Dict("Lambda" => lambda, "Gamma" => gamma, "c" => c)
)

current_jsonfile = JSON.parsefile("src/time_integration/SSP_2Nstar_methods.json")

merge!(current_jsonfile, data_dict)

open("src/time_integration/SSP_2Nstar_methods.json", "w") do f
   JSON.print(f, current_jsonfile)
end
