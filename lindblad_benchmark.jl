include("../buildham/build_benchmark_ham.jl")
using BenchmarkTools

# fix value of tf and solve for different values of n
t_vals = [1, 10, 100, 1000];
n_vals = [2, 3, 4, 5, 6];
gamma_vals = [0.001, 0.01, 0.1];
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5;

function prep_and_solve_lindblad(n, tf, gamma)
    annealing = form_lindblad_annealing_obj(n, tf, gamma)
    sol = solve_lindblad(annealing, tf, alg=Tsit5());

    return sol
end

open("lindblad_benchmark.txt", "w") do f
    write(f, "n,tf,gamma,median_run_time (ns),cumulative_memory (bytes),cumulative_allocs,num_times\n")
end

for tf in t_vals
    for gamma in gamma_vals
        for n in n_vals
            bresult = @benchmark prep_and_solve_lindblad($n, $tf, $gamma);
            exe_times = bresult.times;
            med_time = median(exe_times);
            num_times = length(exe_times);
            memory = bresult.memory;
            allocs = bresult.allocs;
            open("lindblad_benchmark.txt", "a") do f
                write(f, "$(n),$(tf),$(gamma),$(med_time),$(memory),$(allocs),$(num_times)\n")
            end
        end
    end
end
