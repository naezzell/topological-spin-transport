include("../buildham/build_benchmark_ham.jl")
using BenchmarkTools

# fix value of tf and solve for different values of n
t_vals = [1, 10, 100, 1000];
n_vals = [2, 3, 4, 5, 6, 7, 8, 9];
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5;

function prep_and_solve_schrodinger(n, tf)
    annealing = form_closed_annealing_obj(n, tf)
    sol = solve_schrodinger(annealing, tf, alg=Tsit5());

    return sol
end

open("schrodinger_benchmark.txt", "w") do f
    write(f, "n,tf,median_run_time (ns),cumulative_memory (bytes),cumulative_allocs,num_times\n")
end

for tf in t_vals
    for n in n_vals
        bresult = @benchmark prep_and_solve_schrodinger($n, $tf);
        exe_times = bresult.times;
        med_time = median(exe_times);
        num_times = length(exe_times);
        memory = bresult.memory;
        allocs = bresult.allocs;
        open("schrodinger_benchmark.txt", "a") do f
            write(f, "$(n),$(tf),$(med_time),$(memory),$(allocs),$(num_times)\n")
        end
    end
end
