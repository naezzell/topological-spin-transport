include("../buildham/build_benchmark_ham.jl")
using BenchmarkTools

# fix value of tf and solve for different values of n
t_vals = [1, 10, 100, 1000];
#t_vals = [1];
#n_vals = [2];
n_vals = [2, 3, 4, 5, 6, 7, 8, 9];
ηvals = [1e-4, 1e-3, 1e-2, 1e-1];
fcvals = [2*π*4];
Tvals = [12];
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5;

function prep_and_solve_traj_ame(n, tf, η, fc, T)
    annealing = form_traj_ame_annealing_obj(n, tf, 10, η, fc, T)
    prob = build_ensembles(annealing, tf, :ame)
    t_list = range(0,tf,length=200)
    sol = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, reltol=1e-6, saveat=t_list)

    return sol
end

open("traj_ame_benchmark.txt", "w") do f
    write(f, "n,tf,η,fc,T,median_run_time (ns),cumulative_memory (bytes),cumulative_allocs,num_times\n")
end

for T in Tvals
    for fc in fcvals
        for η in ηvals
            for tf in t_vals
                for n in n_vals
                    bresult = @benchmark prep_and_solve_traj_ame($n, $tf, $η, $fc, $T)
                    exe_times = bresult.times;
                    med_time = median(exe_times);
                    num_times = length(exe_times);
                    memory = bresult.memory;
                    allocs = bresult.allocs;
                    open("traj_ame_benchmark.txt", "a") do f
                        write(f, "$(n),$(tf),$(η),$(fc),$(T),$(med_time),$(memory),$(allocs),$(num_times)\n")
                    end
                end
            end
        end
    end
end