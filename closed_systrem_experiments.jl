##
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase
include("helperFuncs.jl")
##

# Form static portion of Hamiltonian
##
J = 1
nt = 2
#Γ = 1.49e-8
Γ = 0.1
#Γ = 1
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
##

# Build time-dependent functions
##
#pause_time = 0.5;
#period = 1 / (Γ);
#amp = Γ;
#cycles = 0.5;
#func, tf, svals = form_pause_into_comb_function(pause_time, period, amp, cycles);
m = 0.001;
total_time = 100;
func, tf, svals = linear_ramp(m, total_time)
plot(svals * tf, func(svals), label=false)
xlabel!("t")
ylabel!("h₁(t)")
title!("m = $(m)")
##
|
# solve Schrodinger equation with local perturbation
##
H = form_dynamic_hamiltonian(x_part, z_part, func);
gs = get_initial_ground_state(H);
sol = solve_closed_annealing_obj(H, gs, tf);
#init_state = q_translate_state("011", normal=true);
#sol = solve_closed_annealing_obj(H, init_state, tf);
##

## plot order parameter as a function of time
##
psuedo_spin1_list = []
psuedo_spin2_list = []
for i=1:length(sol)
    state = sol[i]
    ps1 = compute_psuedo_spin(state, 1, 2, 3)
    push!(psuedo_spin1_list, ps1)
   ps2 = compute_psuedo_spin(state, 2, 3, 4)
   push!(psuedo_spin2_list, ps2)
end

plot(sol.t, real.(psuedo_spin1_list), label="t1 real", xlabel="time (ns)", ylabel="psuedo-spin value", title="Driving with Γ=$(Γ)")
plot!(sol.t, imag.(psuedo_spin1_list), label="t1 imag")
plot!(sol.t, real.(psuedo_spin2_list), label="t2 real")
plot!(sol.t, imag.(psuedo_spin2_list), label="t2 imag")
##

## a test of normalization...
##
H = DenseHamiltonian([(s)->1.0], [-σz], unit=:ħ)
u0 = PauliVec[1][1]
tf = 10000
annealing = Annealing(H, u0)
sol_tsit = solve_schrodinger(annealing, tf, alg=Tsit5(), abstol=1e-6, reltol=1e-6);
#sol_trbdf = solve_schrodinger(annealing, tf, alg=TRBDF2(), abstol=1e-6, reltol=1e-6);
# LinearExponential is a fixed-stepsize method, the user needs to specify the time steps using the keyword argument `tstops`.
#sol_linexp = solve_schrodinger(annealing, tf, alg=LinearExponential(), abstol=1e-6, reltol=1e-6, tstops=range(0,tf,length=100));
# Even though Exprb is an adaptive method, it tends to skip a lot of middle points. So if you want an accurate solution in the middle,
# it is better to manually add more points.
#sol_exprb32 = solve_schrodinger(annealing, tf, alg=Exprb32(), tstops=range(0,tf,length=100));
##