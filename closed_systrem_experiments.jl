##
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase
using Optim
include("helperFuncs.jl")
##

# Form static portion of Hamiltonian
##
J = 1
nt = 3
Γ = 1.0
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
m = 0.01;
total_time = 10;
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
#init_state = get_initial_ground_state(H);
#sol = solve_closed_annealing_obj(H, gs, tf);
init_state = q_translate_state("10000", normal=true);
sol = solve_closed_annealing_obj(H, init_state, tf);
##

## plot order parameter as a function of time
##
psuedo_spin1_list = []
psuedo_spin2_list = []
psuedo_spin3_list = []
for i=1:length(sol)
    state = sol[i]
    if norm(state) != 1
        state = state / norm(state)
    end
    ps1 = compute_psuedo_spin(state, 1, 2, 3)
    push!(psuedo_spin1_list, ps1)
    ps2 = compute_psuedo_spin(state, 2, 3, 4)
    push!(psuedo_spin2_list, ps2)
    ps3 = compute_psuedo_spin(state, 3, 4, 5)
    push!(psuedo_spin3_list, ps3)
end

plot(sol.t, real.(psuedo_spin1_list), label="ψ1 real", xlabel="time (ns)", ylabel="psuedo-spin value", title="Driving with Γ=$(Γ)")
plot!(sol.t, imag.(psuedo_spin1_list), label="ψ1 imag")
plot!(sol.t, real.(psuedo_spin2_list), label="ψ2 real", line=(:dash))
plot!(sol.t, imag.(psuedo_spin2_list), label="ψ2 imag", line=(:dash), legend=:bottom)
plot!(sol.t, real.(psuedo_spin3_list), label="ψ3 real", line=(:dashdot))
plot!(sol.t, imag.(psuedo_spin3_list), label="ψ3 imag", line=(:dashdot))
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

##
function ps_spin_from_str(str)
    psi = compute_psuedo_spin(q_translate_state(str), 1, 2, 3)

    return psi
end
##

##
max_transfer = norm.(psuedo_spin3_list) - norm.(psuedo_spin1_list) - norm.(psuedo_spin2_list)
plot(sol.t, max_transfer, xlabel="time (ns)", ylabel="max transfer")
##

##
max_trans_state = sol[argmax(max_transfer)]
##

##
#ρ = max_trans_state * adjoint(max_trans_state)
ket0 = q_translate_state("0")
proj0 = ket0 * adjoint(ket0)
ket1 = q_translate_state("1")
proj1 = ket1 * adjoint(ket1)
id2 = σi
##

##
states = ["000", "001", "010", "011", "100", "101", "110", "111"]
projectors = [proj0 ⊗ proj0 ⊗ proj0,
              proj0 ⊗ proj0 ⊗ proj1,
              proj0 ⊗ proj1 ⊗ proj0,
              proj0 ⊗ proj1 ⊗ proj1,
              proj1 ⊗ proj0 ⊗ proj0,
              proj1 ⊗ proj0 ⊗ proj1,
              proj1 ⊗ proj1 ⊗ proj0,
              proj1 ⊗ proj1 ⊗ proj1]
##

##
ρ = opt_dynamics[:,end] * adjoint(opt_dynamics[:,end])
tri1 = [real(tr(ρ * (projectors[j] ⊗ id2 ⊗ id2))) for j=1:8]
tri2 = [real(tr(ρ * (id2 ⊗ projectors[j] ⊗ id2))) for j=1:8]
tri3 = [real(tr(ρ * (id2 ⊗ id2 ⊗ projectors[j]))) for j=1:8]
##

##
bar(states, tri1, title="triangle 1")
##

##
bar(states, tri2, title="triangle 2")
##

##
bar(states, tri3, title="triangle 3")
##

##
mn = vcat(tri1, tri2, tri3)
nam = repeat(states, inner=3)

groupedbar(nam, mn)
##