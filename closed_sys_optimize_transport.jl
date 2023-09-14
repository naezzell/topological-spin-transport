##
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase
using Optim, Suppressor
using QuantumControl
using QuantumPropagators.Controls: substitute
include("helperFuncs.jl")
const PROJECTDIR = dirname(Base.active_project())
projectdir(names...) = joinpath(PROJECTDIR, names...)
datadir(names...) = projectdir("data", names...)
##

##
J = 1
nt = 1
Γ = 0.1;
tf = 100;
tstops = range(0, tf, length=100);
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - x_part;
#z1_term = σz ⊗ σi ⊗ σi ⊗ σi ⊗ σi;
z1_term = σz ⊗ σi ⊗ σi;
u = PauliVec[3][1];
d = PauliVec[3][2];
r = (u + d) / sqrt(2);
l = (u - d) / sqrt(2);
#init_state = (u ⊗ u ⊗ r)
init_state = (u ⊗ u ⊗ r)

#init_state = (u ⊗ u ⊗ r ⊗ d ⊗ u)
#init_state = (d ⊗ u ⊗ r)
#init_state = (u ⊗ d ⊗ r ⊗ u ⊗ d + d ⊗ u ⊗ r ⊗ d ⊗ u) / sqrt(2);
#init_state = (d ⊗ d ⊗ r ⊗ u ⊗ d - u ⊗ u ⊗ r ⊗ u ⊗ d) / sqrt(2) ;
print("init spin: $(compute_all_pseudospin(init_state))\n")
print("init energy: $(real(adjoint(init_state) * h0 * init_state))\n")
#targ_state = d ⊗ u ⊗ u ⊗ u ⊗ u;
#print("targ spin: $(compute_all_pseudospin(targ_state))\n")
#print("targ energy: $(real(adjoint(targ_state) * h0 * targ_state))\n")
tilings = enumerate_triangular_clock_tilings(nt);
#init_state = tilings[1];
targ_state = tilings[1];
#desired_spin = compute_all_pseudospin(targ_state);

pert_h = h0 + z1_term;
##

##
function sin_cost(x)
    coeffs = [(s) -> 1, (s) -> -1, (s) -> (x[1] * sin(x[2] * tf * s))]
    ham_terms = [z_part, x_part, z1_term]
    H = DenseHamiltonian(coeffs, ham_terms, unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    #final_spins = compute_all_pseudospin(sol[end]);
    #diff = norm(final_spins - (1.5 + 0.8660254037844383im) * [1.0, 1.0, 1.0])
    diff = norm(sol[end] - targ_state)

    return diff
end

function sin_cost2(x)
    coeffs = [(s) -> 1, (s) -> -1, (s) -> (x[1] * sin(6.56 * tf * s))]
    ham_terms = [z_part, x_part, z1_term]
    H = DenseHamiltonian(coeffs, ham_terms, unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    final_spins = compute_all_pseudospin(sol[end]);
    diff = norm(final_spins - desired_spin);
    #diff = norm(final_spins - (1.5 + 0.8660254037844383im) * [1.0, 1.0, 1.0])
    #diff = norm(sol[end] - targ_state)

    return diff
end

function lin_cost(x)
    coeffs = [(s) -> 1, (s) -> -1, (s) -> x[1] * s * tf]
    ham_terms = [z_part, x_part, z1_term]
    H = DenseHamiltonian(coeffs, ham_terms, unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    #final_spins = compute_all_pseudospin(sol[end]);
    #diff = norm(final_spins - (1.5 + 0.8660254037844383im) * [1.0, 1.0, 1.0])
    diff = norm(sol[end] - targ_state)

    return diff
end
##

##
result = Optim.optimize(lin_cost, rand(1), NelderMead(),
Optim.Options(iterations=100))
##

##
result = Optim.optimize(sin_cost, rand(2), NelderMead(),
Optim.Options(iterations=100))
##

##
function get_sin_spin_dynamics(opt_x)
    coeffs = [(s) -> 1, (s) -> -1, (s) -> opt_x[1] * sin(opt_x[2] * tf * s)]
    ham_terms = [z_part, x_part, z1_term]
    H = DenseHamiltonian(coeffs, ham_terms, unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    times = sol.t
    states = sol.u
    func_vals = opt_x[1] * sin.(opt_x[2] * times)
    spins = [compute_all_pseudospin(states[i]) for i=1:length(times)]

    return times, spins, func_vals
end

function get_sin2_spin_dynamics(opt_x)
    coeffs = [(s) -> 1, (s) -> -1, (s) -> opt_x[1] * sin(6.56 * tf * s)]
    ham_terms = [z_part, x_part, z1_term]
    H = DenseHamiltonian(coeffs, ham_terms, unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    print(compute_fidelity(sol[end], targ_state))
    times = sol.t
    states = sol.u
    func_vals = opt_x[1] * sin.(6.56 * times)
    spins = [compute_all_pseudospin(states[i]) for i=1:length(times)]

    return times, spins, func_vals
end

function get_linear_spin_dynamics(m)
    coeffs = [(s) -> 1, (s) -> -1, (s) -> m * s * tf]
    ham_terms = [z_part, x_part, z1_term]
    H = DenseHamiltonian(coeffs, ham_terms, unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    times = sol.t
    states = sol.u
    func_vals = m * times
    spins = [compute_all_pseudospin(states[i]) for i=1:length(times)]

    return times, spins, func_vals
end

function get_linear_spin_dynamics(m, input_times)
    coeffs = [(s) -> 1, (s) -> -1, (s) -> m * s * tf]
    ham_terms = [z_part, x_part, z1_term]
    H = DenseHamiltonian(coeffs, ham_terms, unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf, Exprb32(), solve_times);
    times = sol.t
    states = sol.u
    func_vals = m * times
    spins = [compute_all_pseudospin(states[i]) for i=1:length(times)]

    return times, spins, func_vals
end

function get_constant_spin_dynamics(c)
    coeffs = [(s) -> 1, (s) -> -1, (s) -> c]
    ham_terms = [z_part, x_part, z1_term]
    H = DenseHamiltonian(coeffs, ham_terms, unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    times = sol.t
    states = sol.u
    func_vals = [c for i=1:length(times)]
    spins = [compute_all_pseudospin(states[i]) for i=1:length(times)]

    return times, spins, func_vals
end

function get_constant_spin_dynamics(c, solve_times)
    coeffs = [(s) -> 1, (s) -> -1, (s) -> c]
    ham_terms = [z_part, x_part, z1_term]
    H = DenseHamiltonian(coeffs, ham_terms, unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf, Exprb32(), solve_times);
    times = sol.t
    states = sol.u
    func_vals = [c for i=1:length(times)]
    spins = [compute_all_pseudospin(states[i]) for i=1:length(times)]

    return times, spins, func_vals
end
##

# set up as quantum optimal control problem
##
ϵ(t) = QuantumControl.Shapes.flattop(t, T=tf, t_rise=0.3, func=:blackman);
"""Two-level-system Hamiltonian."""
function tls_hamiltonian(ϵ=ϵ)
    H0 = z_part - x_part
    H1 = z1_term
    return hamiltonian(H0, (H1, ϵ))
end;
H = tls_hamiltonian();
tlist = collect(range(0, tf, length=500));
function plot_control(pulse::Vector, tlist)
    plot(tlist, pulse, xlabel="time", ylabel="amplitude", legend=false)
end

plot_control(ϵ::Function, tlist) = plot_control([ϵ(t) for t in tlist], tlist);

#init_state = q_translate_state("10000", normal=true);
#targ_state = q_translate_state("00001", normal=true);

objectives = [Objective(initial_state=init_state, generator=H, target_state=targ_state)]
##

##
guess_dynamics = propagate_objective(
    objectives[1],
    problem.tlist;
    storage=true,
    #observables=(Ψ -> abs.(Ψ) .^ 2,)
)
times = problem.tlist;
spins = [compute_all_pseudospin(guess_dynamics[:, k]) for k=1:500];
fvals = [0 for k=1:500];

compute_fidelity(guess_dynamics[:,end], targ_state)
##

##
problem = ControlProblem(
    objectives=objectives,
    lambda_a=1/10,
    update_shape=(t -> QuantumControl.Shapes.flattop(t, T=tf, t_rise=0.3, func=:blackman)),
    tlist=tlist,
    iter_stop=1000,
    J_T=QuantumControl.Functionals.J_T_ss,
    check_convergence=res -> begin
        ((res.J_T < 1e-3) && (res.converged = true) && (res.message = "J_T < 10⁻³"))
    end
);
##


##
val = rand();
opt_result = @optimize_or_load(
    datadir("TLSOCT#J_T=J_T_ss#iter_stop=1000#method=krotov$(val).jld2"),
    problem,
    method = :krotov,
);
##

##
opt_dynamics = propagate_objective(
    substitute(objectives[1], IdDict(ϵ => opt_result.optimized_controls[1])),
    problem.tlist;
    storage=true
)
print(compute_fidelity(opt_dynamics[:,end], targ_state))
fig = plot_control(opt_result.optimized_controls[1], tlist)
##

##
times = problem.tlist
fvals = opt_result.optimized_controls[1]
spins = [compute_all_pseudospin(opt_dynamics[:,i]) for i=1:length(times)]
##

##
guess_dynamics = propagate_objective(
    objectives[1],
    problem.tlist;
    storage=true,
    #observables=(Ψ -> abs.(Ψ) .^ 2,)
)
times = problem.tlist;
spins = [compute_all_pseudospin(guess_dynamics[:, k]) for k=1:500];
fvals = [0 for k=1:500];

compute_fidelity(guess_dynamics[:,end], targ_state)
##

# Extract dominant frequency by FFT
##
# center the data to remove DC component
c_data = fvals .- mean(fvals)
# get magnitudes of fft
y = abs.(rfft(c_data));
# get frequencies
freqs = fftfreq(length(y), 1 / (times[2] - times[1]))
max_peak = findmax(y);
peaked_freq = freqs[max_peak[2]];
##