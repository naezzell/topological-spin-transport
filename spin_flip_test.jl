##
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase
include("helperFuncs.jl")
##

##
J = 1
nt = 1
Γ = 0.1
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
x_drive = σx ⊗ σi ⊗ σi;
amp = 100;
function driving(amp)
    tf = π / (2 * amp);
    tvals = LinRange(0, tf, 1000)
    func = amp * interval.(tvals, 0, tf)
    slist = tvals / tf;
    i_func = construct_interpolations(slist, func)

    return i_func, tf
end
i_func, tf = driving(amp);
#H = DenseHamiltonian([(s) -> 1, (s) -> i_func(s)], [z_part, x_drive], unit=:ħ)
H = DenseHamiltonian([(s) -> i_func(s)], [x_drive], unit=:h)
##

## 
init_state = q_translate_state("011");
sol = solve_closed_annealing_obj(H, init_state, tf);
##

##
desired = q_translate_state("111");
compute_fidelity(desired, sol[end])
##

##
compute_psuedo_spin(sol[end], 1, 2, 3)
##

# bare bones thing
##
init_state = q_translate_state("0");
H = DenseHamiltonian([(s) -> (π/2)], [σx])
sol = solve_closed_annealing_obj(H, init_state, 1);
desired = q_translate_state("1");
fid = norm(adjoint(desired) * sol[end])^2
##

##
init_state = q_translate_state("0");
H = DenseHamiltonian([(s) -> (π/2)], [σx])
sol = solve_closed_annealing_obj(H, init_state, 1 / (2 * π));
desired = q_translate_state("1");
fid = norm(adjoint(desired) * sol[end])^2
##

##
init_state = q_translate_state("0");
H = DenseHamiltonian([(s) -> (π/2)], [σx], unit=:ħ)
sol = solve_closed_annealing_obj(H, init_state, 1);
desired = q_translate_state("1");
fid = norm(adjoint(desired) * sol[end])^2
##

##
init_state = q_translate_state("0");
amp = 2.2346;
tf = (π / 2)  / amp;
H = DenseHamiltonian([(s) -> amp], [σx], unit=:ħ)
sol = solve_closed_annealing_obj(H, init_state, tf);
desired = q_translate_state("1");
fid = norm(adjoint(desired) * sol[end])^2
##


##
# Adding a latex σz term
##
##
init_state = q_translate_state("0");
desired = q_translate_state("1");
amp_list = 0.1:1:100
fid_list = []
for amp in amp_list    
    tf = (π / 2)  / amp;
    H = DenseHamiltonian([(s) -> 1, (s) -> amp], [σz, σx], unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    fid = norm(adjoint(desired) * sol[end])^2
    push!(fid_list, fid)
end
##

##
# Trying to flip a spin under the influence of two qubits
##
##
init_state = q_translate_state("000");
desired = q_translate_state("100");
amp_list = 0.1:1:100
zzi = σz ⊗ σz ⊗ σi
ziz = σz ⊗ σi ⊗ σz
izz = σi ⊗ σz ⊗ σz
hz = zzi + ziz + izz
xi = σx ⊗ σi ⊗ σi
fid_list = []
for amp in amp_list
    tf = (π / 2)  / amp;
    H = DenseHamiltonian([(s) -> 1, (s) -> amp], [hz, xi], unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    fid = norm(adjoint(desired) * sol[end])^2
    push!(fid_list, fid)
end
##

##
# 
##
##
init_state = q_translate_state("000");
desired = q_translate_state("100");
amp_list = 0.1:1:100
zzi = σz ⊗ σz ⊗ σi
ziz = σz ⊗ σi ⊗ σz
izz = σi ⊗ σz ⊗ σz
hz = zzi + ziz + izz
xi = σx ⊗ σi ⊗ σi
fid_list = []
for amp in amp_list
    tf = (π / 2)  / amp;
    H = DenseHamiltonian([(s) -> 1, (s) -> amp], [hz, xi], unit=:ħ)
    sol = solve_closed_annealing_obj(H, init_state, tf);
    fid = norm(adjoint(desired) * sol[end])^2
    push!(fid_list, fid)
end
##