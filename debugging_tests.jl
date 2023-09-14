##
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase
include("helperFuncs.jl")
##

# testing that hamiltonian's are correct
##
J = -0.75;
nt = 1;
Γ = 0.1;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
# x_part should be Γ (X_1 + X_2 + X_3)
x_correct = Γ * (σx ⊗ σi ⊗ σi + σi ⊗ σx ⊗ σi + σi ⊗ σi ⊗ σx);
x_part == x_correct
# z_part should be J (σ
#z_correct = J * (σz ⊗ σz ⊗ σi + σi ⊗ σz ⊗ σz + σz ⊗ σi ⊗ σz);
#z_part == z_correct
##

##
nt = 1;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
J = 1;
Γ = 0.1;
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
H = DenseHamiltonian([(s)-> (1-s), (s) -> s], [-x_part, z_part], unit=:ħ)
plot(H, 0:0.01:1, 8, linewidth=2, legend=:bottom)
##

##
J = 1
nt = 1
Γ = 0
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
#H_hbar = DenseHamiltonian([(s)-> -1, (s) -> 1], [x_part, z_part], unit=:ħ)
H_h = DenseHamiltonian([(s)-> 1, (s) -> 1], [-x_part, z_part], unit=:h)
# evaluates H_h(0) in units given, h
mat_h = evaluate(H_h, 0)
# evaluates H_h(0) in units of ħ
mat_ħ = H_h(0)
##

##
J = 1
nt = 1
Γ = 0.1
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
#H_hbar = DenseHamiltonian([(s)-> -1, (s) -> 1], [x_part, z_part], unit=:ħ)
H_h = DenseHamiltonian([(s)-> -1, (s) -> 1], [x_part, z_part], unit=:h)
# evaluates H_h(0) in units given, h
mat_h = evaluate(H_h, 0)
# compute evals
evals = eigvals(mat_h)
print(evals)
##

##
J = 1
nt = 1
Γ = 1
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
H = form_dynamic_hamiltonian(x_part, z_part, func);
##

##
nt = 1;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
J = 1;
Γ = 0.1;
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
H = DenseHamiltonian([(s)-> s, (s) -> (1-s)], [-x_part, z_part], unit=:ħ)
# get specrum with HOQST at s 
hoqst_evals, hoqst_evecs = eigen_decomp(H, 1; lvl=8)
hoqst_elist = [(hoqst_evals[j], hoqst_evecs[:, j]) for j=1:8]
# get manual spectrum at s 
mat_H = evaluate(H, 1)
julia_evals, julia_evecs = eigen(mat_H)
julia_elist = [(julia_evals[j], julia_evecs[:, j]) for j=1:8]

# compute difference between spectrums
one_norm_spec = sum(abs.(hoqst_evals - julia_evals))
# compute fidelities between eigenvectors... not the right thing
fids = [compute_fidelity(hoqst_elist[j][2], julia_elist[j][2]) for j=1:8]

hoqst_faithful = sum([norm(mat_H * hoqst_elist[k][2] - hoqst_elist[k][1] * hoqst_elist[k][2]) for k=1:8])
julia_faithful = sum([norm(mat_H * julia_elist[k][2] - julia_elist[k][1] * julia_elist[k][2]) for k=1:8])
#cross_faithful = sum([norm(mat_H * hoqast_elist[k][2] - hoqst_elist[k][1] * hoqst_elist[k][2]) for k=1:8])
##

##
nt = 1;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
J = 1;
Γ = 0.1;
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
H = DenseHamiltonian([(s)-> s, (s) -> (1-s)], [-x_part, z_part], unit=:ħ)
# get specrum with HOQST at s = 0
hoqst_evals, hoqst_evecs = eigen_decomp(H, 1; lvl=8)
hoqst_elist = [(hoqst_evals[j], hoqst_evecs[j, :]) for j=1:8]
# get manual spectrum at s = 0
mat_H = evaluate(H, 1)
julia_evals, julia_evecs = eigen(mat_H)
julia_elist = [(julia_evals[j], julia_evecs[j, :]) for j=1:8]

julia_elist == hoqst_elist
##

##
nt = 1;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
J = 1;
Γ = 0.1;
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
H = DenseHamiltonian([(s)-> s, (s) -> (1-s)], [-x_part, z_part], unit=:ħ)
# get specrum with HOQST at s = 0
hoqst_evals, hoqst_evecs = eigen_decomp(H, 0; lvl=8)
hoqst_elist = [(hoqst_evals[j], hoqst_evecs[j, :]) for j=1:8]
# get manual spectrum at s = 0
mat_H = H(0)
julia_evals, julia_evecs = eigen(mat_H)
julia_elist = [(julia_evals[j], julia_evecs[j, :]) for j=1:8]

julia_elist == hoqst_elist
##