##
using OpenQuantumTools, OpenQuantumBase
include("helperFuncs.jl")
using Statistics
using Combinatorics
##

############################################# 
# Step 1: Understanding single triangle
#############################################
##
## Step 1a: Purely classical frustrated triangle
##
J = 1
nt = 1
Γ = 0;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - Γ * x_part;
# form projector onto ground-state by diagolizing
vals, vecs = eigen(h0);
gs_states = [vecs[:,i] for i=1:6];
gs_proj = sum([gs_states[i] * adjoint(gs_states[i]) for i=1:6]);
# form expected groundstate projector
z = q_translate_state("0");
z3 = q_translate_state("000");
o = q_translate_state("1");
o3 = q_translate_state("111");
excited_proj = z3 * adjoint(z3) + o3 * adjoint(o3);
guess_gs_proj = σi ⊗ σi ⊗ σi - excited_proj;

# compare two approaches
norm(guess_gs_proj - gs_proj)
##

##
## Step 1b: Small Γ effect
##
J = 1
nt = 1
Γ = 0.1;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - Γ * x_part;
# form projector onto ground-state and first excited state subspace
vals, vecs = eigen(h0);
gs = vecs[:, 1];
gs_proj = sum([gs_states[i] * adjoint(gs_states[i]) for i=1:1]);
ex1_states = [vecs[:, k] for k=2:3];
ex1_proj = sum([ex1_states[i] * adjoint(ex1_states[i]) for i=1:length(ex1_states)]);
ts = [vecs[:, k] for k=1:6];
pert_proj = sum([ts[i] * adjoint(ts[i]) for i=1:length(ts)]);

# form predicted projectors
pred_gs_proj = 0;
plus = q_translate_state("0+1", normal=:true);
sc = unique(Combinatorics.permutations([o, z, plus]));
temp_s = [⊗(sc[i][1], sc[i][2], sc[i][3]) for i=1:length(sc)]
energies = [adjoint(temp_s[i]) * h0 * temp_s[i] for i=1:length(sc)]
pred_ex1_proj = sum([temp_s[i] * adjoint(temp_s[i]) for i=1:length(temp_s)]);
norm(pred_ex1_proj - ex1_proj)
norm(pert_proj - pred_ex1_proj)
##

############################################# 
# Step 2: Understanding two triangles
#############################################
##
## Step 2a: Purely classical frustrated triangle
##
J = 1
nt = 2
Γ = 0;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - Γ * x_part;
# form projector onto ground-state by diagolizing
vals, vecs = eigen(h0);
gs_dicts = [vec_to_dict(vecs[:, k]) for k=1:2]
ex1_dicts = [vec_to_dict(vecs[:, k]) for k=3:10]
##

##
## Step 2b: Lightly perturbed 
##
J = 1
nt = 2
Γ = 0.01;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - Γ * x_part;
# form projector onto ground-state by diagolizing
vals, vecs = eigen(h0);
gs_dicts = [vec_to_dict(vecs[:, k]) for k=1:2]
ex1_dicts = [vec_to_dict(vecs[:, k]) for k=3:10]
##

##
## Step 2c: Largely perturbed 
##
J = 1
nt = 2
Γ = 1.25;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - Γ * x_part;
# form projector onto ground-state by diagolizing
vals, vecs = eigen(h0);
gs_dicts = [vec_to_dict(vecs[:, k]) for k=1:1]
ex1_dicts = [vec_to_dict(vecs[:, k]) for k=2:2]
##

############################################# 
# Step 3: Understanding three triangles
#############################################
##
## Step 3a: Purely classical frustrated triangle
##
J = 1
nt = 3
Γ = 0;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - Γ * x_part;
# form projector onto ground-state by diagolizing
vals, vecs = eigen(h0);
gs_dicts = [vec_to_dict(vecs[:, k]) for k=1:8]
ex1_dicts = [vec_to_dict(vecs[:, k]) for k=8:18]
##

## Step 3b: Lightly pertrubed
##
J = 1
nt = 3
Γ = 0.01;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - Γ * x_part;
# form projector onto ground-state by diagolizing
vals, vecs = eigen(h0);
gs_dicts = [vec_to_dict(vecs[:, k]) for k=1:8]
ex1_dicts = [vec_to_dict(vecs[:, k]) for k=8:18]
##

## Step 3b: Largely pertrubed
##
J = 1
nt = 3
Γ = 1.25;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - Γ * x_part;
# form projector onto ground-state by diagolizing
vals, vecs = eigen(h0);
gs_dicts = [vec_to_dict(vecs[:, k]) for k=1:1]
ex1_dicts = [vec_to_dict(vecs[:, k]) for k=2:2]
##

############################################# 
# Step 3: Distribution of tiling energies
#############################################
##
nt = 4
Γ = 1.25;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
h0 = z_part - x_part;
Egs = eigmin(h0);
tilings = enumerate_triangular_clock_tilings(nt);
tiling_energies = [real(adjoint(tilings[i]) * h0 * tilings[i]) for i=1:6];
print("$(Egs), $(tiling_energies)\n")
##