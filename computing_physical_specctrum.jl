#
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase, Optim
include("helperFuncs.jl")
#

##
nt = 1;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
J = 1;
Γ = 0.1;
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
##

##
nt = 1;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
J = 1;
Γ = 0.1;
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
H = DenseHamiltonian([(s)-> s, (s) -> (1-s)], [-x_part, z_part],  unit=:h)
plot(H, 0:0.01:1, 8, linewidth=2, legend=:left)
##


##
nt = 1;
adj_mat = build_triangular_chain_adjacency_matrix(nt);
J = 1;
##

##
Γ = 1;
function compute_spectrum(Γ)
    x_part, z_part = form_static_ham_ops(adj_mat, J, Γ);
    H = DenseHamiltonian([(s)-> s, (s) -> (1-s)], [-x_part, z_part], unit=:ħ)
    z_decomp = eigen_decomp(H, 0, lvl=8)
    Egs = minimum(z_decomp[1])
    degen = sum(z_decomp[1] .== Egs)
    println(degen)
    x_decomp = eigen_decomp(H, 1, lvl=8)
    Δ = x_decomp[1][degen + 1] - x_decomp[1][degen]

    return Δ
end
##