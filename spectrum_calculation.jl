##
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase, Optim
include("helperFuncs.jl")
##

# Step 0: Get a sense of Δ(Γ)
##
J = 1
nt = 1
xaxis = []
yaxis = []
for Γ in 0.0:0.01:1
    push!(xaxis, Γ)
    adj_mat = build_triangular_chain_adjacency_matrix(nt);
    x_part, z_part = form_static_ham_ops(adj_mat, J, Γ)
    h_static = z_part - x_part
    spec = eigvals(h_static)
    E0 = spec[1]
    j = 2
    E1 = spec[j]
    while isapprox(E1, E0)
        j += 1
        E1 = spec[j]
    end
    Δ = E1 - E0
    push!(yaxis, Δ)
end
plot(xaxis, yaxis, label="nt=$(nt)")
title!("Triangular plaquette")
xlabel!("Γ")
ylabel!("Δ")
##

# Step 1: Find Hamiltonian J and Γ values with
# smallest spectrum by minimization
##
function compute_delta(Γ, nt=1)
    """
    Computes Δ = E₁-E₀, difference
    in energies between ground
    sub-space and 1st excited sub-space.
    Inputs
    Γ -- float, trasverse field strength
    nt -- number triangular plaquettes
    """
    adj_mat = build_triangular_chain_adjacency_matrix(nt);
    x_part, z_part = form_static_ham_ops(adj_mat, J, Γ)
    h_static = z_part - x_part
    spec = eigvals(h_static)
    E0 = spec[1]
    j = 2
    E1 = spec[j]
    while isapprox(E1, E0)
        j += 1
        E1 = spec[j]
    end
    Δ = E1 - E0;

    return Δ
end
##

##
for Γ=0:0.1:2
    adj_mat = build_triangular_chain_adjacency_matrix(nt);
    x_part, z_part = form_static_ham_ops(adj_mat, J, Γ)
    h_static = z_part - x_part
    E = eigvals(h_static)
    push!(energy_list, E)
end
##

##
nt = 1
adj_mat = build_triangular_chain_adjacency_matrix(nt);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ)
h_static = z_part - x_part
H = DenseHamiltonian([(s)-> s, (s)->1-s], [z_part, x_part])
##

##
Γ = 0.1
adj_mat = build_triangular_chain_adjacency_matrix(1);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ)
h_static = z_part - x_part
E = eigvals(h_static)
##

## perform optimizations
##
nt1_optim = optimize(x -> compute_delta(x, 1), 0.0, 1.0)
nt2_optim = optimize(x -> compute_delta(x, 2), 0.0, 1.0)
nt3_optim = optimize(x -> compute_delta(x, 3), 0.0, 1.0)
##

## step 2a: make combined plot data
##
xaxis=0.0:0.01:1
yaxis = [[], [], []]
for nt=1:3
    for Γ in xaxis
        push!(yaxis[nt], compute_delta(Γ, nt))
    end
end
plot(xaxis, yaxis, label=["nₜ=$(k)" for k=(1:3)'])
xlabel!("Γ")
ylabel!("Δ")
title!("Triangular plaquette")
#for k=2:3
#    plot!(xaxis, yaxis[k], label="nₜ = $(k)")
#end
##

## step 2b: make combined plot


##
function find_min_spec(Γmin, Γmax, Δi, iter, Γprev)
    """
    Finds minimum spectrum of h_static by 
    bisection search.
    """
    iter += 1
    Γ = (Γmax + Γmin) / 2
    Δ_iplus1 = compute_delta(Γ)
    if iter > 100
        if Δ_iplus1 < Δi
            return Γ, Δ_iplus1
        else
            return Γprev, Δi
        end
    elseif Δ_iplus1 < Δi
        Γmax = Γ
        find_min_spec(Γmin, Γmax, Δ_iplus1, iter, Γ)
    else
        Γmin = Γ
        find_min_spec(Γmin, Γmax, Δ_iplus1, iter, Γ)
    end
end

Γ, Δ = find_min_spec(0, 1, 100, 0, 100)
##