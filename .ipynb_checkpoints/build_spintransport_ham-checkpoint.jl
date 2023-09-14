using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase

function build_triangular_chain_dict(t::Int64)
   """
   Builds a chain of [t]
   triangular plaquettes.
   """
   n = t + 2
   adj_dict = Dict()
   for i=0:n-2
       if i < n - 2
           adj_dict[(i, i + 1)] = 1
           adj_dict[(i, i + 2)] = 1
       else
           adj_dict[(i, i + 1)] = 1
       end
   end

   return adj_dict
end

function hexagonal_grid_qubit_count(m::Int64, n::Int64)
    """
    Gets number of qubits in [n] x [m]
    grid of hexagons.
    """
    return (2 * m + 1) * (n + 2) - (m + 1)
end
    
function form_occupancy_matrix(m, n)
    """
    Forms matrix with 1 if a qubit is in position
    (i, j) and 0 if not.
    """
    occ_mat = ones((2 * m + 1, n + 2))
    for j=1:2:2 * m + 1
      occ_mat[j, n + 2] = 0
    end

    return occ_mat
end

function build_hexagonal_grid_dict(n::Int64, m::Int64)
    
    occ_mat = form_occupancy_matrix(n, m)
    r, c = size(occ_mat)
    
    adj_dict = Dict()
    map_from_ij_to_k = Dict()
    running_index = 0
    for i=1:r
       for j=1:c
           neighbors = Vector{Tuple{Int64, Int64}}([])
           # make sure a qubit is actually there
           if occ_mat[i, j] != 0
                map_from_ij_to_k[(i, j)] = running_index
                running_index += 1
                # last row is easy
                if i == r
                    if j < r
                        push!(neighbors, (i, j + 1))
                    end
                # for even rows before the end
                elseif (i + 1) % 2 == 0 && i < r
                    push!(neighbors, (i + 1, j))
                    push!(neighbors, (i + 1, j + 1))
                    if j < c - 1
                        push!(neighbors, (i, j + 1))
                    end
                              # for odd rows
                else
                    # last column is easy
                    if j == c
                        push!(neighbors, (i + 1, j - 1))
                    elseif j == 1
                        push!(neighbors, (i + 1, j))
                        push!(neighbors, (i, j + 1))
                    else
                        push!(neighbors, (i + 1, j - 1))
                        push!(neighbors, (i + 1, j))
                        push!(neighbors, (i, j + 1))
                    end
                end
           end
            if neighbors != []
                adj_dict[(i, j)] = neighbors
            end
       end
    end
    return adj_dict, map_from_ij_to_k
end

    """
    Given adjacency dictionary, builds the static
    portion of spin Hamiltonian,
    H = J \sum_{<i, j>} Z_i Z_j - \Gamma \sum_{i} X_i
    """

function form_static_ham(nq, adj_dict, qubit_map, J, Γ)
    op_list = []
    coeffs = []
    # add transverse field (X) if Γ != 0
    if !isapprox(Γ, 0)
        x_terms = standard_driver(nq)
        push!(op_list, x_terms)
        push!(coeffs, Γ)
    end
    # add ZZ terms always (if J = 0 what we doing?)
    zz_term_sum = nothing
    for (key, value) in adj_dict
        q0 = qubit_map[key]
        q1 = qubit_map[value[1]]
        zz_term = two_local_term(1, [q0, q1], nq)
        for v in value[2:end]
            zz_term += two_local_term(1, [q0, q1], nq)
        end
        if zz_term_sum == nothing 
            zz_term_sum = zz_term
        else
            zz_term_sum += zz_term
        end
    end 
    push!(op_list, zz_term_sum)
    push!(coeffs, J)

    return op_list, coeffs
end

        
function form_ham_ops(n)
    """
    Forms the three types of Hamiltonian
    terms, 
    H = X_i + Z_i + Z_i Z_j
    """
    # form base matrix operators
    x_terms = standard_driver(n)
    q_idx = 1:n
    coeffs = [1 for _=q_idx]
    z_terms = local_field_term(coeffs, q_idx, n)
    qq_idx = [[j,j+1] for j=1:n-1]
    coeffs = [1 for _=qq_idx]
    zz_terms = two_local_term(coeffs, qq_idx, n)
    op_list = [x_terms, z_terms, zz_terms]

    return op_list
end

function form_op_coefficients(tf)
    """
    Form the operator coefficients as a function
    of s, i.e. f(s(t)) and g(s(t)).
    """
    # evaluate functions at 10000 grid points
    tlist = range(0, tf; length=10000)
    ft = cos.(10 * tlist / tf)
    gt = sin.(10 * tlist / tf)
    # compute interpolations of functions h(s_k) = f(t_k)
    slist = tlist / tf
    hfs = construct_interpolations(slist, ft)
    hgs = construct_interpolations(slist, gt)

    return hfs, hgs
end

function form_annealing_ham(n, tf)
    """
    Forms annealing Hamiltonian object,
    H(s) = f(s) Z_i + g(s) X_i - Z_i Z_j. 
    """
    ham_ops = form_ham_ops(n)
    hfs, hgs = form_op_coefficients(tf)
    op_coeffs = [hfs, hgs, (s) -> -1.0]
    H = DenseHamiltonian(op_coeffs, ham_ops, unit=:ħ)

    return H
end

function form_initial_state(n)
    """
    Forms initial ground-state of the Hamiltonian,
    |1>^{⊗n}
    """
    u0 = q_translate_state("1"^n)

    return u0
end

function form_initial_densityop(n)
    """
    Forms initial ground-state of the Hamiltonian,
    |1>^{⊗n}
    """
    u0 = q_translate_state("1"^n)
    ρ0 = u0 * u0'

    return ρ0
end

function form_closed_annealing_obj(n, tf)
    """
    Forms an annealing object for the benchmark
    Hamiltonian for a closed system (Schrodinger)
    solver.
    """
    H = form_annealing_ham(n, tf)
    u0 = form_initial_state(n)
    annealing = Annealing(H, u0)

    return annealing
end

function form_lindblad_annealing_obj(n, tf, γ)
    """
    Forms annealing object for Lindblad equation
    evolution of density matrix.
    """
    H = form_annealing_ham(n, tf)
    ρ0 = form_initial_densityop(n)
    # form lindblad set
    z1 = local_field_term(1, 1, n)
    lind_list = [Lindblad(γ, z1)]
    for j=2:n
        zj = local_field_term(1, j, n)
        lj = Lindblad(γ, zj)
        push!(lind_list, lj)
    end
    int_set = InteractionSet(lind_list...)
    lind_annealing = Annealing(H, ρ0, interactions = int_set)

    return lind_annealing
end

function form_traj_ame_annealing_obj(n, tf, num_fluc, η, fc, T)
    """
    Forms trajectory AME Lindblad equation Annealing
    object ready for solving.
    """
    # form H and u0
    H = form_annealing_ham(n, tf)
    u0 = form_initial_state(n)
    # build fluctuators params (AME with spin fluctuators tutorial)
    bvec = 0.01 * ones(num_fluc)
    γvec = log_uniform(0.01, 1, num_fluc)
    ## Form iid Ohmic Dephasing coupling operators as interaction set
    # also add iid telegraph noise to each term
    #η=1e-4; fc=2*π*4; T=12;
    ohmic_bath = Ohmic(η, fc, T)
    int_list = []
    for j=1:n
        str = ""
        for k=1:j-1
            str *= "I"
        end
        str *= "Z"
        for k=j+1:n
            str *= "I"
        end
        z_coupling = ConstantCouplings([str], unit=:ħ)
        dephasing_int = Interaction(z_coupling, ohmic_bath)
        push!(int_list, dephasing_int)
        # add iid telegraph noise term
        fluctuator_ensemble = EnsembleFluctuator(bvec, γvec);
        fluc_int = Interaction(z_coupling, fluctuator_ensemble)
        push!(int_list, fluc_int)

    end
    int_set = InteractionSet(int_list...)
    annealing = Annealing(H, u0, interactions=int_set)

    return annealing
end