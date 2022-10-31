using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase

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