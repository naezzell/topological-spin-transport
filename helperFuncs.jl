##
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase
##

## ======================================
# Misc like building bitstrings
## ======================================
##
function get_bitstrings(n)
    """
    Generates all bitstrings of length n.
    """
    x = [bitstring(0)[end-n+1:end]]
    for i = 1:2^n-1
       s = bitstring(i)
       push!(x, s[end-n+1:end])
    end
    
    return x
end

function generate_prob_plot(psi)
    """
    Generates a barchart with bitstring
    probabilities.
    """
    probs = real(psi .* conj(psi))
    d = length(probs)
    n = trunc(Int, log2(d))
    x = get_bitstrings(n)

    return bar(x, probs)
end
##

## ======================================
# Building Hamiltonian
## ======================================
##
function get_hexagonal_qubit_count(n, m)
    """
    Gets number of qubits in [n] x [m]
    grid of hexagons.
    """
    return (2 * n + 1) * (m + 2) - (n + 1)
end

function form_occupancy_matrix(n, m)
    """
    Forms matrix with 1 if a qubit is in position
    (i, j) and 0 if not.
    """
    occ_mat = ones(Int64, (2 * n + 1, m + 2))
    for j=1:2:2 * n + 1
      occ_mat[j, m + 2] = 0
    end

    return occ_mat
end

function build_hexagonal_adjacency_matrix(n, m)
    """
    Builds adjacency matrix encoding
    a hexagonal tiling made up of triangular plaquettes.
    """
    # get number of qubits in tiling to build adjacency size
    nq = get_hexagonal_qubit_count(n, m)
    adj_mat = zeros((nq, nq))
    # build occupancy matrix to make (i, j) positions make sense
    occ_mat = form_occupancy_matrix(n, m)
    num_rows, num_cols = size(occ_mat)
    # the type of connections are determined by the node type which
    # is determined by (i, j) position in matrix
    k = 1
    ab1 = m + 1
    ab2 = m + 2
    for i=1:num_rows
        for j=1:num_cols
            if occ_mat[i, j] != 0
                # redefine i and j to 0 based indexing conditions
                i0 = i - 1
                j0 = j - 1
                # upper left
                if i0 == 0 && j0 == 0
                    adj_mat[k, k + 1] = 1
                    adj_mat[k, k + ab1] = 1
                    adj_mat[k, k + ab2] = 1
                # upper right
                elseif i0 == 0 && j0 == (num_cols - 2)
                    adj_mat[k, k - 1] = 1
                    adj_mat[k, k + ab1] = 1
                    adj_mat[k, k + ab2] = 1
                # bottom left
                elseif i0 == (num_rows - 1) && j0 == 0
                    adj_mat[k, k - ab2] = 1
                    adj_mat[k, k - ab1] = 1
                    adj_mat[k, k + 1] = 1
                # bottom right
                elseif i0 == (num_rows - 1) && j0 == (num_cols - 2)
                    adj_mat[k, k - ab2] = 1
                    adj_mat[k, k - ab1] = 1
                    adj_mat[k, k - 1] = 1
                # top middle (aka C)
                elseif i0 == 0 && j0 > 0 && j0 < (num_cols - 2)
                    adj_mat[k, k - 1] = 1
                    adj_mat[k, k + 1] = 1
                    adj_mat[k, k + ab1] = 1
                    adj_mat[k, k + ab2] = 1
                # bottom middle (aka C')
                elseif i0 == (num_rows - 1) && 0 < j0 && j0 < (num_cols - 2)
                    adj_mat[k, k - 1] = 1
                    adj_mat[k, k + 1] = 1
                    adj_mat[k, k - ab2] = 1
                    adj_mat[k, k - ab1] = 1
                # left-most hexagonal points (aka B)
                elseif i0 % 2 == 1 && j0 == 0
                    adj_mat[k, k - ab1] = 1
                    adj_mat[k, k + 1] = 1
                    adj_mat[k, k + ab2] = 1
                # right-most hexagonal points (aka B')
                elseif i0 % 2 == 1 && j0 == (num_cols - 1)
                    #print(k)
                    adj_mat[k, k - ab2] = 1
                    adj_mat[k, k - 1] = 1
                    adj_mat[k, k + ab1] = 1
                # bottom left extremal hexagonal points (aka D)
                elseif 0 < i0 && i0 < (num_rows - 1) && i0 % 2 == 0 && j0 == 0
                    adj_mat[k, k - ab2] = 1
                    adj_mat[k, k - ab1] = 1
                    adj_mat[k, k + 1] = 1
                    adj_mat[k, k + ab1] = 1
                    adj_mat[k, k + ab2] = 1
                # bottom right extremal hexagonal points (aka D')
                elseif 0 < i0 && i0 < (num_rows - 1) && i0 % 2 == 0 && j0 == (num_cols - 2)
                    #print(k)
                    adj_mat[k, k - ab2] = 1
                    adj_mat[k, k - ab1] = 1
                    adj_mat[k, k - 1] = 1
                    adj_mat[k, k + ab1] = 1
                    adj_mat[k, k + ab2] = 1
                # bulk, center of hexagons (aka A)
                #elif 0 < i && i < 2 * n && 0 < j && j < m + 1:
                else
                    #print(f"else also triggered: {k}")
                    adj_mat[k, k - ab2] = 1
                    adj_mat[k, k - ab1] = 1
                    adj_mat[k, k - 1] = 1
                    adj_mat[k, k + 1] = 1
                    adj_mat[k, k + ab1] = 1
                    adj_mat[k, k + ab2] = 1
                end
                k += 1
            end
        end
    end
   
    return adj_mat, occ_mat
end

function form_static_ham_ops(adj_mat, J = 1, gamma = 0.1)
    """
    Form static porition of Hamiltonian.
    """
    n = size(adj_mat)[1]
    # form x part
    x_part = gamma * standard_driver(n)
    # form zz part
    couplings_list = []
    for i=1:n
        for j=i:n
            if adj_mat[i, j] == 1
                c = [i, j]
                push!(couplings_list, c)
            end
        end
    end
    coeffs = [J for _=1:length(couplings_list)]
    zz_part = two_local_term(coeffs, couplings_list, n)
    
    return x_part, zz_part
end

function build_triangular_chain_dict(t)
    """
    Builds a chain of [t]
    triangular plaquettes.
    """
    n = t + 2
    adj_dict = Dict()
    for i=1:n - 1
        if i < n - 1
            adj_dict[(i, i + 1)] = 1
            adj_dict[(i, i + 2)] = 1
        else
            adj_dict[(i, i + 1)] = 1
        end
    end
            
    return adj_dict
end

function build_triangular_chain_adjacency_matrix(t)
    """
    Builds adjacency matrix given an adjacency
    dictionary.
    """
    n = t + 2
    adj_mat = zeros(Int64, (n, n))
    adj_dict = build_triangular_chain_dict(t)
    for key in keys(adj_dict)
        adj_mat[key[1], key[2]] = 1
        adj_mat[key[2], key[1]] = 1
    end
        
    return adj_mat         
end

function form_dynamic_hamiltonian(x_part, z_part, func)
    """
    Forms a Hamiltonian of the form
    """
    n = log2(size(x_part)[1])
    z1_drive = local_field_term(1, 1, n)
    ham_ops = [x_part, z_part, z1_drive]
    op_coeffs = [(s) -> -1.0, (s) -> 1.0, (s) -> func(s)]
    H = DenseHamiltonian(op_coeffs, ham_ops, unit=:ħ)

    return H
end
##

## ======================================
# Amplitude function builders
## ======================================
##
function heaviside(t)
    0.5 * (sign(t) + 1)
 end

function interval(t, a, b)
    heaviside(t - a) - heaviside(t - b)
 end

 function form_pause_function(pause_time)
    """
    Form the operator coefficients as a function
    of s, i.e. f(s(t)) and g(s(t)).
    """
    # build pause 
    pause_times = range(0, pause_time; length=100)
    pause = zeros(100)
    # compute interpolations of functions h(s_k) = f(t_k)
    slist = pause_times / pause_time
    i_func = construct_interpolations(slist, pause)

    return i_func, pause_time, slist
end

 function form_pause_into_comb_function(pause_time, period, amp, cycles)
    """
    Form the operator coefficients as a function
    of s, i.e. f(s(t)) and g(s(t)).
    """
    # build pause 
    pause_times = range(0, pause_time; length=100)
    pause = zeros(100)
    # evaluate functions at 10000 grid points
    comb_time = period * cycles
    tf = pause_time + comb_time
    start_t = pause_time + period/1000
    comb_times = range(start_t, tf; length=1000)
    comb = sum([(-1)^j * amp * interval.(comb_times, j * period / 2 + start_t, (j + 1) * period / 2 + start_t) for j=0:2*cycles])
    # compute interpolations of functions h(s_k) = f(t_k)
    tlist = vcat(pause_times, comb_times)
    func = vcat(pause, comb)
    slist = tlist / tf
    i_func = construct_interpolations(slist, func)

    return i_func, tf, slist
end

function linear_ramp(slope, total_time, y0)
    """
    Form the operator coefficients as a function
    of s, i.e. f(s(t)) and g(s(t)).
    """
    # build linear equation
    x_pts = LinRange(0, total_time, 10000)
    y_pts = y0 .+ slope * x_pts
    # compute interpolations of functions h(s_k) = f(t_k)
    slist = x_pts / total_time
    i_func = construct_interpolations(slist, y_pts)

    return i_func, total_time, slist
end

function linear_ramp(slope, total_time)
    """
    Form the operator coefficients as a function
    of s, i.e. f(s(t)) and g(s(t)).
    """
    # build linear equation
    x_pts = LinRange(0, total_time, 10000)
    y_pts = slope * x_pts
    # compute interpolations of functions h(s_k) = f(t_k)
    slist = x_pts / total_time
    i_func = construct_interpolations(slist, y_pts)

    return i_func, total_time, slist
end
##

## ======================================
# Solve with various master equations
## ======================================
##
function get_initial_ground_state(H)
    gs_info = eigen_decomp(H, 0; lvl=1)
    gs_vec = vec(gs_info[2])

    return gs_vec
end

function compute_groundstate(H)
    """
    Finds the lowest energy eigenvector.
    Takes an equal superposition of degenerate
    ground-states.
    """
    evals, evecs = eigen(H)
    minE = evals[1]
    newE = evals[2]
    ground_state_list = [evecs[:,1]]
    j = 2
    while isapprox(minE, newE)
        push!(ground_state_list, evecs[:,j])
        j += 1
        newE = evals[j]
    end
    n = length(ground_state_list)
    gs = sum(ground_state_list) / sqrt(n)

    return gs
end

function solve_closed_annealing_obj(H, init_state, tf)
    """
    Solves for closed-system dyanmics via Schrödinger equation.
    """
    annealing = Annealing(H, init_state)
    sol = solve_schrodinger(annealing, tf, alg=Exprb32());

    return sol
end

function solve_closed_annealing_obj(H, init_state, tf, alg)
    """
    Solves for closed-system dyanmics via Schrödinger equation.
    """
    annealing = Annealing(H, init_state)
    sol = solve_schrodinger(annealing, tf, alg=alg);

    return sol
end

function solve_closed_annealing_obj(H, init_state, tf, alg, sample_pts)
    """
    Solves for closed-system dyanmics via Schrödinger equation.
    """
    annealing = Annealing(H, init_state)
    sol = solve_schrodinger(annealing, tf, alg=alg, tstops=sample_pts);

    return sol
end

function solve_lindblad_dephasing(H, init_state, tf, γ)
    """
    Solves for closed-system dyanmics via Schrödinger equation.
    """
    # form initial state 
    ρ0 = init_state * init_state'
    # form dephasing operators
    n = log2(length(init_state))
    z1 = local_field_term(1, 1, n)
    lind_list = [Lindblad(γ, z1)]
    for j=2:n
        zj = local_field_term(1, j, n)
        lj = Lindblad(γ, zj)
        push!(lind_list, lj)
    end
    int_set = InteractionSet(lind_list...)
    lind_annealing = Annealing(H, ρ0, interactions = int_set)
    sol = solve_lindblad(lind_annealing, tf, alg=Tsit5());

    return sol
end

function solve_traj_ame_dephasing(H, init_state, tf, η, fc, T, num_fluc, num_traj)
    """
    Solves for closed-system dyanmics via Schrödinger equation.
    """
    # form dephasing operators
    n = log2(length(init_state))
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
    annealing = Annealing(H, init_state, interactions=int_set)
    prob = build_ensembles(annealing, tf, :ame)
    #t_list = range(0,tf,length=200)
    #sol = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, reltol=1e-6, saveat=t_list)
    sol = solve(prob, Tsit5(), EnsembleSerial(), trajectories=num_traj, reltol=1e-6)

    return sol
end

function solve_ame_dephasing(H, init_state, tf, η, fc, T)
    """
    Solves for closed-system dyanmics via Schrödinger equation.
    """
    # form dephasing operators
    n = log2(length(init_state))
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
    end
    int_set = InteractionSet(int_list...)
    annealing = Annealing(H, init_state, interactions=int_set)
    #t_list = range(0,tf,length=200)
    #sol = solve(prob, Tsit5(), EnsembleSerial(), trajectories=1, reltol=1e-6, saveat=t_list)
    sol = solve_ame(annealing, tf; alg=Tsit5(), ω_hint=range(-6, 6, length=100), reltol=1e-4)

    return sol
end
##

## ======================================
## Compute order parameters, fidelities
## ======================================
##
function compute_local_z_mag_field(phi::Vector{ComplexF64}, i::Int64)
    """
    Computes local z expectation value on qubit [i].
    """
    n = log2(length(phi))
    sz_i = local_field_term(1, i, n)
    exp_val = adjoint(phi) * sz_i * phi

    return exp_val
end

function compute_psuedo_spin(phi::Vector{ComplexF64}, i::Int64, j::Int64, k::Int64)
    """
    Computes ψⱼ 
    """
    m1 = compute_local_z_mag_field(phi, i)
    m2 = compute_local_z_mag_field(phi, j)
    m3 = compute_local_z_mag_field(phi, k)
    ψⱼ = m1 + m2 * exp(2*pi*im / 3) + m3 * exp(4*pi*im / 3)

    return ψⱼ
end

function compute_fidelity(phi, psi)
    amp = adjoint(phi) * psi
    return real(amp * conj(amp))
end

function compute_fidelity(ρ::Matrix{ComplexF64}, σ::Matrix{ComplexF64})
    """
    Computes fidelity between two density matrices.
    """
    sqrt_rho = sqrt(ρ)
    fid = real(tr(sqrt(sqrt_rho * σ * sqrt_rho))^2)

    return fid
end

function compute_fidelity(ψ::Vector{ComplexF64}, σ::Matrix{ComplexF64})
    """
    Computes fidelity between density matrix and pure state.
    """
    fid = real(ψ' * σ * ψ)

    return fid
end
##