##
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase
##

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
                print(c)
                push!(couplings_list, c)
            end
        end
    end
    coeffs = [J for _=1:length(couplings_list)]
    zz_part = two_local_term(coeffs, couplings_list, n)
    
    return x_part, zz_part
end
##

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
##

##
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
##

## 
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

function get_initial_ground_state(H)
    gs_info = eigen_decomp(H, 0; lvl=1)
    gs_vec = vec(gs_info[2])

    return gs_vec
end
##

##
function solve_closed_annealing_obj(H, init_state, tf)
    """
    Solves for closed-system dyanmics via Schrödinger equation.
    """
    annealing = Annealing(H, init_state)
    sol = solve_schrodinger(annealing, tf, alg=Tsit5());

    return sol
end
##

##
# Run example
##

## Build up all the pre-solution moving parts
##
J = 1
Γ = 0.1
adj_mat = build_triangular_chain_adjacency_matrix(1);
# build static Hamiltonian
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ)
# build time-depndent function
pause_time = 10;
period = 1;
amp = 1;
cycles = 5;
func, tf, svals = form_pause_into_comb_function(pause_time, period, amp, cycles);
# get full Hamiltonian
H = form_dynamic_hamiltonian(x_part, z_part, func)
gs = get_initial_ground_state(H)
##

## Solve the Schrodinger equation
##
sol = solve_closed_annealing_obj(H, gs, tf)
##

## Compute local z magnetizations
##
function compute_local_z_mag_field(phi::Vector{ComplexF64}, i::Int64)
    """
    Computes local z expectation value on qubit [i].
    """
    n = log2(length(gs))
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

phi_gs_up = q_translate_state("(001)+(101)", normal=true)
compute_psuedo_spin(phi_gs_up, 1, 2, 3)

phi_es_right = q_translate_state("011")
compute_psuedo_spin(phi_es_right, 1, 2, 3)

# down = 1; up = 0
# hex states lie on outer hexagon
# down, up, up
hex0 = q_translate_state("100")
psi_hex0 = compute_psuedo_spin(hex0, 1, 2, 3)
# down up down
hex1 = q_translate_state("101")
psi_hex1 = compute_psuedo_spin(hex1, 1, 2, 3)
# up up down 
hex2 = q_translate_state("001")
psi_hex2 = compute_psuedo_spin(hex2, 1, 2, 3)
# up down down
hex3 = q_translate_state("011")
psi_hex3 = compute_psuedo_spin(hex3, 1, 2, 3)
# up down up
hex4 = q_translate_state("101")
psi_hex4 = compute_psuedo_spin(hex4, 1, 2, 3)
# down down up
hex5 = q_translate_state("110")
psi_hex5 = compute_psuedo_spin(hex5, 1, 2, 3)

# down = 1; up = 0; right = +; left = -
# circ states lie on inner circle
# down up right
circ0 = q_translate_state("(100)+(101)", normal=true)
psi_circ0 = compute_psuedo_spin(circ0, 1, 2, 3)
# right up down
circ1 = q_translate_state("(001)+(101)", normal=true)
psi_circ1 = compute_psuedo_spin(circ1, 1, 2, 3)
# up right down
circ2 = q_translate_state("(001)+(011)", normal=true)
psi_circ2 = compute_psuedo_spin(circ2, 1, 2, 3)
# up down right
circ3 = q_translate_state("(010)+(011)", normal=true)
psi_circ3 = compute_psuedo_spin(circ3, 1, 2, 3)
# right down up
circ4 = q_translate_state("(010)+(110)", normal=true)
psi_circ4 = compute_psuedo_spin(circ4, 1, 2, 3)
# down right up
circ5 = q_translate_state("(100)+(110)", normal=true)
psi_circ5 = compute_psuedo_spin(circ5, 1, 2, 3)
##
##
psuedo_spin_list = []
for i=1:length(sol)
    state = sol[i]
    ps = compute_psuedo_spin(state, 1, 2, 3)
    push!(psuedo_spin_list, ps)
end
##

##
plot(sol.t, real.(psuedo_spin_list), label="real", xlabel="time (ns)", ylabel="psuedo-spin value", title="AC field drives psuedo-spin")
plot!(sol.t, imag.(psuedo_spin_list), label="imag")
##

##
plot(svals * tf, func(svals), legend=false, xlabel="time (ns)", ylabel="h_1 driving amplitude")
##

##
guess_gs = sqrt(1 / 12) * (circ0 + circ1 + circ2 + circ3 + circ4 + circ5)
##

## 
function compute_fid(phi, psi)
    amp = adjoint(phi) * psi
    return real(amp * conj(amp))
end
##

## prepare thermal state by direct matrix exp
##
temp = 10
β = 2 * π / temp
mat_exp = Matrix(exp(-β * H(0)))
ρ_gibbs = mat_exp / tr(mat_exp)

d = size(ρ_gibbs)[1]
ρ_mixed = (1 / d) * diagm(ones(d)) .+ 0im
##

##
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
##

## prepare thermal state by relaxation
##
J = 1
Γ = 0.1
adj_mat = build_triangular_chain_adjacency_matrix(1);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ)
pause_time = 1000;
func_pause, tf_pause, svals_pause = form_pause_function(pause_time);
H_pause = form_dynamic_hamiltonian(x_part, z_part, func_pause);
gs_pause = get_initial_ground_state(H_pause)
##

## Solve for Lindblad dynamics
γ = 10
sol_lind = solve_lindblad_dephasing(H_pause, gs_pause, tf_pause, γ)
##

##
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

##
fid_list = [compute_fidelity(ρ_mixed, sol_lind[j]) for j=1:length(sol_lind)]
plot(sol_lind.t, fid_list, legend=false, xlabel="time (ns)", ylabel="fidelity with mixed state", title="Lindblad prepares mixed state")
##

##
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

## prepare thermal state by relaxation
##
J = 1
Γ = 0.1
adj_mat = build_triangular_chain_adjacency_matrix(1);
x_part, z_part = form_static_ham_ops(adj_mat, J, Γ)
pause_time = 100;
func_pause, tf_pause, svals_pause = form_pause_function(pause_time);
H_pause = form_dynamic_hamiltonian(x_part, z_part, func_pause);
gs_pause = get_initial_ground_state(H_pause)
##

## Solve for AME dynamics
η = 1e-2;
fc = 4;
sol_ame = solve_ame_dephasing(H, gs_pause, tf_pause, η, fc, temp)
##

##
fid_list = [compute_fidelity(ρ_gibbs, sol_ame[j]) for j=1:length(sol_ame)]
plot(sol_ame.t, fid_list, legend=false, xlabel="time (ns)", ylabel="fidelity with Gibbs state", title="AME prepares Gibbs state")
##