# import statements
##
using OpenQuantumTools, OrdinaryDiffEq, Plots, LaTeXStrings
using OpenQuantumBase
##

## Define input
##
# Initial state: the ground state of $H(0) = -Z1Z2-Z2Z3+Z1+Z2+Z3$, which is $\lvert 111 \rangle$
u0 = PauliVec[3][2] ⊗ PauliVec[3][2] ⊗ PauliVec[3][2]
##

## Read in Hamiltonian input data
##
coeff_array = [[] for x in 1:7]
open("/Users/worknic/Programming/AnnealingQGANs/code_from_leeseok/hamiltonian_coefficients.txt", "r") do f
    # read till end of file
    line = 1
    h_idx = 1
    while ! eof(f)
        # read a new / next line for every iteration          
        coeff_str = readline(f)
        coeffs = [parse(Float64, x) for x in split(coeff_str, ",")]
        push!(coeff_array[h_idx], coeffs)
        if mod(line, 6) == 0
            h_idx += 1
        end
        line += 1
    end
end
##

# construct interpolation for the 6 functions
##
tf = 5
tlist = range(0, tf, 50) / tf
ham_functions = []
for learned_coeffs in coeff_array
    interp_list = []
    for samp_f in learned_coeffs
        f = construct_interpolations(tlist, samp_f)
        push!(interp_list, f)
    end
    push!(ham_functions, interp_list)
end
##

## make plots
h_funcs = ham_functions[4]
labels = [L"f_1(t)" L"f_2(t)" L"f_3(t)" L"g_1(t)" L"g_2(t)" L"g_3(t)"]
data = vcat([h_funcs[j] for j=4:6] , [h_funcs[j] for j=1:3])
plot(tlist, data, label=labels)
xlabel!(L"s = t / t_f")
ylabel!("Coupling Strength")
savefig("theta_4_Hamps_result.pdf")
##

## Actually build the Hamiltonian 
# Hamiltonian: $H(t) = -Z1Z2 - Z2Z3 + e1(t)X1 + e2(t)X2 + e3(t)X3 + e4(t)Z1 + e5(t)Z2 + e6(t)Z3$
##
h_funcs = ham_functions[1]
# make time-dependent coefficient (function) array
tdc_array = vcat([(s) -> -1.0, (s) -> -1.0], h_funcs)
# make Hamiltonian
Z1Z2 = σz ⊗ σz ⊗ σi
Z2Z3 = σi ⊗ σz ⊗ σz
X1 = σx ⊗ σi ⊗ σi
X2 = σi ⊗ σx ⊗ σi
X3 = σi ⊗ σi ⊗ σx
Z1 = σz ⊗ σi ⊗ σi
Z2 = σi ⊗ σz ⊗ σi
Z3 = σi ⊗ σi ⊗ σz
op_list = [Z1Z2, Z2Z3, X1, X2, X3, Z1, Z2, Z3]
H = DenseHamiltonian(tdc_array, op_list, unit=:ħ)
##

## Define Annealing object and solve
##
tf = 5
annealing = Annealing(H, u0)
closed_sol = solve_schrodinger(annealing, tf, alg=Tsit5(), abstol=1e-6, reltol=1e-6);
##

## solve von-neumann
##
ρ0 = u0*u0'
annealing = Annealing(H, ρ0)
vn_sol = solve_von_neumann(annealing, tf, alg = Tsit5(), abstol=1e-6, reltol=1e-6);
##

## Read in supposed final state
##
state_array = []
open("/Users/worknic/Programming/AnnealingQGANs/code_from_leeseok/test_states.txt", "r") do f
    # read till end of file
    while ! eof(f)
        # read a new / next line for every iteration          
        state_str = readline(f)
        coeffs = [parse(ComplexF64, x) for x in split(state_str, ",")]
        print(coeffs)
        push!(state_array, coeffs)
    end
end
##

## Compare closed system solution to desired state
##
closed_fid = abs(dot(state_array[1], closed_sol[end]))^2
##

## See what happens when noise added via Redfield
##
coupling = ConstantCouplings(["ZII", "IZI", "IIZ"], unit=:ħ)
bath = Ohmic(1e-4, 4, 16)
ρ0 = u0*u0'
annealing = Annealing(H, ρ0; coupling=coupling, bath=bath)
U = solve_unitary(annealing, tf, alg=Tsit5(), abstol=1e-8, retol=1e-8);
red_sol = solve_redfield(annealing, tf, U; alg=Tsit5(), abstol=1e-8, retol=1e-8)
##

## Compare desired to Redfield
##
red_fid = real(state_array[1]' * red_sol[end] * state_array[1])
##

## See what happens if we add some noise via Lindblad equation
##
dephasing = (Lindblad(0.1, σz⊗σi⊗σi), Lindblad(0.1, σi⊗σz⊗σi), Lindblad(0.1, σi⊗σi⊗σz))
int_set = InteractionSet(dephasing)
# combine them into an Annealing object
lind_annealing = Annealing(H, ρ0, interactions = int_set)
# open dynamics solution
lind_sol = solve_lindblad(lind_annealing, tf, alg=Tsit5());
##

## Compare desired to Lindblad
##
lind_fid = real(state_array[1]' * lind_sol[end] * state_array[1])
##


##
## Define Noise 1 / f spectral density
##
function build_custom_1f_spectrum(amplitude, low_cut_off, high_cut_off, temperature, α)
    β = temperature_2_β(temperature)
    function spectrum(ω)
        if ω>high_cut_off
            0
        elseif ω>low_cut_off
            amplitude / (ω)^α
        elseif ω>0
            amplitude / (low_cut_off)^α
        elseif ω>-low_cut_off
            exp(β*ω) * amplitude / (low_cut_off)^α
        elseif ω>-high_cut_off
            -exp(β*ω) * amplitude / (ω)^α
        else
            0
        end
    end
end;

amp = 1e-3
fl = 1e-4
fh = 1
T = 12
α = 1

Sf = build_custom_1f_spectrum(amp, 2*π*fl, 2*π*fh, T, α);
##

##
## Combine Ohmic + 1 / f into single bath spectrum
##
η=1e-4; fc=2*π*4; T=12;
ohmic_bath = Ohmic(η, fc, T)
function combine_spectrum(ohmic_object, Sf)
    (ω) -> spectrum(ω, ohmic_object) + Sf(ω)
end
combined_spectrum = combine_spectrum(ohmic_bath, Sf);

w_axis = 2*π*range(-10, 75, length=10000)
plot(w_axis, combined_spectrum.(w_axis), xlabel=L"\omega", ylabel=L"S(\omega)", label="", linewidth=2)
#savefig("spectral_density.pdf")
##

## Blah
##
c1 = ConstantCouplings(["ZII"], unit=:ħ)
b1 = CustomBath(spectrum=combined_spectrum);
interaction_1 = Interaction(c1, b1)
c2 = ConstantCouplings(["IZI"], unit=:ħ)
b2 = CustomBath(spectrum=combined_spectrum);
interaction_2 = Interaction(c2, b2)
c3 = ConstantCouplings(["IIZ"], unit=:ħ)
b3 = CustomBath(spectrum=combined_spectrum);
interaction_3 = Interaction(c3, b3);
intset = InteractionSet(interaction_1, interaction_2, interaction_3)
##

##
## Define custom bath object and solve AME
##
annealing = Annealing(H, u0, interactions=intset)
ame_sol2 = solve_ame(annealing, tf, alg=Tsit5(), ω_hint=range(-2,2,length=100), abstol=1e-6, reltol=1e-6)
##

# Pre-compute Lamb-shift and Davies Operators 
##
L = OpenQuantumBase.davies_from_interactions(intset, range(-2,2,length=100), true, nothing)
##

##
## Compare to AME
##
ame_fid2 = real(state_array[1]' * ame_sol2[end] * state_array[1])
##
##



##
## Define custom bath object and solve AME
##
coupling = ConstantCouplings(["ZII", "IZI", "IIZ"], unit=:ħ)
bath = CustomBath(spectrum=combined_spectrum);
annealing = Annealing(H, u0, coupling=coupling, bath=bath)
ame_sol = solve_ame(annealing, tf, alg=Tsit5(), ω_hint=range(-2,2,length=100), abstol=1e-6, reltol=1e-6)
##

##
## Compare to AME
##
ame_fid = real(state_array[1]' * ame_sol[end] * state_array[1])
##
##

## See what happens if we add some noise via Lindblad equation
##
γ = 0.35
dephasing = (Lindblad(γ, σz⊗σi⊗σi), Lindblad(γ, σi⊗σz⊗σi), Lindblad(γ, σi⊗σi⊗σz))
int_set = InteractionSet(dephasing)
ρ0 = u0*u0'
# combine them into an Annealing object
lind_annealing = Annealing(H, ρ0, interactions = int_set)
# open dynamics solution
lind_sol = solve_lindblad(lind_annealing, tf, alg=Tsit5());
lind_fid = real(state_array[1]' * lind_sol[end] * state_array[1])
##

## Compare desired to Lindblad
##
lind_fid = real(state_array[1]' * lind_sol[end] * state_array[1])
##



##
## Define custom bath object and solve Redfield (didn't work--correlation not specified)
##
coupling = ConstantCouplings(["ZII", "IZI", "IIZ"], unit=:ħ)
bath = CustomBath(spectrum=combined_spectrum);
annealing = Annealing(H, u0, coupling=coupling, bath=bath)
U = solve_unitary(annealing, tf, alg=Tsit5(), abstol=1e-8, retol=1e-8);
red_sol = solve_redfield(annealing, tf, U; alg=Tsit5(), ω_hint=range(-2,2,length=100), abstol=1e-8, retol=1e-8)
##
##
## Compare to Redfield
##
red_fid = real(state_array[1]' * ame_sol[end] * state_array[1])
##
##

## Wrap together
# Hamiltonian: $H(t) = -Z1Z2 - Z2Z3 + e1(t)X1 + e2(t)X2 + e3(t)X3 + e4(t)Z1 + e5(t)Z2 + e6(t)Z3$
##
tf = 5
tlist = range(0, tf, 50) / tf
ham_functions = []
for learned_coeffs in coeff_array
    interp_list = []
    for samp_f in learned_coeffs
        f = construct_interpolations(tlist, samp_f)
        push!(interp_list, f)
    end
    push!(ham_functions, interp_list)
end
h_s_idx = 7
state = state_array[h_s_idx]
h_funcs = ham_functions[h_s_idx]
# make time-dependent coefficient (function) array
tdc_array = vcat([(s) -> -1.0, (s) -> -1.0], h_funcs)
# make Hamiltonian
Z1Z2 = σz ⊗ σz ⊗ σi
Z2Z3 = σi ⊗ σz ⊗ σz
X1 = σx ⊗ σi ⊗ σi
X2 = σi ⊗ σx ⊗ σi
X3 = σi ⊗ σi ⊗ σx
Z1 = σz ⊗ σi ⊗ σi
Z2 = σi ⊗ σz ⊗ σi
Z3 = σi ⊗ σi ⊗ σz
op_list = [Z1Z2, Z2Z3, X1, X2, X3, Z1, Z2, Z3]
H = DenseHamiltonian(tdc_array, op_list, unit=:ħ)
annealing = Annealing(H, u0)
closed_sol = solve_schrodinger(annealing, tf, alg=Tsit5(), abstol=1e-6, reltol=1e-6);
closed_fid = abs(dot(state, closed_sol[end]))^2
print("Closed system fidelity for state_$(h_s_idx) is: $(closed_fid)\n")

coupling = ConstantCouplings(["ZII", "IZI", "IIZ"], unit=:ħ)
bath = CustomBath(spectrum=combined_spectrum);
annealing = Annealing(H, u0, coupling=coupling, bath=bath)
ame_sol = solve_ame(annealing, tf, alg=Tsit5(), ω_hint=range(-2,2,length=100), abstol=1e-6, reltol=1e-6)
ame_fid = real(state' * ame_sol[end] * state)
print("AME fidelity for state_$(h_s_idx) is: $(ame_fid)\n")

c1 = ConstantCouplings(["ZII"], unit=:ħ)
b1 = CustomBath(spectrum=combined_spectrum);
interaction_1 = Interaction(c1, b1)
c2 = ConstantCouplings(["IZI"], unit=:ħ)
b2 = CustomBath(spectrum=combined_spectrum);
interaction_2 = Interaction(c2, b2)
c3 = ConstantCouplings(["IIZ"], unit=:ħ)
b3 = CustomBath(spectrum=combined_spectrum);
interaction_3 = Interaction(c3, b3);
intset = InteractionSet(interaction_1, interaction_2, interaction_3)
annealing = Annealing(H, u0, interactions=intset)
ame_sol2 = solve_ame(annealing, tf, alg=Tsit5(), ω_hint=range(-2,2,length=100), abstol=1e-6, reltol=1e-6)
ame_fid2 = real(state' * ame_sol2[end] * state)
print("AME fidelity separate baths for state_$(h_s_idx) is: $(ame_fid2)\n")

γ = 0.35
dephasing = (Lindblad(γ, σz⊗σi⊗σi), Lindblad(γ, σi⊗σz⊗σi), Lindblad(γ, σi⊗σi⊗σz))
int_set = InteractionSet(dephasing)
# combine them into an Annealing object
lind_annealing = Annealing(H, ρ0, interactions = int_set)
# open dynamics solution
lind_sol = solve_lindblad(lind_annealing, tf, alg=Tsit5());
lind_fid = real(state' * lind_sol[end] * state)
print("Lindblad fidelity for state_$(h_s_idx) is: $(lind_fid)\n")
##
##