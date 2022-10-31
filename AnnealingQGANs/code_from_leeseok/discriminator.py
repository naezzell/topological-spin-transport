import numpy as np
import picos
import mosek
import gc
#import pennylane_utils
import itertools
import math
import cvxopt

def discriminator(rho_A, rho_B, norm, nq):
    
    P = picos.Constant(rho_A)
    Q = picos.Constant(rho_B)

    F = picos.Problem()
    O = picos.HermitianVariable("O", 2**nq)

    F.set_objective("max", picos.trace(O*(P-Q)))
    
    if norm == 'trace':
        F.add_constraint(1>=picos.NuclearNorm(O))
    elif norm == 'l2':
        F.add_constraint(1>=picos.norm(O))
    elif norm == 'operator':
        F.add_constraint(1>=picos.expressions.exp_specnorm.SpectralNorm(O))

    F.solve(solver = "mosek")
    print("maximum Trace is :", round(F,4))
        
    return np.array(F.get_valued_variable("O"))