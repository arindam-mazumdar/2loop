import numpy as np
import dynesty
from F_kernels_interface import pyF2sym

# Define the dimensionality of our problem.
ndim = 3

def ptform(u):
    return 1.* u


# Sample from our distribution.
sampler = dynesty.NestedSampler(pyF2sym, ptform, ndim,
                                bound='single', nlive=100)
sampler.run_nested(dlogz=0.1)
res = sampler.results
