from __future__ import division
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plot
import itertools

T = 3   # Temperature
J = 1   # Coupling constant
L = 4   # System size
N = L * L

n_thermalization = 10000
n_sweeps = 2 ** 14

def hamiltonian(configuration):
    total_contribution = 0
    for (x, y), spin in np.ndenumerate(configuration):
        total_contribution += np.sum( spin * neighbors(configuration, x, y) )

    return -J * total_contribution

def magnetization(configuration):
    return np.sum(configuration)

def order_parameter(configuration):
     np.abs(magnetization(configuration)) / configuration.size

def energy_per_site(configuration):
    hamiltonian(configuration) / configuration.size

def magnetization_squared(configuration):
    magnetization(configuration) ** 2

# Returns a 1D array of neighbours of given element. (Does not include
# the element itself.) Expects config to be a 2D array-like structure.
# Use n = 1 for nearest neighbors, n = 2 to include next-to-nearest, etc.
def neighbors(config, x, y, n=1):
    n_rows, n_columns = config.shape
    # Do not include the element itself.
    indices = [ index for index in itertools.product(range(-n, n+1), range(-n, n+1)) if index != (0, 0) ]
    return np.array([ config[(x+a) % n_rows, (y+b) % n_columns] for a, b in indices ])

def random_configuration():
    return 2 * random.random_integers(0, 1, (L, L)) - 1

def flip_random_spin(config):
    x, y = random.randint(L), random.randint(L)
    config[x, y] *= -1

def p(energy_difference):
    return np.exp(- energy_difference / T)

def next_config(config):
    # This can be optimized, heavily. Now I'm copying an array and looping over all neighbors twice.
    new_config = np.copy(config)
    flip_random_spin(new_config)
    energy_difference = hamiltonian(new_config) - hamiltonian(config)
    if energy_difference < 0 or random.uniform() < p(energy_difference):
        return new_config
    else:
        return config

def take_steps(config, n):
    for _ in range(n):
        config = next_config(config)

    return config

### MAIN LOOP ###
config = random_configuration()
