from __future__ import division
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plot
import itertools
import pdb

T = 3   # Temperature
J = 1   # Coupling constant
L = 4   # System size
N = L * L

n_thermalization = 10000
n_sweeps = 2 ** 14

def hamiltonian(configuration):
    total_contribution = 0
    # Sum over each spin's neighbors. Need to multiple by 1/2 because I count every
    # interaction twice (?). Doesn't seem to be the case. I removed factor 1/2.
    for (x, y), spin in np.ndenumerate(configuration):
        total_contribution += np.sum( spin * neighbors(configuration, x, y) )

    return -J * total_contribution

def magnetization(configuration):
    return np.sum(configuration)

def order_parameter(configuration):
    return np.abs(magnetization(configuration)) / N

def energy_per_site(configuration):
    return hamiltonian(configuration) / N

def magnetization_per_site(configuration):
    return magnetization(configuration) / N

def magnetization_squared(configuration):
    return magnetization(configuration) ** 2

# Returns a 1D array of neighbours of given element. (Does not include
# the element itself.) Expects config to be a 2D array-like structure.
def neighbors(config, x, y, n=1):
    return [
        config[(x + 1) % L, y], config[x - 1, y], config[x, (y + 1) % L], config[x, y - 1]
    ]

def random_configuration():
    return 2 * random.random_integers(0, 1, (L, L)) - 1

def random_coordinate():
    return (random.randint(L), random.randint(L))

def flip_spin(config, coordinate):
    config[coordinate] *= -1

# Calculates energy difference between config and config with spin flipped
# at coordinate.
def energy_difference_with_flipped_spin(config, coordinate):
    return 2 * J * config[coordinate] * np.sum(neighbors(config, *coordinate))

def boltzmann_weight(energy_difference):
    return np.exp(-energy_difference / T)

def next_config(config):
    coordinate = random_coordinate()
    energy_difference = energy_difference_with_flipped_spin(config, coordinate)

    if energy_difference < 0 or random.uniform() < boltzmann_weight(energy_difference):
        flip_spin(config, coordinate)

    return config

def take_steps(config, n):
    for _ in range(n):
        config = next_config(config)

    return config



### MAIN LOOP ###
# Thermalize the markov chain
config = random_configuration()
config = take_steps(config, n_thermalization)

# Perform measurements
measurements_order_parameter = np.empty(n_sweeps)
measurements_energy_per_site = np.empty(n_sweeps)
measurements_magnetization_squared = np.empty(n_sweeps)
measurements_magnetization_per_site = np.empty(n_sweeps)
for i in range(n_sweeps):
    # Measure after performing N = L * L metropolis steps
    config = take_steps(config, N)
    measurements_order_parameter[i] = order_parameter(config)
    measurements_energy_per_site[i] = energy_per_site(config)
    measurements_magnetization_squared[i] = magnetization_squared(config)
    measurements_magnetization_per_site[i] = magnetization_per_site(config)

print "<E/N>", np.mean(measurements_energy_per_site)
print "<m>  ", np.mean(measurements_order_parameter)
print "<M^2>", np.mean(measurements_magnetization_squared)
print "<M/N>", np.mean(measurements_magnetization_per_site)
