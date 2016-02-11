from __future__ import division
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plot
import itertools
import pdb

temperatures = np.array([0.125, 0.25, 0.5, 1, 2, 4])
L = 4
N = L * L
J = 1   # Coupling constant

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

def binning_analysis(samples):
    """
    Assumes samples.size is a power of 2. Returns an array with the naive monte carlo error estimate sigma / sqrt(N - 1)
    (assumption that samples are uncorrelated) for each binning order, starting with order 0 (i.e.
    all samples in their own bin) down to the minimum number of bins.
    """
    minimum_number_of_bins = 2 ** 6
    binning_steps = int( np.log2(samples.size / minimum_number_of_bins) ) + 1
    bins = np.copy(samples)

    errors = np.empty(binning_steps)
    for l in range(binning_steps):
        errors[l] = np.std(bins) / np.sqrt(bins.size - 1)
        # Take average of consecutive bins for next binning order
        bins = np.array( (bins[::2] + bins[1::2]) / 2 )

    return errors

def convergence_rate(errors):
    """
    Determines the average of the last few relative changes of an array of error estimates of a binning analysis.
    """
    relative_changes = errors[1:] - errors[:-1] / errors[1:]
    return np.mean(relative_changes[-3:])

def integrated_autocorrelation_time(bin_errors):
    return 0.5 * ( (errors[-1] / errors[0]) ** 2 - 1 )

def measurement(observable, temperatures, initial_config):
    measurements = np.empty(temperatures.size)

    for temp in temperatures:
        samples = np.empty(n_sweeps)

        for i in range(n_sweeps):
            config = take_steps(config, N)
            samples[i] = observable(config)

        return measurements

### MAIN LOOP ###
# Thermalize the markov chain
config = random_configuration()
config = take_steps(config, n_thermalization)

# Perform measurements
order_parameters = np.empty(temperatures.size)
# for T in temperatures:
#     measurements = np.empty(n_sweeps)
#
#     for i in range(n_sweeps):
#         # Measure after performing N = L * L metropolis steps
#         config = take_steps(config, N)
#         measurements[i] = order_parameter(config)
#         # measurements_energy_per_site[i] = energy_per_site(config)
#         # measurements_magnetization_squared[i] = magnetization_squared(config)
#         # measurements_magnetization_per_site[i] = magnetization_per_site(config)
#
#     bin_errors = binning_analysis(measurements)
#     plot.plot(bin_errors, 'ko')
#     plot.xlabel("binning step")
#     plot.ylabel("estimated error")
#     print "<m>: {m} +/- {error}.".format(m = np.mean(measurements), error = bin_errors[-1])
#     if convergence_rate(bin_errors) > 0.05:
#         print """Warning: binning analysis did not converge on error estimate.
#         Mean relative change of last few estimates was only {change}.""".format(change = convergence_rate(bin_errors))
#
# plot.show()

plot.plot([(1, 0), (2, 1), (3, 10), (4, 20)])
plot.show()

# print "<E/N>", np.mean(measurements_energy_per_site)
# print "<M^2>", np.mean(measurements_magnetization_squared)
# print "<M/N>", np.mean(measurements_magnetization_per_site)
