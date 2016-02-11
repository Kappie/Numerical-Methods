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

class Measurement:
    # Coupling constant
    J = 1
    n_samples = 2 ** 14

    def __init__(self, observable, T, L, markov_chain):
        self.observable = observable
        self.T = T
        self.L = L
        self.N = L * L
        self.markov_chain = markov_chain

    def samples(self):
        samples = np.empty(self.n_samples)
        for i in range(self.n_samples):
            self.markov_chain.perform_sweep(self.N)
            samples[i] = self.observable(self.markov_chain.config)

        return samples

class MarkovChain:
    n_thermalization = 10000

    def __init__(self, initial_config, T):
        self.config = initial_config
        self.T = T
        self.thermalize()

    def thermalize(self):
        self.perform_sweep(n_thermalization)

    def perform_sweep(self, N):
        for i in range(N):
            self.next_config()

    def next_config(self):
        coordinate = self.random_coordinate()
        energy_difference = self.energy_difference_with_flipped_spin(coordinate)

        if energy_difference < 0 or random.uniform() < self.boltzmann_weight(energy_difference):
            self.flip_spin(coordinate)

    def random_coordinate(self):
        return (random.randint(L), random.randint(L))

    def flip_spin(self, coordinate):
        self.config[coordinate] *= -1

    # Calculates energy difference between config and config with spin flipped
    # at coordinate.
    def energy_difference_with_flipped_spin(self, coordinate):
        return 2 * J * self.config[coordinate] * np.sum(self.neighbors(*coordinate))

    def boltzmann_weight(self, energy_difference):
        return np.exp(-energy_difference / self.T)

    def neighbors(self, x, y):
        return [
            self.config[(x + 1) % L, y], self.config[x - 1, y], self.config[x, (y + 1) % L], self.config[x, y - 1]
        ]


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
    return np.abs(magnetization(configuration)) / configuration.size

def energy_per_site(configuration):
    return hamiltonian(configuration) / N

def magnetization_per_site(configuration):
    return magnetization(configuration) / N

def magnetization_squared(configuration):
    return magnetization(configuration) ** 2

def random_configuration():
    return 2 * random.random_integers(0, 1, (L, L)) - 1

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


### MAIN LOOP ###
T = 3
markov_chain = MarkovChain(initial_config = random_configuration(), T = T)
samples = Measurement(observable = order_parameter, T = T, L = 4, markov_chain = markov_chain).samples()

print np.mean(samples)
