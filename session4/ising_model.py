from __future__ import division
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plot
import itertools
import pdb

class IsingModel:
    """
    Two-dimensional Ising model.
    """

    J = 1   # Coupling constant

    def __init__(self, T, L, markov_chain):
        self.T = T
        self.L = L
        self.N = L * L
        self.lattice = Lattice(L)
        self.markov_chain = markov_chain(self)


    def perform_sweep(self):
        self.markov_chain.perform_sweep()

    # Calculates energy difference between config and config with spin flipped
    # at coordinate.
    def energy_difference_with_flipped_spin(self, coordinate):
        return 2 * self.J * self.lattice.spin_at(coordinate) * np.sum(self.lattice.neighbors(coordinate))

    def boltzmann_weight(self, energy_difference):
        return np.exp(-energy_difference / self.T)


class Lattice:
    """
    Two-dimensional spin lattice with periodic boundary conditions. Initializes randomly.
    Spins are represented by 1 and -1.
    """

    # Coupling constant
    J = 1

    def __init__(self, L):
        self.L = L
        self.N = L * L
        self.config = self.random_configuration()

        self.order_parameter = self.calculate_order_parameter()
        self.energy_per_site = self.calculate_energy_per_site()
        self.magnetization_per_site = self.calculate_magnetization_per_site()
        self.magnetization = self.calculate_magnetization()

    def spin_at(self, coordinate):
        return self.config[coordinate]

    def flip_spin(self, coordinate):
        self.config[coordinate] *= -1
        self.update_observables(coordinate)

    def neighbors(self, coordinate):
        x, y = coordinate
        return [
            self.config[(x + 1) % self.L, y],
            self.config[x - 1, y],
            self.config[x, (y + 1) % self.L],
            self.config[x, y - 1]
        ]

    def update_observables(self, coordinate):
        self.update_order_parameter(coordinate)
        self.update_energy_per_site(coordinate)
        self.update_magnetization_per_site(coordinate)
        self.update_magnetization(coordinate)

    def update_order_parameter(coordinate):

    def update_energy_per_site(coordinate):


    def random_configuration(self):
        return random.choice([-1, 1], [self.L, self.L])

    def calculate_order_parameter(self):
        return np.abs(self.magnetization()) / self.N

    def calculate_energy_per_site(self):
        return self.hamiltonian() / self.N

    def calculate_magnetization_per_site(self):
        return self.magnetization() / self.N

    def calculate_magnetization_squared(self):
        return self.magnetization() ** 2

    def calculate_hamiltonian(self):
        total_contribution = 0

        for coordinate, spin in np.ndenumerate(self.config):
            total_contribution += np.sum( spin * self.neighbors(coordinate) )

        return -self.J * total_contribution

    def calculate_magnetization(self):
        return np.sum(self.config)

class Wolff:
    """
    Implements #perform_sweep method as only public method. Uses the Wolff algorithm to decide the next configuration of the
    ising model configuration.
    """

    def __init__():
        return "hoi"

class Metropolis:
    """
    Implements #perform_sweep method as only public method. Uses the metropolis algorithm to decide the next configuration of the ising model configuration.
    """

    thermalization_sweeps = 5000

    def __init__(self, model):
        """
        This class has a lot of "train wrecks", e.g. self.model.lattice.flip_spin. How to modularize this?
        Maybe it is unavoidable, since the markov chain and the model are so tightly bound.
        """
        self.model = model

    def thermalize(self):
        for _ in range(self.thermalization_sweeps):
            self.perform_sweep()

    def perform_sweep(self):
        self.take_steps(self.model.N)

    def take_steps(self, n):
        for _ in range(n):
            self.take_step()

    def take_step(self):
        coordinate = self.random_coordinate()
        energy_difference = self.model.energy_difference_with_flipped_spin(coordinate)

        if energy_difference < 0 or random.uniform() < self.model.boltzmann_weight(energy_difference):
            self.model.lattice.flip_spin(coordinate)

    def random_coordinate(self):
        return (random.randint(self.model.L), random.randint(self.model.L))


class Measurement:
    """
    Samples the Monte Carlo Ising model and performs binning analysis.
    """
    n_samples = 2 ** 14

    def __init__(self, model):
        self.model = model
        self.get_samples()
        self.binning_analysis()

    # Only sample order parameter for now.
    def get_samples(self):
        self.samples = np.empty(self.n_samples)

        for i in range(self.n_samples):
            self.model.perform_sweep()
            self.samples[i] = self.model.order_parameter()

        self.expectation_value = np.mean(self.samples)
        return self.samples

    def binning_analysis(self):
        """
        Assumes samples.size is a power of 2. Returns an array with the naive monte carlo error estimate sigma / sqrt(N - 1)
        (assumption that samples are uncorrelated) for each binning order, starting with order 0 (i.e.
        all samples in their own bin) down to the minimum number of bins.
        """
        minimum_number_of_bins = 2 ** 6
        binning_steps = int( np.log2(self.samples.size / minimum_number_of_bins) ) + 1
        bins = np.copy(self.samples)

        self.errors = np.empty(binning_steps)
        for l in     range(binning_steps):
            self.errors[l] = np.std(bins) / np.sqrt(bins.size - 1)
            # Take average of consecutive bins for next binning order
            bins = np.array( (bins[::2] + bins[1::2]) / 2 )

        self.std = self.errors[-1]

        if self.convergence_rate() > 0.05:
            self.error_message = """Error estimates of binning analysis did not converge.
            Convergence rate was only {rate}.""".format(rate = self.convergence_rate())
            print self.error_message

        return self.errors

    def convergence_rate(self):
        """
        Determines the average of the last few relative changes of an array of error estimates of a binning analysis.
        """
        relative_changes = self.errors[1:] - self.errors[:-1] / self.errors[1:]
        return np.mean(relative_changes[-3:])

    def integrated_autocorrelation_time(self):
        return 0.5 * ( (self.errors[-1] / self.errors[0]) ** 2 - 1 )

temperatures = [3]
L = 4
expectation_values = []
stds = []

for T in temperatures:
    model = IsingModel(T, L, Metropolis)
    measurement = Measurement(model)
    expectation_values.append(measurement.expectation_value)
    stds.append(measurement.std)

print expectation_values
print stds
plot.errorbar(temperatures, expectation_values, yerr = stds)
plot.show()
