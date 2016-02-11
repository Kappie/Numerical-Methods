from __future__ import division
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plot

delta = 0.02
n_thermalization = 10000
n_steps = 2 ** 18
# Denoted by l in the lecture notes
max_bin_order = 17

left_bound = 0
right_bound = 1
normalization = 1 / 0.632121

correct_value = 0.418023

def A(x):
    return x

def p(x):
    return np.exp(-x)

def change():
    return random.uniform(-delta, delta)

def next_point(x):
    x_new = x + change()

    # Reject new points that fall outside the integration domain
    while x_new < left_bound or x_new > right_bound:
        x_new = x + change()

    ratio = p(x_new) / p(x)

    if ratio > 1 or random.uniform() < ratio:
        return x_new
    else:
        return x

def perform_steps(x_init, steps):
    x = x_init
    for _ in range(steps):
        x = next_point(x)
    return x

# Initialize point randomly in the integration interval
x = random.uniform(left_bound, right_bound)

# Allow the markov chain to converge to the desired distribution p(x).
# Do not make any m easurements yet.
x = perform_steps(x, n_thermalization)

# Calculate expectation value of A.
# Caution: np.empty(n) initializes with random values. Every entry should be
# manually overwritten.
expectation_values = []
for _ in range(5):
    measurements = np.empty(n_steps)
    for i in range(n_steps):
        x = next_point(x)
        measurements[i] = A(x)

    expectation_value = np.mean(measurements)
    expectation_values.append(expectation_value)


print "standard deviation measured:"
print np.std(expectation_values)
print "what standard deviation should be naively: {std}".format(std = np.std(measurements) / np.sqrt(n_steps - 1))
print "expectation value measured: {exp}, correct_value: {correct}".format(exp = expectation_value, correct = correct_value)

# Split measurements into increasing amount of bins and calculate standard
# deviation of bin averages.
estimated_errors = np.empty(max_bin_order)
for bin_order in range(max_bin_order):
    bins = np.split(measurements, n_steps / 2 ** bin_order)
    bin_averages = np.mean(bins, axis=1)
    # ddof = 1 implements Bessel's correction (?)
    estimated_errors[bin_order] = np.std(bin_averages) / np.sqrt(bin_averages.size)

print "converging values of the standard deviation of the bin averages:"
print estimated_errors
