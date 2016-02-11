from __future__ import division
import numpy as np
import matplotlib.pyplot as plot

def gaussian(x):
    return np.exp(-x*x)

def monte_carlo(f, a, b, N=1000):
    return ((b - a) / N) * np.sum( f( (b-a) * np.random.random_sample(N) + a ) )

def monte_carlo_importance_sampling(f, a, b, N=1000):
    # Sample randomly from e^-x
    # left and right bounds a = 0 and b = infinity are assumed for calculating the average value of
    # f(x)e^x. Left bound should be zero, and b large for this guess to be about right.
    random_samples = -np.log(1 - np.random.random_sample(N))
    average_value = (1 / N) * np.sum( f(random_samples) * np.exp(random_samples) )
    return average_value

def error(estimate, true_value):
    return estimate - true_value

errors = [ error(monte_carlo(gaussian, 0, 100, 1000), 0.88623) for _ in range(1000) ]
errors_importance_sampling = [ error(monte_carlo_importance_sampling(gaussian, 0, 100, 1000), 0.88623) for _ in range(1000) ]

print "mean: {mean}, std: {std}".format(mean = np.mean(errors), std = np.std(errors))
print "(with importance sampling) mean: {mean}, std: {std}".format(mean = np.mean(errors_importance_sampling), std = np.std(errors_importance_sampling))
