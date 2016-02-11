from __future__ import division
import numpy as np
import matplotlib.pyplot as plot

def riemann_old(f, a, b, N=1000):
    delta_x = (b - a) / N
    sampling_points = np.linspace(a, b, N, endpoint=False)
    return np.sum( [delta_x * f(x) for x in sampling_points] )

def trapezoid_old(f, a, b, N=1000):
    delta_x = (b - a) / N
    return ( (delta_x/2) * (f(a)+f(b)) + np.sum( [delta_x * f(a + i*delta_x) for i in range(1, N) ] ) )

def riemann(f, a, b, N=1000):
    dx = (b - a) / N
    sampling_points = np.linspace(a, b, N, endpoint=False)
    return dx * np.sum(np.vectorize(f)(sampling_points))


def trapezoid(f, a, b, N=1000):
    delta_x = (b - a) / N
    return ( (delta_x/2) * (f(a)+f(b)) + np.sum( [delta_x * f(a + i*delta_x) for i in range(1, N) ] ) )

def monte_carlo(f, a, b, N=1000):
    return ((b - a) / N) * np.sum( f( (b-a) * np.random.random_sample(N) + a ) )

def plot_errors(method, N_values=np.arange(1000, 1000000, 1000)):
    errors = [relative_error(method(gaussian, -5, 5, N), 2.506626838) for N in N_values]
    plot.plot(N_values, errors, 'r^')
    plot.show()

def relative_error(estimate, true_value):
    return np.abs((estimate - true_value) / true_value)

def gaussian(x):
    return np.exp(-x*x/2)
