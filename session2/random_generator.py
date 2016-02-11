from __future__ import division
import numpy as np
import matplotlib.pyplot as plot

def generator(a, c, m, seed):
    next_seed = (a * seed + c) % m
    return next_seed

def generate_random_numbers(amount, a, c, m, initial_seed):
    random_numbers = []
    next_seed = generator(a, c, m, initial_seed)

    for _ in range(amount):
        random_numbers.append(next_seed / m)
        next_seed = generator(a, c, m, next_seed)

    return random_numbers

# expects an even number of numbers
def square_test(numbers):
    x = numbers[::2]
    y = numbers[1::2]
    plot.scatter(x, y)

random_numbers = generate_random_numbers(amount = 100000, a = 16807, m = 2**31 - 1, c = 0, initial_seed = 667790)

#plot.hist(random_numbers)
square_test(random_numbers)
plot.show()
