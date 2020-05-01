import numpy as np
import random
from matplotlib import pyplot as plt

'''
TASK 1
'''
print('Starting task 1')

def uniform(N):
    sample = [random.uniform(0, 1) for i in range(N)]
    x = np.asarray(sample)
    f = 1 / (1 + x ** 3)
    f2 = (1 / (1 + x ** 3)) ** 2
    I_N = np.sum(f) / N
    variance = np.sum(f2) / N - I_N ** 2
    error = np.sqrt(variance)
    print('N=', N, ': \tI_N=', I_N, ', error=', error)
print('--UNIFORM DISTRIBUTION--')
uniform(10)
uniform(100)
uniform(1000)
uniform(10000)

def linear(N, a=1.0, b=0.5):
    sample = [random.uniform(0, np.sqrt((2*a+b)/b)-b) for i in range(N)]    # not quite sure about np.sqrt((2*a+b)/b)-b, got there by trial and error
    for i in range(len(sample)):
        sample[i] = (np.sqrt(b ** 2 + a * sample[i]) - b) / a               # inverse of primitive fctn of ax+b
    x = np.asarray(sample)
    f = 1 / (1 + x ** 3)
    f2 = (1 / (1 + x ** 3)) ** 2
    I_N = np.sum(f) / N
    variance = np.sum(f2) / N - I_N ** 2
    error = np.sqrt(variance)
    print('N=', N, ': \tI_N=', I_N, ', error=', error)
print('--LINEAR DISTRIBUTION--')
linear(10)
linear(100)
linear(1000)
linear(10000)

def metropolis(N, a=1.0, b=0.5):
    sample = []
    x_i = random.uniform(0, 1)
    sample.append(x_i)
    # half is to be accepted: delta=4*(1-x_i) if x_i>2/3, delta=2-x_i otherwise
    delta = 2-x_i               # see above
    # populate sample
    while len(sample) < N:
        if x_i > 2/3:
            delta = 4*(1-x_i)   # see above
        # generate trial change
        r = random.uniform(0,1)
        x_j = x_i+delta*(r-0.5)
        # determine whether x_j is to be accepted
        if x_j < 0 or x_j > 1:
            continue
        r = random.uniform(0,1)
        if r < (a*x_j+b)/(a*x_i+b): # note that this entails x_j > x_i
            x_i = x_j
            sample.append(x_i)
            delta = 2-x_i
    # calculate integral
    x = np.asarray(sample)
    f = 1 / (1 + x ** 3)
    I_N = np.sum(f) / N
    # calculate autocorrelation
    k = -1
    psi_k = 1                           # arbitrary value larger than 0.1
    while np.abs(psi_k) > 0.1 and k < N-1: # np.abs really needed??
        k += 1
        w = a*x+b
        g = f/w
        corr = []
        for i in range(len(x)-k):
            corr.append(g[i]*g[i+k])
        corr_mean = np.mean(corr)
        g_mean = np.mean(g)
        g_square_mean = np.mean(g*g)
        psi_k = (corr_mean-g_mean*g_mean)/(g_square_mean-g_mean*g_mean)
    # print results
    print('N=', N, ': \tI_N=', I_N, ', psi_k=', psi_k, 'for k=', k)
print('--METROPOLIS--')
metropolis(10)
metropolis(100)
metropolis(1000)
metropolis(10000)

'''
bins = [0.1 * i for i in range(12)]
sample = [random.uniform(0,2/a+b) for i in range(25000)]    # not quite sure about 2/a+b, got there by trial and error
for i in range(len(sample)):
    sample[i] = (np.sqrt(b**2+a*sample[i])-b)/a             # inverse of primitive fctn of ax+b
plt.hist(sample, bins)
plt.savefig("1-triangular.png")
plt.clf()
'''

