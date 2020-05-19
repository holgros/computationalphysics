
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Fixing random state for reproducibility
np.random.seed(19680801)


def randrange(n, vmin, vmax):
    '''
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    '''
    return (vmax - vmin)*np.random.rand(n) + vmin

'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

n = 100

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
for m, zlow, zhigh in [('o', -50, -25), ('^', -30, -5)]:
    xs = randrange(n, 23, 32)
    ys = randrange(n, 0, 100)
    zs = randrange(n, zlow, zhigh)
    ax.scatter(xs, ys, zs, marker=m)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
'''
'''
bins = [0.1 * i for i in range(12)]
sample = [random.uniform(0,2/a+b) for i in range(25000)]    # not quite sure about 2/a+b, got there by trial and error
for i in range(len(sample)):
    sample[i] = (np.sqrt(b**2+a*sample[i])-b)/a             # inverse of primitive fctn of ax+b
plt.hist(sample, bins)
plt.savefig("1-triangulardistribution.png")
plt.clf()
'''

def addition(a, b):
    out = a
    for i in range(b):
        out += 1
    return out

def multiplication(a, b):
    out = 0
    for i in range(b):
        out = addition(out, a)
    return out

def exponent(a, b):
    out = 1
    for i in range(b):
        out = multiplication(out, a)
    return out

def hyperexponent(a, b):
    out = 1 # identitet för exponent
    for i in range(b):
        out = exponent(out, a)
    return out

test = hyperexponent(2, 2)
print(test)

#  if x_i>2/3, delta=2-x_i otherwise
# delta = np.sqrt(8*(a*x_i+b)*(1-x_i)/a)  # choose delta such that half is to be accepted

if x_i < delta / 2:
    # delta = (a*x_i**2+2*b*x_i-2*x_i+2)/(a*x_i+b)
    pass
