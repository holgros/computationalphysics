
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
#import montecarlo as mc

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
'''
if x_i < delta / 2:
    # delta = (a*x_i**2+2*b*x_i-2*x_i+2)/(a*x_i+b)
    pass
'''
#x = mc.psi(0.125, np.asarray([0,0,3]),np.asarray([0,0,0]))
#print(x)
import numpy as np
array = np.random.random((5,5))
print(array)

density = np.histogram(array, density=True)
print(density)


'''
print('len(h):', len(h))
print('h[steps//2][steps//2]:',h[steps//2][steps//2])
max_ycol = np.argmax(sum(h))
max_xcol = np.argmax(h[max_ycol])
print('max:',np.max(h[max_ycol]))
print('max:',h[max_ycol][max_xcol])
print('coords x-y:',max_xcol, max_ycol)
'''

# Big bins
#plt.hist2d(x, y, bins=(50, 50), cmap=plt.cm.jet)
# plt.show()

# Small bins

# If you do not set the same values for X and Y, the bins aren't square !
#plt.hist2d(x, y, bins=(300, 30), cmap=plt.cm.jet)

'''
# parameter to set discretization steps for the distribution functions
h = 0.1
if test_mode:
    h=1
steps = np.arange(-xyz_max, xyz_max, h)
# from output of task 2, choose sample for some intermediate value of alpha
alpha = len(alphas)//2
X = np.asarray(output[alpha][0][0])                         # NOTE EXTRA DIMENSION NEEDED HERE??
# calculate rho as above
psi_2 = []
for i in range(N):
    psi_2.append(psi(alpha, X[i])*psi(alpha, X[i]))
psi_2 = np.asarray(psi_2)
normalization_factor = np.sum(psi_2)/N
rho = psi_2/normalization_factor                            # NOT SURE HOW THIS AFFECTS THE CORRELATION??
# get distributions as one cube for each electron
r_1 = []
r_2 = []
for i in range(N):
    r_1.append(X[i][0])
    r_2.append(X[i][1])
nbr_steps = len(steps)
cube_1 = np.zeros((nbr_steps,nbr_steps,nbr_steps))
cube_2 = np.zeros((nbr_steps,nbr_steps,nbr_steps))
distributions = ((r_1,cube_1),(r_2,cube_2))
for distr in distributions:
    r = distr[0]
    cube = distr[1]
    for i in range(N):
        #print(r[i])
        for x in range(nbr_steps):
            if not (x == nbr_steps-1 or (steps[x] < r[i][0] and steps[x+1] >= r[i][0])):
                continue
            for y in range(nbr_steps):
                if not (y == nbr_steps - 1 or (steps[y] < r[i][1] and steps[y + 1] >= r[i][1])):
                    continue
                for z in range(nbr_steps):
                    if not (z == nbr_steps - 1 or (steps[z] < r[i][2] and steps[z + 1] >= r[i][2])):
                        continue
                    cube[x][y][z] += 1
                    #print(x,' ',y,' ',z)
                    break
                break
            break
# identify most frequent position for electron 1
max_freq = -1
max_pos = ()
for x in range(nbr_steps):
    for y in range(nbr_steps):
        for z in range(nbr_steps):
            if cube_1[x][y][z] > max_freq:
                max_freq = cube_1[x][y][z]
                max_pos = (x,y,z)
# build n_2 as a cube of conditional probabilities where the coordinates of the first electron is held fixed at max_pos
conditional_r_2 = []            # identify all r_2 for which r_1 = max_pos
for i in range(N):
    if (steps[max_pos[0]] < r_1[i][0] and steps[max_pos[0] + 1] >= r_1[i][0] and
            steps[max_pos[1]] < r_1[i][1] and steps[max_pos[1] + 1] >= r_1[i][1] and
            steps[max_pos[2]] < r_1[i][2] and steps[max_pos[2] + 1] >= r_1[i][2]):
            conditional_r_2.append(r_2[i])
n_2 = np.zeros((nbr_steps,nbr_steps,nbr_steps))
for i in range(len(conditional_r_2)):
    for x in range(nbr_steps):
        if not (x == nbr_steps - 1 or (steps[x] < conditional_r_2[i][0] and steps[x + 1] >= conditional_r_2[i][0])):
            continue
        for y in range(nbr_steps):
            if not (y == nbr_steps - 1 or (steps[y] < conditional_r_2[i][1] and steps[y + 1] >= conditional_r_2[i][1])):
                continue
            for z in range(nbr_steps):
                if not (z == nbr_steps - 1 or (steps[z] < conditional_r_2[i][2] and steps[z + 1] >= conditional_r_2[i][2])):
                    continue
                n_2[x][y][z] += 1
                # print(x,' ',y,' ',z)
                break
            break
        break
# calculate g
g = np.zeros((nbr_steps,nbr_steps,nbr_steps))
for x in range(nbr_steps):
    for y in range(nbr_steps):
        for z in range(nbr_steps):
            if cube_2[x][y][z] != 0:
                g[x][y][z] = n_2[x][y][z]/(max_freq*cube_2[x][y][z])
            else:
                g[x][y][z] = 0
# roll-up the z dimension for visualization
g_twodim = np.zeros((nbr_steps,nbr_steps))
for x in range(nbr_steps):
    for y in range(nbr_steps):
        for z in range(nbr_steps):
            g_twodim[x][y] += g[x][y][z]
# visualize g
from mpl_toolkits import mplot3d
x = np.outer(np.linspace(-xyz_max, xyz_max, int(2*xyz_max/h)), np.ones(int(2*xyz_max/h)))
y = x.copy().T # transpose
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(x, y, g_twodim,cmap='viridis', edgecolor='none')
ax.set_title('Helium atom distribution')
plt.savefig("3-Helium-atom-distribution.png")
'''
