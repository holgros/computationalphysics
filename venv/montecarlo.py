import numpy as np
import random
from matplotlib import pyplot as plt

test_mode = True

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
print('---TASK 2 OUTPUT---')
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

def metropolis_task1(N, a=1.0, b=0.5):
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
            sample.append(x_i)                                      # NOT SURE??
            continue
        r = random.uniform(0,1)
        if r < (a*x_j+b)/(a*x_i+b): # note that this entails x_j > x_i
            x_i = x_j
        sample.append(x_i)                                          # NOT SURE ABOUT ADDING OLD VALUE AGAIN??
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
metropolis_task1(10)
metropolis_task1(100)
metropolis_task1(1000)
metropolis_task1(10000)

'''
bins = [0.1 * i for i in range(12)]
sample = [random.uniform(0,2/a+b) for i in range(25000)]    # not quite sure about 2/a+b, got there by trial and error
for i in range(len(sample)):
    sample[i] = (np.sqrt(b**2+a*sample[i])-b)/a             # inverse of primitive fctn of ax+b
plt.hist(sample, bins)
plt.savefig("1-triangulardistribution.png")
plt.clf()
'''

'''
TASK 2
'''
print('Starting task 2')

# the wave function - input consists of two vectors in cartesian coordinates and a value of alpha
def psi(alpha, vec_1, vec_2=None):
    if vec_2 is None:
        vec_2 = vec_1[1]
        vec_1 = vec_1[0]
    r_1 = np.linalg.norm(vec_1)
    r_2 = np.linalg.norm(vec_2)
    r_12 = np.linalg.norm(vec_2-vec_1)
    return np.exp(-2*r_1)*np.exp(-2*r_2)*np.exp(r_12/(2*(1+alpha*r_12)))

# get samples for monte carlo integration
def get_samples(alpha, N=100000, N_eq=10000):
    first_sample = [np.zeros(3), np.zeros(3)]
    for i in range(6):
        first_sample[i%2][i%3] = random.uniform(-xyz_max, xyz_max)
    #x_i = np.array([random.uniform(-xyz_max, xyz_max),random.uniform(-xyz_max, xyz_max),random.uniform(-xyz_max, xyz_max)])
    x_i = sample(alpha, first_sample, N_eq)[-1]       # sample N_eq times first
    return sample(alpha, x_i, N)
def sample(alpha, starting_points, nbr_samples):
    out = []
    dim = 0  # x, y or z dimension
    out.append(starting_points)
    while len(out) < nbr_samples:
        r = random.uniform(-xyz_max, xyz_max)
        next_coordinate = starting_points[dim%2][dim%3] + delta * (r - 0.5)
        if next_coordinate < -xyz_max or next_coordinate > xyz_max:
            out.append(starting_points)                  # NOT SURE??
            continue
        # create two new points where one coordinate of one point differs from starting points
        next_points = [np.zeros(3), np.zeros(3)]
        for i in range(6):
            next_points[i % 2][i % 3] = starting_points[i % 2][i % 3]
        next_points[dim%2][dim%3] = next_coordinate
        # evaluate probability and decide whether to accept the next points
        prob = (psi(alpha,next_points)/psi(alpha,starting_points))**2
        r = random.uniform(0, 1)
        if r < prob:  # note that this entails 1 < prob
            starting_points = next_points
        out.append(starting_points)                      # NOT SURE ABOUT ACCEPTING OLD VALUE AGAIN??
        # alternate between x, y and z dimensions for two points
        dim += 1
        dim %= 6
    return out

# local energy - input same as for wavefunction
def local_energy(alpha, positions):
    vec_1 = positions[0]
    vec_2 = positions[1]
    r_1 = np.linalg.norm(vec_1)
    r_2 = np.linalg.norm(vec_2)
    r_12 = np.linalg.norm(vec_2-vec_1)
    r_1_hat = vec_1/r_1
    r_2_hat = vec_2/r_2
    r_12_hat = (vec_2-vec_1)/r_12
    term_1 = -4
    term_2 = np.dot(r_1_hat-r_2_hat, r_12_hat)/(1+alpha*r_12)**2
    term_3 = alpha/(1+alpha*r_12)
    term_4 = alpha/(1+alpha*r_12)**2
    term_5 = alpha/(1+alpha*r_12)**3
    term_6 = -1/(4*(1+alpha*r_12)**4)
    return term_1+term_2+term_3+term_4+term_5+term_6

def metropolis_task2(alpha, sample):
    N = len(sample)
    # calculate integral
    psi_2 = []
    E_L = []
    for i in range(N):
        psi_2.append(psi(alpha, sample[i])*psi(alpha, sample[i]))
        E_L.append(local_energy(alpha, sample[i]))
    psi_2 = np.asarray(psi_2)
    normalization_factor = np.sum(psi_2)/N
    rho = psi_2/normalization_factor             # NOT SURE HOW THIS AFFECTS THE CORRELATION??
    E_L = np.asarray(E_L)
    I_N = np.sum(E_L*rho)/N
    # calculate autocorrelation
    psi_array = []      # append to output
    k = -1
    psi_k = 1           # arbitrary value larger than 0.1
    while np.abs(psi_k) > 0.1 and k < N - 1:  # np.abs really needed??
        k += 1
        corr = []
        # g = f/w = (E_L*rho)/rho = E_L
        g = E_L
        x = sample
        # same as above
        for i in range(len(x) - k):
            corr.append(g[i] * g[i + k])
        corr_mean = np.mean(corr)
        g_mean = np.mean(g)
        g_square_mean = np.mean(g * g)
        psi_k = (corr_mean - g_mean * g_mean) / (g_square_mean - g_mean * g_mean)
        psi_array.append(psi_k)       # append to output
    # print results
    return I_N, psi_array

# main process
delta = 0.5                     # NOT SURE??
xyz_max = 5                     # NOT SURE??
alphas = (0.100, 0.125, 0.150)
N=100000
N_eq=10000
M = 10
if test_mode:
    N=1000
    M = 1

# perform calculations M times for each value of alpha
output = []
for i in range(len(alphas)):
    print('Working with alpha=',alphas[i])      # track progress
    all_samples = []
    integral_results = []
    correlation_results = []
    for j in range(M):
        sample_array = get_samples(alphas[i], N, N_eq)
        all_samples.append(sample_array)
        I_N, psi_array = metropolis_task2(alphas[i], sample_array)
        integral_results.append(I_N)
        correlation_results.append(psi_array)
        print('\tCompleted M=',j)
    output.append((all_samples, integral_results, correlation_results))

print('---TASK 2 OUTPUT---')
print('Mean values for I_N:')
for i in range(len(alphas)):
    print('alpha=', alphas[i], ':\t\tI_N=',np.mean(output[i][1]))

'''
TASK 3
'''
print('Starting task 3')

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
# TODO: calculate g=n_2/(max_freq*cube_2)


print(max_pos)
print(max_freq)


'''
for i in range(N):
    #print(r_1[i])
    for x in range(nbr_steps):
        if not (x == nbr_steps-1 or (steps[x] < r_1[i][0] and steps[x+1] >= r_1[i][0])):
            continue
        for y in range(nbr_steps):
            if not (y == nbr_steps - 1 or (steps[y] < r_1[i][1] and steps[y + 1] >= r_1[i][1])):
                continue
            for z in range(nbr_steps):
                if not (z == nbr_steps - 1 or (steps[z] < r_1[i][2] and steps[z + 1] >= r_1[i][2])):
                    continue
                cube_1[x][y][z] += 1
                #print(x,' ',y,' ',z)
                break
            break
        break
cube_2 = np.zeros((nbr_steps,nbr_steps,nbr_steps))
for i in range(N):
    #print(r_2[i])
    for x in range(nbr_steps):
        if not (x == nbr_steps-1 or (steps[x] < r_2[i][0] and steps[x+1] >= r_2[i][0])):
            continue
        for y in range(nbr_steps):
            if not (y == nbr_steps - 1 or (steps[y] < r_2[i][1] and steps[y + 1] >= r_2[i][1])):
                continue
            for z in range(nbr_steps):
                if not (z == nbr_steps - 1 or (steps[z] < r_2[i][2] and steps[z + 1] >= r_2[i][2])):
                    continue
                cube_2[x][y][z] += 1
                #print(x,' ',y,' ',z)
                break
            break
        break
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
