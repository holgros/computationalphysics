import numpy as np
import random
from matplotlib import pyplot as plt

test_mode = False

'''
TASK 1
'''
# main process begins after function definitions

def uniform(N):
    sample = [random.uniform(0, 1) for i in range(N)]
    x = np.asarray(sample)
    f = 1 / (1 + x ** 3)
    f2 = (1 / (1 + x ** 3)) ** 2
    I_N = np.sum(f) / N
    variance = np.sum(f2) / N - I_N ** 2
    sigma = np.sqrt(variance)
    error = sigma/np.sqrt(N)
    print('N=', N, ': \tI_N=', I_N, ', error=', error)

def linear(N, a=-0.5, b=1.25):
    #sample = [random.uniform(0, np.sqrt((2*a+b)/b)-b) for i in range(N)]    # not quite sure about np.sqrt((2*a+b)/b)-b, got there by trial and error
    sample = [random.uniform(0, 1) for i in range(N)]
    for i in range(len(sample)):
        # xs=-b/a - sqrt(2/a*rand(1,N) + (b/a)^2 )
        #sample[i] = (np.sqrt(b ** 2 + a * sample[i]) - b) / a               # inverse of primitive fctn of ax+b
        sample[i] = -b/a-np.sqrt(2/a*sample[i]+(b/a)**2)
    x = np.asarray(sample)
    f = 1 / (1 + x ** 3)
    f = f/(a*x+b)
    f2 = f*f
    I_N = np.sum(f) / N
    variance = np.sum(f2) / N - I_N ** 2
    sigma = np.sqrt(variance)
    error = sigma/np.sqrt(N)
    print('N=', N, ': \tI_N=', I_N, ', error=', error)

def metropolis_task1(N, a=-0.5, b=1.25):
    sample = []
    x_i = random.uniform(0, 1)  # initial state
    sample.append(x_i)
    nbr_accepts = 0
    nbr_rejects = 0
    delta = 1.9
    # populate sample
    while len(sample) < N:
        nbr_rejects += 1
        '''
        delta = np.sqrt(8*(a*x_i+b)*(1-x_i)/a)              # choose delta such that half is to be accepted
        if x_i < delta/2:
            delta = (a*x_i**2+2*b*x_i-2*x_i+2)/(a*x_i+b)    # see above
        '''
        # generate trial change
        r = random.uniform(0,1)
        x_j = x_i+delta*(r-0.5)
        # determine whether x_j is to be accepted
        if x_j < 0 or x_j > 1:
            sample.append(x_i)                                      # NOT SURE??
            continue
        r = random.uniform(0,1)
        if r < (a*x_j+b)/(a*x_i+b): # note that this entails x_j < x_i
            x_i = x_j
            nbr_rejects -= 1
            nbr_accepts += 1
        sample.append(x_i)                                          # NOT SURE ABOUT ADDING OLD VALUE AGAIN??
    # calculate integral
    x = np.asarray(sample)
    f = 1 / (1 + x ** 3)
    f = f/(a*x+b)
    I_N = np.sum(f) / N
    # calculate standard error
    f2 = f*f
    variance = np.sum(f2) / N - I_N ** 2
    sigma = np.sqrt(variance)
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
    error = sigma/np.sqrt(N/k)
    # print results
    acceptance_rate = nbr_accepts/(nbr_accepts+nbr_rejects)
    print('N=', N, ': \tI_N=', I_N, '\terror=', error, '\tpsi_k=', psi_k, 'for k=', k, '\tacceptance_rate=',acceptance_rate)
    #print('nbr_accepts:',nbr_accepts,'\tnbr_rejects:',nbr_rejects)

# MAIN PROCESS
if __name__ == '__main__':
    print('Starting task 1')
    print('---TASK 1 OUTPUT---')
    print('--UNIFORM DISTRIBUTION--')
    uniform(10)
    uniform(100)
    uniform(1000)
    uniform(10000)
    print('--LINEAR DISTRIBUTION--')
    linear(10)
    linear(100)
    linear(1000)
    linear(10000)
    print('--METROPOLIS--')
    metropolis_task1(10)
    metropolis_task1(100)
    metropolis_task1(1000)
    metropolis_task1(10000)

'''
TASK 2
'''
print('Starting task 2')
# main process begins after function definitions

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
    samples, acceptance_rate = sample(alpha, first_sample, N_eq)       # sample N_eq times first
    x_i = samples[-1]
    return sample(alpha, x_i, N)
def sample(alpha, starting_points, nbr_samples):
    # alpha = corresponding parameter in trial wave function
    # starting_points = a pair of coordinate triples ((x1,y1,z1),(x2,y2,z2)) corresponding to the initial positions
    # nbr_samples = number of items to be returned; each item is a pair of coordinate triples
    out = []
    dim = 0  # x, y or z dimension
    out.append(starting_points)
    nbr_accepted = 0
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
            nbr_accepted += 1
        out.append(starting_points)                      # NOT SURE ABOUT ACCEPTING OLD VALUE AGAIN??
        # alternate between x, y and z dimensions for two points
        dim += 1
        dim %= 6
    acceptance_rate = nbr_accepted/nbr_samples
    return out, acceptance_rate

# local energy - input same as for wavefunction
def local_energy(alpha, positions):
    vec_1 = positions[0]
    vec_2 = positions[1]
    r_1 = np.linalg.norm(vec_1)
    r_2 = np.linalg.norm(vec_2)
    r_12 = np.linalg.norm(vec_1-vec_2)
    r_1_hat = vec_1/r_1
    r_2_hat = vec_2/r_2
    r_12_hat = (vec_1-vec_2)/r_12
    term_1 = -4
    term_2 = np.dot(r_1_hat-r_2_hat, r_12_hat)/((1+alpha*r_12)**2)
    term_3 = alpha/(1+alpha*r_12)
    term_4 = alpha/((1+alpha*r_12)**2)
    term_5 = alpha/((1+alpha*r_12)**3)
    term_6 = -1/(4*((1+alpha*r_12)**4))
    return term_1+term_2+term_3+term_4+term_5+term_6

def calculate_integral(alpha, sample):
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
    #I_N = np.sum(E_L*rho)/N                    # SHOULD DIVIDE INTEGRAND WITH WEIGHT FUNCTION??
    for i in range(len(E_L)):
        #print(E_L[i])
        pass
    I_N = np.sum(E_L) / N
    # calculate autocorrelation
    psi_array = []      # append to output
    k = -1
    psi_k = 1           # arbitrary value larger than 0.1
    while np.abs(psi_k) > 0.1 and k < N - 1:  # np.abs really needed??
        k += 1
        corr = []
        # g = f/w = ((E_L*rho)/rho)/rho = E_L/rho
        g = E_L/rho
        x = sample
        # same as above
        for i in range(len(x) - k):
            corr.append(g[i] * g[i + k])
        corr_mean = np.mean(corr)
        g_mean = np.mean(g)
        g_square_mean = np.mean(g * g)
        psi_k = (corr_mean - g_mean * g_mean) / (g_square_mean - g_mean * g_mean)
        psi_array.append(psi_k)       # append to output
    # calculate variance
    variance = np.sum(E_L*E_L) / N - I_N ** 2
    sigma = np.sqrt(variance)
    error = sigma/np.sqrt(N/len(psi_array))
    # return results
    return I_N, psi_array, E_L

# MAIN PROCESS
if __name__ == '__main__':
    # main process
    delta = 0.3                     # got this value by trial and error, but doesn't appear to affect the result much
    xyz_max = 5                     # at this distance from the origin, probabilities are very small
    alphas = (0.100, 0.125, 0.150)
    N=100000
    N_eq=10000
    M = 10
    if test_mode:
        M = 1
        alphas = [0.125]

    # perform calculations M times for each value of alpha
    output = []
    for i in range(len(alphas)):
        print('Working with alpha=',alphas[i])      # track progress
        all_samples = []
        integral_results = []
        k = []
        local_energies = []
        errors = []
        for j in range(M):
            sample_array, acceptance_rate = get_samples(alphas[i], N, N_eq)
            all_samples.append(sample_array)
            I_N, psi_array, E_L = calculate_integral(alphas[i], sample_array)
            integral_results.append(I_N)
            local_energies.append(E_L)
            k.append(len(psi_array))
            variance = np.sum(E_L * E_L) / N - np.mean(E_L) ** 2
            sigma = np.sqrt(variance)
            errors.append(sigma / np.sqrt(N / len(psi_array)))
            print('\tCompleted M=',j,'\tacceptance_rate:', acceptance_rate)
            if (test_mode):
                for i in range(len(sample_array)):
                    #print(sample_array[i])
                    pass
        '''
        # calculate error
        local_energies = np.asarray(local_energies)
        variance = np.sum(local_energies * local_energies) / (N*M) - np.mean(local_energies) ** 2
        sigma = np.sqrt(variance)
        error = sigma / np.sqrt(N*M / np.sum(k))          # NOT SURE ABOUT THE LAST DENOMINATOR??
        '''
        # output results for this value of alpha
        output.append((all_samples, integral_results, k, errors))

    print('---TASK 2 OUTPUT---')
    print('Mean values for I_N:')
    for i in range(len(alphas)):
        print('alpha=', alphas[i], ':\t\tI_N=',np.mean(output[i][1]),'\t\terror=',np.sum(output[i][3]))

'''
TASK 3

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
