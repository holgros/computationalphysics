import numpy as np
from matplotlib import pyplot as plt

'''
TASK 1
'''
print('Starting task 1')
# input parameters - change as appropriate
r_min = 0.0 		            # lower range limit
r_max = 5.0 		            # upper range limit, should be larger than r_min
h = 0.01 			            # step, should be small compared to r_max-r_min
r = np.arange(r_min, r_max, h)  # array of r

# calculate discretization of u0
n = int((r_max-r_min)/h)        #number of steps
A = np.zeros( (n-2, n-2) )  # initialize A
for i in range(n-2):      # calculate A, see eq.43
  A[i][i] = 2/h
  if 0 < i:
    A[i][i-1] = -1/h
  if i < n-3:
    A[i][i+1] = -1/h
'''
f is calculated as follows:
f = u^2/r = 4*pi*r^2*phi(r)^2/r = 4*pi*r*phi(r)^2,
where phi(r) = 2*e^-x/sqrt(4*pi)
see eq.28 and 'useful quantities' eq.2
This together gives:
f = 4*r*e^-2r
'''
f = 4*r*np.exp(-2*r)
r = r[1:n - 1]  # remove first and last elements
def getHartree(f_in, r_vec):
    b = np.zeros(len(f_in) - 2)     # initialize b
    for i in range(len(f_in) - 2):  # calculate b, see eq.45
        b[i] = (h / 6) * (f_in[i] + 4 * f_in[i + 1] + f_in[i + 2])
    A_copy = A
    while len(A_copy[0]) > len(b):
        A_copy = np.delete(A_copy, slice(1), 0) # slice first row
        A_copy = np.delete(A_copy, slice(1), 1)  # slice first row
    A_inv = np.linalg.inv(A_copy)  # invert A
    c_out = np.dot(A_inv, b)  # calculate c=u0

    # calculate potential
    U_out = c_out + np.divide(r_vec, r_max)  # calculate U=u0+r/r_max
    V_out = np.divide(U_out, r_vec)  # calculate V=U/r
    return (c_out, U_out, V_out)
(c, U, Vh_numerical) = getHartree(f, r)
Vh_analytical = 1/r-(1+1/r)*np.exp(-2*r) # calculate V from eq.52

# plot solution to differential equation
plt.plot(r, U, label='U')
plt.plot(r, c, label='U_0')
plt.xlabel('r')
plt.suptitle('FEM: Task 1 - solution to differential equation')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.savefig("1-wavefunction.png")
plt.clf()

# plot Hartree potential
plt.plot(r, Vh_analytical, label='Hartree potential (analytical)')
plt.plot(r, Vh_numerical, 'ro', label='Numerical solution')
plt.xlabel('r')
plt.ylabel('V_H')
plt.suptitle('FEM: Task 1 - comparison between numerical and analytical solutions')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.savefig("1-potential.png")
plt.clf()

'''
TASK 3
'''
print('Starting task 3')
# generate matrix S
S = np.zeros( (n-2, n-2) )  # initialize S
for i in range(n-2):      # calculate S, see eq.50
  S[i][i] = 2*h/3
  if 0 < i:
    S[i][i-1] = h/6
  if i < n-4:
    S[i][i+1] = h/6

# generate matrix Q
# q(r_{i}) is defined for 2<i<n; Python uses zero-indexed arrays, so r[0]=r_{1} and so on
Q = np.zeros( (n-1, n-1) ) # initialize Q
r = np.arange(r_min, r_max, h)  # array of r
for i in range(1, n-1):     # calculate Q

  # limiting cases
  if r[i-1] == 0:
    Q[i][i] = 4 * r[i] / h + (2 * (r[i + 1] / h) ** 2) * np.log(r[i] / r[i + 1])
    Q[i][i - 1] = 1 - 2 * r[i] / h # equals -1 for r[i]=h
    Q[i][i + 1] = -1 - 2 * r[i] / h + (2 * r[i] * r[i + 1] / (h ** 2)) * np.log(r[i + 1] / r[i])
    continue

  # case i=j
  Q[i][i] = 4*r[i] / h - (2 * (r[i - 1] / h) ** 2) * np.log(r[i] / r[i - 1]) + (2 * (r[i + 1] / h) ** 2) * np.log(r[i] / r[i + 1])
  # See WolframAlpha:
  # integrate -2(x-(j-h))^2/xh^2 dx from j-h to j
  # integrate -2(j+h-x)^2/xh^2 dx from j to j+h

  # case i=j+1
  Q[i][i-1] = 1-2*r[i]/h+(2*r[i]*r[i-1]/(h**2))*np.log(r[i]/r[i-1])
  # See WolframAlpha: integrate -2(x-(j-h))*(j-x)/xh^2 dx from j-h to j

  # case i=j-1
  if i < n-2:
    Q[i][i+1] = -1 - 2 * r[i] / h + (2 * r[i] * r[i + 1] / (h ** 2)) * np.log(r[i+1] / r[i])
    # See WolframAlpha: integrate -2(x-j)*(j+h-x)/xh^2 dx from j to j+h

# q(r_{1}) is undefined, since r_{0} is undefined
Q = np.delete(Q, slice(1), 0) # slice first row
Q = np.delete(Q, slice(1), 1) # slice first column

# generate matrix H=S^{-1}(A+Q) and obtain its eigenvectors
A_plus_Q = A+Q                        # calculate matrix A+Q
S_inv = np.linalg.inv(S)            # invert S
H = np.dot(S_inv, A_plus_Q)           # calculate S^{-1}(A+Q)
lam, eig = np.linalg.eig(H)         # get eigenvalues and eigenvectors
# identify eigenvector with smallest eigenvalue
index_min = np.argmin(lam)
eig = np.transpose(eig)          # the eig function returns eigenvalues arranged column-wise
norm = 2*h*np.exp(-h)/eig[index_min][0]    # ... and normalized
eig[index_min] = norm*eig[index_min]

# print output
print('---TASK 3 OUTPUT---')
print("smallest eigenvalue:",lam[index_min])
print("ground state energy:",lam[index_min]/2)

# plot numerical and analytical solutions
psi = 2*r[1:n-1]*np.exp(-r[1:n-1])  # analytical solution
plt.plot(r[1:n-1], eig[index_min], 'ro', label='Calculated')
plt.plot(r[1:n-1], psi, 'b', label='Analytical solution 2re^-r')
plt.xlabel('r')
plt.ylabel('u')
#plt.suptitle('FEM: Task 3 - ground state wave function')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.savefig("3-wavefunction.png")
plt.clf()

'''
TASK 4
'''
print('Starting task 4')

# input from previous tasks
Q_0 = Q.copy()
u_0 = eig[index_min]
Vh = 1 / r[1:-1] - (1 + 1 / r[1:-1]) * np.exp(-2 * r[1:-1])
epsilon_old = lam[index_min]/2

# parameters
params = {'A': 0.0311,
          'B': -0.048,
          'C': 0.002,
          'D': -0.0116,
          'gamma': -0.1423,
          'beta_1': 1.0529,
          'beta_2': 0.3334}

# loop until eigenvalue converges
nbr_loops = 0
convergence = False                 # condition that abs(1-epsilon_new/epsilon_old) < conv_ratio
conv_ratio = 0.000001               # 10^-6 requires about 5 loops for h=0.01
while not convergence:
    # update wavefunction with limit values
    u_0 = np.insert(u_0, 0, 0.0)            # insert 0.0 at beginning
    u_0 = np.append(u_0, 0.0)               # insert 0.0 at end

    '''
    # get Vh
    b = np.zeros(len(u_0))                  # initialize b
    for i in range(len(b)-2):               # calculate b, see eq.45
        b[i] = (h / 6) * (u_0[i] + 4 * u_0[i + 1] + u_0[i + 2])
    # append approximations at the end where u_0 is approximately zero
    b[-2] = (h / 6) * (u_0[-2] + 4 * u_0[-1])
    b[-1] = (h / 6) * u_0[-1]
    A_inv = np.linalg.inv(A)                # invert A
    Vh = np.dot(A_inv, b[1:-1])             # solve Poisson equation
    Vh = Vh + np.divide(r[1:-1], r_max)     # calculate corresponding wavefunction for single orbital
    Vh = np.divide(Vh, r[1:-1])             # get Hartree potential
    '''

    # get densities
    n_densities = []
    for i in range(1, len(u_0)-1):
        n_densities.append((np.abs(u_0[i]/r[i])**2)/(2*np.pi))

    # get Vx
    Vx = []
    for i in range(0, len(n_densities)):
        Vx.append(-(3*n_densities[i]/np.pi)**(1/3))

    # get Vc
    Vc = []
    for i in range(0, len(n_densities)):
        if n_densities[i] > 3 / (4 * np.pi):
            term1 = (params['A'] / 3) * np.log(3 / (4 * np.pi * n_densities[i]))
            term2 = params['B']
            term3 = -params['A'] / 3
            term4 = (2 / 9) * params['C'] * (3 / 4 * np.pi * n_densities[i]) ** (1 / 3) * np.log(3 / 4 * np.pi * n_densities[i])
            term5 = ((2 * params['D'] - params['C']) / 3) * (3 / (4 * np.pi * n_densities[i])) ** (1 / 3)
            Vc.append(term1 + term2 + term3 + term4 + term5)
        else:
            nominator = 1 + (7 / 6) * params['beta_1'] * (3 / (4 * np.pi * n_densities[i])) ** (1 / 6) + params['beta_2'] * (
                        3 / (4 * np.pi * n_densities[i])) ** (1 / 3)
            denominator = (1 + params['beta_1'] * (3 / 4 * np.pi * n_densities[i]) ** (1 / 6) + params['beta_2'] * (
                        3 / (4 * np.pi * n_densities[i])) ** (1 / 3)) ** 2
            Vc.append(params['gamma'] * nominator / denominator)

    # calculate matrix Q
    Q = Q_0.copy()      # first term copied from task 3
    V_tot = Vh+Vx+Vc    # total potential
    for i in range(len(V_tot)):
        for j in range(len(V_tot)):
            # take care of two limit cases by assigning V_tot[-1]=V_tot[0] and V_tot[r_max+1]=V_tot[r_max]
            if i == 0 and j == 0:
                Q[i][j] += (h / 6) * (V_tot[0] + 2 * V_tot[0] + V_tot[1])
            elif i == len(V_tot)-1 and j == len(V_tot)-1:
                Q[i][j] += (h / 6) * (V_tot[len(V_tot)-2] + 2 * V_tot[len(V_tot)-1] + V_tot[len(V_tot)-1])
            # remaining cases as in eq.51
            elif i == j:
                Q[i][j] += (h/6)*(V_tot[i-1]+2*V_tot[i]+V_tot[i+1])
            elif np.abs(i-j) == 1:
                Q[i][j] += (h/12)*(V_tot[i]+V_tot[j])
    Q = 2*Q             # note multiplication with 2

    # repeat procedure from task 3
    A_plus_Q = A+Q
    S_inv = np.linalg.inv(S)
    H = np.dot(S_inv, A_plus_Q)          # calculate S^{-1}(A+Q)
    eigvals, eigvecs = np.linalg.eig(H)  # get eigenvalues and eigenvectors
    eigvecs = np.transpose(eigvecs)      # the eig function returns eigenvalues arranged column-wise
    '''
    # normalize wavefunctions
    for i in range(len(eigvecs)):
        sum = 0
        for y in eigvecs[i]:
            term = h * y * y
            sum += term
        eigvecs[i] /= np.sqrt(sum)
    '''
    # get index of minimum energy and reset u_0 as the corresponding wavefunction
    index_min = np.argmin(eigvals)
    u_0 = eigvecs[index_min]
    epsilon_new = eigvals[index_min]/2

    # print to debug and track process
    nbr_loops += 1
    print('After', nbr_loops, 'loops:\t\tEPSILON =', epsilon_new)
    #print('u_0:',np.ndarray.tolist(u_0))

    # check convergence condition
    if np.abs(1 - epsilon_new / epsilon_old) < conv_ratio:
        convergence = True
    epsilon_old = epsilon_new

# update copy of wavefunction with limit values
wavefctn = u_0.copy()
wavefctn = np.insert(wavefctn, 0, 0.0)  # insert 0.0 at beginning
wavefctn = np.append(wavefctn, 0.0)     # ... and 0.0 at end
wavefctn = np.abs(wavefctn)

'''
# normalize wavefunction
sum = 0
for y in wavefctn:
    term = h * y * y
    sum += term
wavefctn /= np.sqrt(sum)
#print('sum:',sum)
'''

# plot wavefunction
plt.plot(r, wavefctn, 'ro', markersize=0.2, label='Wavefunction')
plt.xlabel('r')
plt.ylabel('u')
plt.suptitle('FEM: Task 4 - helium ground state wavefunction')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
plt.savefig("4-wavefunction.png")
plt.clf()

# calculate ground state energy
energy = eigvals[index_min]          # first term is 2*epsilon = eigenvalue of the solution above
# get densities - same as above
n_densities = []
for j in range(0, len(wavefctn)-1):
    n_densities.append((np.abs(wavefctn[j] / r[j+1]) ** 2) / (2 * np.pi))
n_densities.append(0.0)
# add the Vh term to the energy
Vh_term = 0
'''
# get Vh
b = np.zeros(len(wavefctn))                  # initialize b
for i in range(len(b)-2):               # calculate b, see eq.45
    b[i] = (h / 6) * (wavefctn[i] + 4 * wavefctn[i + 1] + wavefctn[i + 2])
# append approximations at the end where u_0 is approximately zero
b[-2] = (h / 6) * (wavefctn[-2] + 4 * wavefctn[-1])
b[-1] = (h / 6) * wavefctn[-1]
A_inv = np.linalg.inv(A)                # invert A
Vh = np.dot(A_inv, b[1:-1])             # solve Poisson equation
Vh = Vh + np.divide(r[1:-1], r_max)     # calculate corresponding wavefunction for single orbital
Vh = np.divide(Vh, r[1:-1])             # get Hartree potential
'''
for j in range(len(n_densities[1:-1])):
    term = h * n_densities[j] * Vh[j]
    Vh_term += term
energy -= 0.5*Vh_term          # note minus sign!
# get the Exc term and add to the energy
Exc_term = 0.0
for n in n_densities:  # calculating integral of n*epsilon_x
    term = h * n * (-3 / 4) * ((3 / np.pi) ** (1 / 3)) * n ** (1 / 3)
    Exc_term += term
for n in n_densities:  # calculating integral of n*epsilon_c
    epsilon_c = 0.0
    if n > 3 / (4 * np.pi):
        term1 = (params['A'] / 3) * np.log(3 / (4 * np.pi * n))
        term2 = params['B']
        term3 = (1 / 3) * params['C'] * (3 / 4 * np.pi * n) ** (1 / 3) * np.log(3 / 4 * np.pi * n)
        term4 = params['D'] * (3 / (4 * np.pi * n)) ** (1 / 3)
        epsilon_c = term1 + term2 + term3 + term4
    else:
        epsilon_c = params['gamma'] / (1 + params['beta_1'] * (3 / 4 * np.pi * n) ** (1 / 6) + params['beta_2'] * (
                    3 / 4 * np.pi * n) ** (1 / 3))
    term = h * n * epsilon_c
    Exc_term += term
energy += Exc_term
# get the Vxc energy term and add to the energy
# get Vx - same as above
Vx = []
for j in range(0, len(n_densities)):
    Vx.append(-(3 * n_densities[j] / np.pi) ** (1 / 3))
# get Vc - same as above
Vc = []
for j in range(0, len(n_densities)):
    if n_densities[j] > 3 / (4 * np.pi):
        term1 = (params['A'] / 3) * np.log(3 / (4 * np.pi * n_densities[j]))
        term2 = params['B']
        term3 = -params['A'] / 3
        term4 = (2 / 9) * params['C'] * (3 / 4 * np.pi * n_densities[j]) ** (1 / 3) * np.log(
            3 / 4 * np.pi * n_densities[j])
        term5 = ((2 * params['D'] - params['C']) / 3) * (3 / (4 * np.pi * n_densities[j])) ** (1 / 3)
        Vc.append(term1 + term2 + term3 + term4 + term5)
    else:
        if n_densities[j] == 0:         # limiting case
            Vc.append(0.0)
            continue
        nominator = 1 + (7 / 6) * params['beta_1'] * (3 / (4 * np.pi * n_densities[j])) ** (1 / 6) + params[
            'beta_2'] * (
                            3 / (4 * np.pi * n_densities[j])) ** (1 / 3)
        denominator = (1 + params['beta_1'] * (3 / 4 * np.pi * n_densities[j]) ** (1 / 6) + params['beta_2'] * (
                3 / (4 * np.pi * n_densities[j])) ** (1 / 3)) ** 2
        Vc.append(params['gamma'] * nominator / denominator)
# calculate the Vxc term
Vxc = Vx+Vc
Vxc_term = 0.0
for j in range(len(n_densities)):
    term = h * n_densities[j] * Vxc[j]
    Vxc_term -= term  # note minus sign!
energy += Vxc_term
print('---TASK 4 OUTPUT---')
print('Ground state energy of Helium:',energy)
