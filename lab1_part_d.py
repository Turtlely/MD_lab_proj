# Constants
ep = 1
sig = 1
convergence_criteria=0.01

def force(r):
    return -4*ep*(-12*(sig/r)**13+6*(sig/r)**7)

def second_deriv(r):
    return 4*ep*(12*13*(sig/r)**14-6*7*(sig/r)**8)

def newtons_method(r,convergence_criteria):
    n=0
    while(abs(force(r))>convergence_criteria):
        n+=1
        print(n,force(r))
        r=r+force(r)/second_deriv(r)
    return r

print(newtons_method(0.95,convergence_criteria))

# Converges to r=1.1225 with r0=0.95
# Converges to r=1.1224 with r0=1.2
# Converges to r=3.0619 with r0=1.3. This fails because the algorithm does not converge to a minima, it just converges to the asymptotically decreasing force at infinity
