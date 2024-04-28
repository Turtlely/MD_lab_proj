
# Constants
ep = 1
sig = 1
alpha = 0.01
convergence_criteria=0.01

def force(r):
    return -4*ep*(-12*(sig/r)**13+6*(sig/r)**7)

def gradient_descent(r,alpha,convergence_criteria):
    n=0
    while(force(r)>convergence_criteria):
        r=r+alpha*force(r)
        n+=1
        print(n)
    return r

print(gradient_descent(1.05,alpha,convergence_criteria))

# Converges to r=1.1223