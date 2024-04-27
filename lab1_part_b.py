"""
Part A: Damped molecular dynamics

This is a practice that introduces you to iterative algorithms.
To get started, inside your directory for lab1, create a new script named lab1_part_b.py.
add the shebang to the topline of this script.
The shebang is #!/usr/bin/env python.
Save and exit the script. Then make the script an executable with the command chmod like this: chmod +x lab1_part_b.py
This algorithm will simulate two atoms separated at a certain distance and find the optimum bond length. 
"""

# Constants
r0 = 0.3
v0 = 0
force_thresh = 0.01
m = 1
ep = 1
sig = 1
dt = 0.01
dx_thresh = 0.1

# LJ potential
def V_lj(r):
    return 4*ep*((sig/r)**12-(sig/r)**6)

# LJ force
def F_lj(r):
    return -4*ep*(-12*(sig/r)**13+6*(sig/r)**7)

# F=ma
def a(F,m):
    return F/m

# Move
def move(v,a):
    # Return dx and dv
    return v*dt + 1/2 * a * dt**2, a*dt

# Damper for convergence on potential minima
def damp(dx,v):
    if dx > dx_thresh:
        return dx_thresh, 0
    else:
        return dx, 0.9*v

# Main loop
def main():
    # Center of mass will be at the origin x=0, 
    # so the initial position of the atoms is r0/2

    # Set initial position, velocity, and force
    r = r0
    x = r0/2
    v = 0
    F = F_lj(r0)

    # Optimization epoch counter
    n=0
    # Iterate until minima in force is found
    while abs(F)>force_thresh:
        n+=1
        # Calculate force between particles
        F = F_lj(r)
        # Calculate resulting acceleration
        ax = a(F,m)
        # Calculate change in position and velocity
        dx, dv = move(v,ax)
        # Update velocity
        v+=dv
        # Damp movements
        dx, v = damp(dx,v)
        # Update positions using a sympletic euler iterative algorithm
        # Positions are updated using the new velocity
        x += dx
        # Due to the symmetry of the problem, 
        # we can constrain the distance between the particles to be twice their x position 
        # as the CM position must be stationary.
        r = 2*x

        print(f"Iteration {n}. F = {F}")
    
    print("Simulation complete: ")
    print(f"Final r: {r}")
    print(f"Final F: {F}")

if __name__ == "__main__":
    main()

# Results:
# Convergence to r=1.12 within 210 iterations, not bad! 