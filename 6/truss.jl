# Used for math
using LinearAlgebra

# Used to output CSV files
using DataFrames
using CSV

# ==========
# Input parameters
m = 16; # Number of members in the system
n = 9; # Number of nodes
dof = 2*n; # Degrees of freedom, two directions for each nodes
a = 2; # Half length of system
h1 = 4; # Height of second level
h2 = 2; # Height of third level relative to second
h3 = h1+h2; # Total height of third level

# Matrix of node coordinates
# For each node, gives x and y coordinates
#    1  2   3    4    5     6      7   8       9
x = [0, a, 2*a, a/3,  a,   10/3, a/2,  a,  a+(a/2)];
y = [0, 0,   0,  h1, h1,      h1,  h3, h3,       h3];
# Changing x of 6 to incorrect value from homework solution, should be (a/3)+a

# Matrix that gives the beam connections, each item is a node value
C = [[1, 2, 1, 2, 2, 2, 3, 4, 5, 4, 5, 5, 5, 6, 7, 8] [2, 3, 4, 4, 5, 6, 6, 5, 6, 7, 7, 8, 9, 9, 8, 9]];
#=
  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16
[[1, 2, 1, 2, 2, 2, 3, 4, 5, 4, 5, 5, 5, 6, 7, 8] 
 [2, 3, 4, 4, 5, 6, 6, 5, 6, 7, 7, 8, 9, 9, 8, 9]]
=#

# ==========
# Matrix of member length L (m) and angle θ (rad)
# L = √((xⱼ - xᵢ)² + (yⱼ - yᵢ)²)
L = zeros(m);
# θ = tan⁻¹ ( (yⱼ - yᵢ)/(xⱼ - xᵢ) )
θ = zeros(m);
for i ∈ 1:m
    ci = C[i, 1];
    cj = C[i, 2];
    L[i] = sqrt( (x[cj]-x[ci])^2 + (y[cj] - y[ci])^2 );
    θ[i] = atan((y[cj] - y[ci]), (x[cj] - x[ci]));
end

function k_local(e)
    a = θ[e];
    k1 = [cos(a)^2, cos(a)*sin(a), -cos(a)^2, -cos(a)*sin(a)];
    k2 = [cos(a)*sin(a), sin(a)^2, -cos(a)*sin(a), -sin(a)^2];
    k3 = [-cos(a)^2, -cos(a)*sin(a), cos(a)^2, cos(a)*sin(a)];
    k4 = [-cos(a)*sin(a), -sin(a)^2, cos(a)*sin(a), sin(a)^2];
    ktemp = [k1 k2 k3 k4];
    k = (1/L[e]) * ktemp;
    return k
end

global K = zeros(dof, dof);

for ℵ ∈ 1:m
    ke = k_local(ℵ);
    Ke = zeros(dof,dof);
    i = C[ℵ, 1];
    j = C[ℵ, 2];
    # k11:k22
    Ke[2(i-1)+1,2(i-1)+1] = ke[1,1];
    Ke[2(i-1)+1,2(i-1)+2] = ke[1,2];
    Ke[2(i-1)+2,2(i-1)+1] = ke[2,1];
    Ke[2(i-1)+2,2(i-1)+2] = ke[2,2];
    # k13:k24
    Ke[2(i-1)+1,2(j-1)+1] = ke[1,3];
    Ke[2(i-1)+1,2(j-1)+2] = ke[1,4];
    Ke[2(i-1)+2,2(j-1)+1] = ke[2,3];
    Ke[2(i-1)+2,2(j-1)+2] = ke[2,4];
    # k31:k42
    Ke[2(j-1)+1,2(i-1)+1] = ke[3,1];
    Ke[2(j-1)+1,2(i-1)+2] = ke[3,2];
    Ke[2(j-1)+2,2(i-1)+1] = ke[4,1];
    Ke[2(j-1)+2,2(i-1)+2] = ke[4,2];
    # k33:k44
    Ke[2(j-1)+1,2(j-1)+1] = ke[3,3];
    Ke[2(j-1)+1,2(j-1)+2] = ke[3,4];
    Ke[2(j-1)+2,2(j-1)+1] = ke[4,3];
    Ke[2(j-1)+2,2(j-1)+2] = ke[4,4];

    # Add to global K
    global K += Ke;
end

# ==========
# Boundary Conditions
kₛ = 0; # Spring constant 24e6
P1 = 100e3; # Right-facing force on nodes 4 and 7
P2 = 200e3; # Downward force on nodes 7 and 9
P3 = 300e3; # Downward force on node 2
EA = 120e6; # Axial stiffness for every member

# Inclusion of EA in K 
K = EA .* K;

q = zeros(dof); # nodal displacement vector
Q = zeros(dof); # nodal force vector

# Boundary conditions for Q
Q[4] = -P3;
Q[7] = P1;
Q[13] = P1;
Q[14] = -P2;
Q[18] = -P2;

# Modify K to fit boundary conditions for q 
# q1
K[1,:] .= 0;
K[:,1] .= 0;
K[1,1] = 1;
# q2
K[2,:] .= 0;
K[:,2] .= 0;
K[2,2] = 1;
# q5 (only for when pinned)
# K[5,:] .= 0;
# K[:,5] .= 0;
# K[5,5] = 1;
K[6,:] .= 0;
K[:,6] .= 0;
K[6,6] = 1;

# ==========
# Solving for displacement q
q = K \ Q;

# Boundary conditions for q
q[1] = 0; # Node is pinned, no displacement
q[2] = 0; # Node is pinned, no displacement
# q[5] = 0; # 0 when pinned, nonzero when roller
q[6] = 0; # 0 for both pin and roller, does not move in vertical direction

# Inclusion of spring constant on node 9
Q[17] = -kₛ*q[17];

# ==========
# Compute deformation coordinates using q
x₊ = zeros(n);
y₊ = zeros(n);
for i ∈ 1:dof
    if isodd(i)
        xi = Int64((i+1)/2);
        x₊[xi] = x[xi] + q[i];
    elseif iseven(i)
        yi = Int64(i/2);
        y₊[yi] = y[yi] + q[i];
    end
end

# Compute the new member lengths
L₊ = zeros(m);
for i ∈ 1:m
    ci = C[i, 1];
    cj = C[i, 2];
    L₊[i] = sqrt( (x₊[cj]-x₊[ci])^2 + (y₊[cj] - y₊[ci])^2 );
end

ϵ = (L₊ .- L) ./ L; # Axial strain for every member
F = EA .* ϵ; # Axial force for every member

Fₖ = F .* 1e-3; # Convert to kN

println(" ")

deform = DataFrame(A=x, B=y, C=x₊, D=y₊);

# CSV.write("deform3.csv", deform)

member = DataFrame(A=L, B=L₊, C=ϵ*1e3, D=Fₖ);

# CSV.write("member3.csv", member)

println("Done!")