# Used for math
using LinearAlgebra

# Used to output CSV files
using DataFrames
using CSV

# ==========
# Input parameters
m = 21; # Number of members in the system
n = 12; # Number of nodes
dof = 2*n; # Degrees of freedom, two directions for each node
b = 2; # Node distance
h = 3; # Height

# Matrix of node coordinates
# For each node, gives x and y coordinates
#    1  2  3    4    5    6    7    8    9   10   11   12
x = [0, b, b, 2*b, 2*b, 3*b, 3*b, 4*b, 4*b, 5*b, 5*b, 6*b];
y = [0, h, 0, 2*h,  0,  2*h,   0, 2*h,   0,   h,   0,   0];
# Changing x of 6 to incorrect value from homework solution, should be (a/3)+a

# Matrix that gives the beam connections, each item is a node value
C = [[1, 1, 2, 3, 2, 2, 4, 5, 4, 4, 6, 7, 7, 6, 8,  9,  9,  8, 10, 11, 10] [3, 2, 3, 5, 5, 4, 5, 7, 7, 6, 7, 9, 8, 8, 9, 11, 10, 10, 11, 12, 12]];
#=
  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15  16  17  18  19  20  21
[[1, 1, 2, 3, 2, 2, 4, 5, 4, 4, 6, 7, 7, 6, 8,  9,  9,  8, 10, 11, 10] 
 [3, 2, 3, 5, 5, 4, 5, 7, 7, 6, 7, 9, 8, 8, 9, 11, 10, 10, 11, 12, 12]]
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
kₛ = 30e6;  # Spring constant
P1 = 200e3; # 
P2 = 100e3; # 
P3 = 300e3; # 
EA = 120e6; # Axial stiffness for every member

# Inclusion of EA in K 
K = EA .* K;

q = zeros(dof); # nodal displacement vector
Q = zeros(dof); # nodal force vector

# Boundary conditions for Q
Q[3] = P1;
Q[4] = -P2;
Q[5] = -P2;
Q[7] = P1;
Q[8] = -P2;
Q[11] = -P2;
Q[13] = -P3;
Q[15] = -P2;
Q[21] = -P2;


# Modify K to fit boundary conditions for q 
# q1
K[1,:] .= 0;
K[:,1] .= 0;
K[1,1] = 1;
# q2
K[2,:] .= 0;
K[:,2] .= 0;
K[2,2] = 1;

K[end-1,:] .= 0;
K[:,end-1] .= 0;
K[end-1,end-1] = 1;
K[end,:] .= 0;
K[:,end] .= 0;
K[end,end] = 1;

# ==========
# Solving for displacement q
q = K \ Q;

# Boundary conditions for q
q[1] = 0; # Node is pinned, no displacement
q[2] = 0; # Node is pinned, no displacement
q[end-1] = 0;
q[end] = 0;

# Inclusion of spring constant on node 9
Q[9] = -kₛ*q[9];
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

deform = DataFrame(x=x, y=y, xp=x₊, yp=y₊);

# CSV.write("deform3.csv", deform)

member = DataFrame(L=L, Lp=L₊, ϵ=ϵ*1e3, F=Fₖ);

# CSV.write("member3.csv", member)

println("Done!")