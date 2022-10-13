# Plotting packages
using Makie
using AbstractPlotting
using GLMakie

# Guassian Quadrature package to compute integrals
using QuadGK

# Used for math
using LinearAlgebra

# ==========
# Input parameters
m = 11;       # Number of members in the system
n = 7;        # Number of nodes
dof = 2*n;    # Degrees of freedom, two directions for each nodes

# Matrix of node coordinates
# For each node, gives x and y coordinates
x = [0, 5, (3/2), 5, (7/2), 5, 5];
y = [0, 0, 3, 3, 7, 7, 10];

# Matrix that gives the beam connections, each item is a node value
C = [[1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6] [2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7]];

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

# println("K1,1 = ", K[1,1])
# println("K2,2 = ", K[2,2])
# println("K3,3 = ", K[3,3])
# println("K4,4 = ", K[4,4])
# println("K5,5 = ", K[5,5])
# println("K6,6 = ", K[6,6])
# println("K7,7 = ", K[7,7])
# println("K8,8 = ", K[8,8])
# println("K9,9 = ", K[9,9])
# println("K10,10 = ", K[10,10])
# println("K11,11 = ", K[11,11])
# println("K12,12 = ", K[12,12])
# println("K13,13 = ", K[13,13])
# println("K14,14 = ", K[14,14])

EA = 120e6;
k = 24e6;
P1 = 300e3;
P2 = 200e3;

K = EA .* K;

q = zeros(dof);
Q = zeros(dof);

Q[5] = P1;
Q[9] = P1;
Q[13]= P1;
Q[6]  = -P2;
Q[10] = -P2;
Q[14] = -P2;
Q[11] += -k*q[11];

K[1,:] .= 0;
K[:,1] .= 0;
K[1,1] = 1;
K[2,:] .= 0;
K[:,2] .= 0;
K[2,2] = 1;
K[4,:] .= 0;
K[:,4] .= 0;
K[4,4] = 1;

K[11,11] += k;

q = K \ Q;

q[1] = 0;
q[2] = 0;
q[4] = 0;

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

L₊ = zeros(m);
for i ∈ 1:m
    ci = C[i, 1];
    cj = C[i, 2];
    L₊[i] = sqrt( (x₊[cj]-x₊[ci])^2 + (y₊[cj] - y₊[ci])^2 );
end

ϵ = (L₊ .- L) ./ L;
F = EA .* ϵ;

Fₖ = F .* 1e-3;

println(" ")

deform = DataFrame(A=x, B=y, C=x₊, D=y₊);

# CSV.write("deform3.csv", deform)

member = DataFrame(A=L, B=L₊, C=ϵ*1e3, D=Fₖ);

# CSV.write("member3.csv", member)

println("Done!")