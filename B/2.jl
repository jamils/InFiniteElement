# Guassian Quadrature package to compute integrals
using QuadGK

# Used for math
using LinearAlgebra

using DataFrames

# ======
# Input parameters
a = 2; # (m)
λ = 4; # Number of nodes per element
m = 8; # Number of elements
L = m*a; # (m) Full length of the bar with 8 beam elements
hₑ = L/m; # Step size
n = m+1; # Degrees of freedom, two for every node in the system
dof = 2*n;
w₁ = 8e3; # (N/m)
w₂ = 6e3; # (N/m)
P = 100e3; # (N)
M = 50e3; # (N⋅m)
k = 20e3; # (N/m)

# (N⋅m²) Bending stiffness function
EI(x) = 30*(1-((3*x)/(4*L))) * 10e6; 

# ==========
# Matrix of node coordinates
xᵥ = zeros(n)
for i ∈ 1:(n-1) 
    xᵥ[i+1] = xᵥ[i]+hₑ;
end
y = zeros(n);

# Matrix of beam connections, each item a node value
C = [[1, 2, 3, 4, 5, 6, 7, 8] [2, 3, 4, 5, 6, 7, 8, 9]];

# =====
# Lagrange shape function for 4 nodes
function N(ℸ, s)
    if ℸ == 1
        (-1/16) * (1 - s - 9*s^2 + 9*s^3);
    elseif ℸ == 2
        (9/16) * (1 - 3*s - s^2 + 3*s^3);
    elseif ℸ == 3
        (9/16) * (1 + 3*s - s^2 - 3*s^3);
    elseif ℸ == 4
        (-1/16) * (1 + s - 9*s^2 - 9*s^3);
    end
end

function dN(ℸ, s)
    if ℸ == 1
        (1/16) * (1 + 18*s - 27*s^2);
    elseif ℸ == 2
        (-9/16) * (3 + 2*s - 9*s^2);
    elseif ℸ == 3
        (9/16) * (3 - 2*s - 9*s^2);
    elseif ℸ == 4
        (-1/16) * (1 - 18*s - 27*s^2);
    end
end

# ==========
# Matrix of member length L (m) and angle θ (rad)
ℓ = zeros(m);
θ = zeros(m);
# L = √((xⱼ - xᵢ)² + (yⱼ - yᵢ)²)
# θ = tan⁻¹ ( (yⱼ - yᵢ)/(xⱼ - xᵢ) )
for i ∈ 1:m
    ci = C[i, 1];
    cj = C[i, 2];
    ℓ[i] = sqrt( (xᵥ[cj]-xᵥ[ci])^2 + (y[cj] - y[ci])^2 );
    θ[i] = atan((y[cj] - y[ci]), (xᵥ[cj] - xᵥ[ci]));
end

# EIvec = zeros(length(xᵥ));
# @. EIvec[:] = EI(xᵥ[:]);

function compute(xL, xR, λ, h, i)

    xs(s) = (1/2) * ((2 * xR) + h * (1+s));

    k = zeros(λ,λ);
    β = zeros(λ);

    for ℶ ∈ 1:λ
        for ℷ ∈ 1:λ
            kf(s) = (2/h) * EI(xs(s)) * ((16/h^4)*(dN(ℶ,s))^2 * (dN(ℷ,s)))^2;
            ktemp, err = quadgk(s -> kf(s), -1, 1);

            k[ℶ,ℷ] = ktemp;
        end

        β_int(s) = N(ℶ,s) * (h/2);
        βtemp, err = quadgk(s -> β_int(s), -1, 1);

        β[ℶ] = βtemp;
    end

    println("Element $i")
    display(k)
    println(β)
    println(" ")

    return k, β
end

# ==========
K = zeros(dof,dof);
β = zeros(dof);

for i ∈ 1:m;
    ib = i + (i-1);
    ie = ib+3;
    kmat, βvec = compute(xᵥ[i], xᵥ[i+1], λ, hₑ, i);
    K[ib:ie,ib:ie] = @. K[ib:ie,ib:ie] + kmat;
    β[ib:ie] = @. β[ib:ie] + βvec;
end

Q = zeros(dof);
fvec = zeros(dof);

Q[18] = -P;

Q[2] = -w₁;
Q[4] = -w₁;
Q[6] = -w₁;
Q[8] = -w₁;
Q[12] = -w₂;
Q[14] = -w₂;
Q[16] = -w₂;
Q[18] += -w₂;

@. fvec = β + Q;

# Modify K to fit boundary conditions for q
# Left node is fixed in place
K[begin,:] .= 0;
K[:,begin] .= 0;
K[begin,begin] = 1;
fvec[begin] = 0;

K[begin+1,:] .= 0;
K[:,begin+1] .= 0;
K[begin+1,begin+1] = 1;
fvec[begin+1] = 0;

K[10,:] .= 0;
K[:,10] .= 0;
K[10,10] = 1;
fvec[10] = 0;

q = K \ fvec;


# Include spring constant
Q[4]  = k * q[4];
Q[6] = k * q[6];
Q[8] = k * q[8];
Q[12] = k * q[12];
Q[14] = k * q[14];
Q[16] = k * q[16];
Q[18] = k * q[18];

K[4,4]   += k;
K[6,6] += k;
K[8,8] += k;
K[12,12] += k;
K[14,14] += k;
K[16,16] += k;
K[18,18] += k;

@. fvec += Q;

q = K \ fvec;

q[begin] = 0;
q[begin+1] = 0;
q[10] = 0;


# ==========
# Compute deformation coordinates using q
x₊ = zeros(n);
y₊ = zeros(n);
for i ∈ 1:dof
    if isodd(i)
        xi = Int64((i+1)/2);
        x₊[xi] = xᵥ[xi] + q[i];
    elseif iseven(i)
        yi = Int64(i/2);
        y₊[yi] = y[yi] + q[i];
    end
end

# Compute the new member lengths
L₊ = zeros(m);
θ₊ = zeros(m);
for i ∈ 1:m 
    ci = C[i, 1];
    cj = C[i, 2];
    L₊[i] = sqrt( (x₊[cj]-x₊[ci])^2 + (y₊[cj] - y₊[ci])^2 );
    θ₊[i] = atan( (y₊[cj] - y₊[ci]), (x₊[cj] - x₊[ci]) );
end

ϵ = (L₊ .- L) ./ L; # Axial strain for every member
F = @. EI(x₊[2:end]) * ϵ; # Axial force for every member

Fₖ = F .* 1e-3; # Convert to kN

deform = DataFrame(x=xv, y=y, xp=x₊, yp=y₊)

member = DataFrame(L=L, Lp=L₊, ϵ=ϵ*1e3, F=Fₖ)

