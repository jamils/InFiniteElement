using LinearAlgebra
using QuadGK
using DataFrames

# Solving for the second problem

#=
    This problem gives a uniform bar of length a+b
    with a non-uniform downward force f(x). There
    are three marked locations A, B, and C. A is 
    pinned to a wall at x=0. B is on a roller, at
    x=a. C is the very end of the bar with no
    attachment. 

    For this, we will use m=5 elements, splitting
    such at there are three elements in AB and two
    in BC. With the given lengths a=6m and b=4m,
    each element is 2m in length.
=#

# Keeping all units in N and m

# ==========
# Parameters given in the problem statement

a = 6; # (m) Length of AB
b = 4; # (m) Length of BC
L = a+b; # (m) Total lengh of the bar
EI = 5e6; # (N⋅m²) Axial stiffness is uniform

f(x) = 20*(1+3*((a+b-x)/(a+b))^2) * 1e3; # (N/m)

λ = 4; # Number of nodes per element
m = 5; # Number of elements
hₑ = L/m; # Step size
n = m+1; # Degrees of freedom, two for every node in the system
dof = 2*n;

# ==========
# Matrix of node coordinates
x = zeros(n)
for i ∈ 1:(n-1) 
    x[i+1] = x[i]+hₑ;
end
y = zeros(n);

# Matrix of beam connections, each item a node value
C = [[1, 2, 3, 4, 5] [2, 3, 4, 5, 6]];

# ==========
# ==========
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

function compute(xL, xR, EI, λ, h)

    xs(s) = (1/2) * ((2 * xR) + h * (1+s));

    k = zeros(λ,λ);
    β = zeros(λ);

    for ℶ ∈ 1:λ
        for ℷ ∈ 1:λ
            kf(s) = (2/h) * EI * ((16/h^4)*(dN(ℶ, s))^2 * (dN(ℷ, s)))^2;
            ktemp, err = quadgk(s -> kf(s), -1, 1);

            k[ℶ,ℷ] = ktemp;
        end

        # β_int(s) = f(xs(s)) * N(ℶ,s) * (h/2);
        β_int(s) = N(ℶ,s) * (h/2);
        βtemp, err = quadgk(s -> β_int(s), -1, 1);

        β[ℶ] = βtemp;
    end

    return k, β
end

# ==========
K = zeros(dof,dof);
β = zeros(dof);

for i ∈ 1:m;
    ib = i + (i-1);
    ie = ib+3;
    kmat, βvec = compute(x[i], x[i+1], EI, λ, hₑ);
    K[ib:ie,ib:ie] = @. K[ib:ie,ib:ie] + kmat;
    β[ib:ie] = @. β[ib:ie] + βvec;
end

#=
for i ∈ 1:dof
    if isodd(i) == true
        β[i] = 0;
    end
end
=#

# β .= -β;

Q = zeros(dof);
fvec = zeros(dof);

for i ∈ 1:dof
    if iseven(i) == true
        temp = Int64((i)/2);
        Q[i] = -f(x[temp]);
    end
end

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

# Right node is fixed in place
K[8,:] .= 0;
K[:,8] .= 0;
K[8,8] = 1;
fvec[8] = 0;

q = K \ fvec;

q[begin] = 0;
q[begin+1] = 0;
q[8] = 0;

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
θ₊ = zeros(m);
for i ∈ 1:m
    ci = C[i, 1];
    cj = C[i, 2];
    L₊[i] = sqrt( (x₊[cj]-x₊[ci])^2 + (y₊[cj] - y₊[ci])^2 );
    θ₊[i] = atan( (y₊[cj] - y₊[ci]), (x₊[cj] - x₊[ci]) );
end

ϵ = (L₊ .- L) ./ L; # Axial strain for every member
F = EI .* ϵ; # Axial force for every member

Fₖ = F .* 1e-3; # Convert to kN

deform = DataFrame(x=x, y=y, xp=x₊, yp=y₊)

# CSV.write("deform3.csv", deform)

member = DataFrame(L=L, Lp=L₊, ϵ=ϵ*1e3, F=Fₖ)

# CSV.write("member3.csv", member)

println("Done!")