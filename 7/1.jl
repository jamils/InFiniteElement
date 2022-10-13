using LinearAlgebra
using QuadGK
using DataFrames

# Solving for the first problem

#=
    This problem gives a non-uniform bar with three distinct sections: 
    AC, CD, and DB. AC and DB can be considered the same, just on
    opposite sides of CD. This will be evaluated with a total of five
    elements, one for AC, one for DB, and three for CD. The local
    parameters for each of these sections will be evaluated separately,
    and joined together for a final computation.

    Here-on, we will use the following designations for each section:

    AC → ζ
    CD → η
    DB → Θ
=#

# Keeping all units in N and m

# ==========
# Parameters given in the problem statement

a = 2; # (m) Length of ζ and Θ in meters
b = 6; # (m) Length of η in meters
w = 6000; # (N/m) Downward force on η, orginially given in kN/m
EI = 1e6; # (N⋅m²) Axial stiffness for ζ and Θ
EIη = 2.5e6; # (N⋅m²) Axial stiffness for η which is equal to 2.5 ⋅ EI

λ = 4; # Number of nodes per element
mη = 3; # m = 5, but 1 element for ζ and Θ each, 3 for η
m = 5; # Total number of elements

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
h = a; # Step size of ζ and Θ
hη = b/mη; # Step size of η
xη = a:hη:(a+b); # x vector of nodes in η
n = m + 1; # Number of nodes in the system
dof = 2*n; # Degrees of freedom

# ==========
# Matrix of node coordinates
x = [0, 2, 4, 6, 8, 10];
y = zeros(6);

# Matrix of beam connections, each item a node value
C = [[1, 2, 3, 4, 5] [2, 3, 4, 5, 6]];

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

# ==========
function ζ(a,EI,λ,h)
    xs(s) = (1/2) * (h*(1+s));

    kζ = zeros(λ,λ);
    βζ = zeros(λ);

    for ℶ ∈ 1:λ
        for ℷ ∈ 1:λ
            k(s) = (2/h) * EI * ((16/hη^4)*(dN(ℶ, s))^2 * (dN(ℷ, s)))^2;
            ktemp, err = quadgk(s -> -k(s), -1, 1);
            kζ[ℶ,ℷ] = ktemp;
        end

        β_int(s) = N(ℶ,s) * (h/2);
        βtemp, err = quadgk(s -> -β_int(s), -1, 1);

        βζ[ℶ] = βtemp;
    end

    return kζ, βζ
end

function η(b,w,mη,EIη,λ,hη,xη)
    dof = λ;
    k1 = zeros(dof,dof);
    k2 = zeros(dof,dof);
    k3 = zeros(dof,dof);

    f1 = zeros(dof);
    f2 = zeros(dof);
    f3 = zeros(dof);

    # Element 1
    for ℶ ∈ 1:λ
        for ℷ ∈ 1:λ
            k(s) = (hη/2) * EIη * ((16/hη^4)*(dN(ℶ, s))^2 * (dN(ℷ, s)))^2;
            ktemp, err = quadgk(s -> -k(s), -1, 1);
            k1[ℶ, ℷ] = ktemp;
        end

        f_int(s) = N(ℶ, s) * (hη/2);
        ftemp, err = quadgk(s -> -f_int(s), -1, 1);

        f1[ℶ] = ftemp;
    end

    # Element 2
    for ℶ ∈ 1:λ
        for ℷ ∈ 1:λ
            k(s) = (hη/2) * EIη * ((16/hη^4)*(dN(ℶ, s))^2 * (dN(ℷ, s)))^2;
            ktemp, err = quadgk(s -> -k(s), -1, 1);
            k2[ℶ, ℷ] = ktemp;
        end

        f_int(s) = N(ℶ, s) * (hη/2);
        ftemp, err = quadgk(s -> -f_int(s), -1, 1);

        f2[ℶ] = ftemp;
    end

    # Element 3
    for ℶ ∈ 1:λ
        for ℷ ∈ 1:λ
            k(s) = (hη/2) * EIη * ((16/hη^4)*(dN(ℶ, s))^2 * (dN(ℷ, s)))^2;
            ktemp, err = quadgk(s -> -k(s), -1, 1);
            k3[ℶ, ℷ] = ktemp;
        end

        f_int(s) = N(ℶ, s) * (hη/2);
        ftemp, err = quadgk(s -> -f_int(s), -1, 1);

        f3[ℶ] = ftemp;
    end

    return k1, f1, k2, f2, k3, f3
end

function Θ(a,EI,λ,h)
    xs(s) = (1/2) * (h*(1+s));

    kΘ = zeros(λ,λ);
    βΘ = zeros(λ);

    for ℶ ∈ 1:λ
        for ℷ ∈ 1:λ
            k(s) = (2/h) * EI * ((16/hη^4)*(dN(ℶ, s))^2 * (dN(ℷ, s)))^2;
            ktemp, err = quadgk(s -> -k(s), -1, 1);
            kΘ[ℶ,ℷ] = ktemp;
        end

        β_int(s) = N(ℶ,s) * (h/2);
        βtemp, err = quadgk(s -> -β_int(s), -1, 1);

        βΘ[ℶ] = βtemp;
    end

    return kΘ, βΘ
end

# ==========
kζ, βζ = ζ(a,EI,λ,h);

kη1, βη1, kη2, βη2, kη3, βη3 = η(b,w,mη,EIη,λ,hη,xη);

kΘ, βΘ = Θ(a,EI,λ,h);

karr = [kζ, kη1, kη2, kη3, kΘ]
βarr = [βζ, βη1, βη2, βη3, βΘ];

K = zeros(dof,dof);
β = zeros(dof);

for i ∈ 1:m;
    ib = i + (i-1);
    ie = ib+3;
    K[ib:ie,ib:ie] = @. K[ib:ie,ib:ie] + karr[i];
    β[ib:ie] = @. β[ib:ie] + βarr[i];
end

Q = zeros(dof);
evec = [4, 6, 8, 10];
for i ∈ evec
    Q[i] = w;
end

fvec = @. β + Q;

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
K[end-1,:] .= 0;
K[:,end-1] .= 0;
K[end-1,end-1] = 1;
fvec[end-1] = 0;

K[end,:] .= 0;
K[:,end] .= 0;
K[end,end] = 1;
fvec[end] = 0;

q = K \ fvec;

q[begin] = 0;
q[begin+1] = 0;
q[end-1] = 0;
q[end] = 0;

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

deform = DataFrame(x=x, y=y, xp=x₊, yp=y₊);

# CSV.write("deform3.csv", deform)

member = DataFrame(L=L, Lp=L₊, ϵ=ϵ*1e3, F=Fₖ);

# CSV.write("member3.csv", member)

println("Done!")