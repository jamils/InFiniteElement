using LinearAlgebra
using QuadGK
using DataFrames
using ModelingToolkit
using Plots

# Solving for the second problem

#=
    This problem gives a uniform bar of length 4a
    with 4 distinct sections of length a. These are
    separated at points A, B, C, D, E. A is a
    pinned node. At all other nodes there is a
    spring of constant k.

    For this, we will use m=8 elements, splitting
    such at there are two elements for each length
    a. With the given length a=2m, each element is
    1m in length.
=#

# Keeping all units in N and m

# ==========
# Parameters given in the problem statement

a = 2; # (m) Length of each element
L = 4*a; # (m) Total lengh of the bar
EI = 30e6; # (N⋅m²) Axial stiffness is uniform
k = 10e-3 * EI; # (N/m) Spring constant

f(x) = (2 + 6*((4*a-x)/4*a)^2) * 1e3; # (N/m)

λ = 4; # Number of nodes per element
m = 8; # Number of elements
hₑ = L/m; # Step size
n = m+1; # Degrees of freedom, two for every node in the system
dof = 2*n;

# ==========
# Matrix of node coordinates
xv = zeros(n)
for i ∈ 1:(n-1) 
    xv[i+1] = xv[i]+hₑ;
end
y = zeros(n);

# Matrix of beam connections, each item a node value
C = [[1, 2, 3, 4, 5, 6, 7, 8] [2, 3, 4, 5, 6, 7, 8, 9]];

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
    L[i] = sqrt( (xv[cj]-xv[ci])^2 + (y[cj] - y[ci])^2 );
    θ[i] = atan((y[cj] - y[ci]), (xv[cj] - xv[ci]));
end

function compute(xL, xR, EI, λ, h)

    xs(s) = (1/2) * ((2 * xR) + h * (1+s));

    k = zeros(λ,λ);
    β = zeros(λ);

    for ℶ ∈ 1:λ
        for ℷ ∈ 1:λ
            kf(s) = (2/h) * EI * ((16/h^4)*(dN(ℶ, s))^2 * (dN(ℷ, s)))^2;
            ktemp, err = quadgk(s -> -kf(s), -1, 1);

            k[ℶ,ℷ] = ktemp;
        end

        β_int(s) = f(xs(s)) * N(ℶ,s) * (h/2);
        βtemp, err = quadgk(s -> -β_int(s), -1, 1);

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
    kmat, βvec = compute(xv[i], xv[i+1], EI, λ, hₑ);
    K[ib:ie,ib:ie] = @. K[ib:ie,ib:ie] + kmat;
    β[ib:ie] = @. β[ib:ie] + βvec;
end

# for i ∈ 1:dof
#     if isodd(i) == true
#         β[i] = 0;
#     end
# end

# β .= -β;

Q = zeros(dof);
fvec = zeros(dof);

# for i ∈ 1:dof
#     if iseven(i) == true
#         temp = Int64((i)/2);
#         Q[i] = f(xv[temp]);
#     end
# end

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

q = K \ fvec;


Q[6]  = k * q[6];
Q[10] = k * q[10];
Q[14] = k * q[14];
Q[18] = k * q[18];

K[6,6]   += k;
K[10,10] += k;
K[14,14] += k;
K[18,18] += k;

@. fvec += Q;

q = K \ fvec;

q[begin] = 0;
q[begin+1] = 0;


# ==========
# Compute deformation coordinates using q
x₊ = zeros(n);
y₊ = zeros(n);
for i ∈ 1:dof
    if isodd(i)
        xi = Int64((i+1)/2);
        x₊[xi] = xv[xi] + q[i];
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

deform = DataFrame(x=xv, y=y, xp=x₊, yp=y₊)

member = DataFrame(L=L, Lp=L₊, ϵ=ϵ*1e3, F=Fₖ)

@variables x 
@derivatives D'~x

N1(x,i) = (1/16) * (-1 + (2*x - 2*xv[i] - xv[i+1])/hₑ + 9*((2*x - 2*xv[i] - xv[i+1])^2)/(hₑ^2) - 9*(2*x - 2*xv[i] - xv[i+1])^3/(hₑ^3));

N2(x,i) = (9/16) * (1 - 3*(2*x - 2*xv[i] - xv[i+1])/hₑ - ((2*x - 2*xv[i] - xv[i+1])^2)/(hₑ^2) + 3*(2*x - 2*xv[i] - xv[i+1])^3/(hₑ^3));

N3(x,i) = (9/16) * (1 + 3*(2*x - 2*xv[i] - xv[i+1])/hₑ - ((2*x - 2*xv[i] - xv[i+1])^2)/(hₑ^2) - 9*(2*x - 2*xv[i] - xv[i+1])^3/(hₑ^3));

N4(x,i) = (1/16) * (-1 - (2*x - 2*xv[i] - xv[i+1])/hₑ + 9*((2*x - 2*xv[i] - xv[i+1])^2)/(hₑ^2) + 9*(2*x - 2*xv[i] - xv[i+1])^3/(hₑ^3));

#=
for i ∈ 1:m
    jvec = 1:2:(dof-2);
    j1 = jvec[i];
    j2 = j1+1;
    j3 = j1+2;
    j4 = j1+3;
    println("w$i(x) = N1(x,$i)*q[$j1] + N2(x,$i)*q[$j2] + N3(x,$i)*q[$j3] + N4(x,$i)*q[$j4];")
end
=#

w1(x) = N1(x,1)*q[1] + N2(x,1)*q[2] + N3(x,1)*q[3] + N4(x,1)*q[4];
w2(x) = N1(x,2)*q[3] + N2(x,2)*q[4] + N3(x,2)*q[5] + N4(x,2)*q[6];
w3(x) = N1(x,3)*q[5] + N2(x,3)*q[6] + N3(x,3)*q[7] + N4(x,3)*q[8];
w4(x) = N1(x,4)*q[7] + N2(x,4)*q[8] + N3(x,4)*q[9] + N4(x,4)*q[10];
w5(x) = N1(x,5)*q[9] + N2(x,5)*q[10] + N3(x,5)*q[11] + N4(x,5)*q[12];
w6(x) = N1(x,6)*q[11] + N2(x,6)*q[12] + N3(x,6)*q[13] + N4(x,6)*q[14];
w7(x) = N1(x,7)*q[13] + N2(x,7)*q[14] + N3(x,7)*q[15] + N4(x,7)*q[16];
w8(x) = N1(x,8)*q[15] + N2(x,8)*q[16] + N3(x,8)*q[17] + N4(x,8)*q[18];

for i ∈ 1:m
    print("simplify(w$i(x)) ")
end

wtemp = [simplify(w1(x)) simplify(w2(x)) simplify(w3(x)) simplify(w4(x)) simplify(w5(x)) simplify(w6(x)) simplify(w7(x)) simplify(w8(x))];

for i ∈ 1:m
    println("w$i(x) = ", wtemp[i])
end

for i ∈ 1:m
    println("θ$i(x) = expand_derivatives(D(w$i(x)));")
end

for i ∈ 1:m
    print("expand_derivatives(D(w$i(x))) ")
end

ttemp = [expand_derivatives(D(w1(x))) expand_derivatives(D(w2(x))) expand_derivatives(D(w3(x))) expand_derivatives(D(w4(x))) expand_derivatives(D(w5(x))) expand_derivatives(D(w6(x))) expand_derivatives(D(w7(x))) expand_derivatives(D(w8(x)))];

for i ∈ 1:m
    println("θ$i(x) = ", ttemp[i])
end

θ1(x) = expand_derivatives(D(w1(x)));
θ2(x) = expand_derivatives(D(w2(x)));
θ3(x) = expand_derivatives(D(w3(x)));
θ4(x) = expand_derivatives(D(w4(x)));
θ5(x) = expand_derivatives(D(w5(x)));
θ6(x) = expand_derivatives(D(w6(x)));
θ7(x) = expand_derivatives(D(w7(x)));
θ8(x) = expand_derivatives(D(w8(x)));

for i ∈ 1:m
    println("M$i(x) = EI * expand_derivatives(D(θ$i(x)));")
end

M1(x) = EI * expand_derivatives(D(θ1(x)));
M2(x) = EI * expand_derivatives(D(θ2(x)));
M3(x) = EI * expand_derivatives(D(θ3(x)));
M4(x) = EI * expand_derivatives(D(θ4(x)));
M5(x) = EI * expand_derivatives(D(θ5(x)));
M6(x) = EI * expand_derivatives(D(θ6(x)));
M7(x) = EI * expand_derivatives(D(θ7(x)));
M8(x) = EI * expand_derivatives(D(θ8(x)));

for i ∈ 1:m
    print("M$i(x) ")
end

mtemp = [M1(x) M2(x) M3(x) M4(x) M5(x) M6(x) M7(x) M8(x)];

for i ∈ 1:m
    println("M$i(x) = ", mtemp[i])
end

for i ∈ 1:m
    println("V$i(x) = -EI * expand_derivatives(D(M$i(x)));")
end

V1(x) = -EI * expand_derivatives(D(M1(x)));
V2(x) = -EI * expand_derivatives(D(M2(x)));
V3(x) = -EI * expand_derivatives(D(M3(x)));
V4(x) = -EI * expand_derivatives(D(M4(x)));
V5(x) = -EI * expand_derivatives(D(M5(x)));
V6(x) = -EI * expand_derivatives(D(M6(x)));
V7(x) = -EI * expand_derivatives(D(M7(x)));
V8(x) = -EI * expand_derivatives(D(M8(x)));

for i ∈ 1:m
    print("V$i(x) ")
end

vtemp = [V1(x) V2(x) V3(x) V4(x) V5(x) V6(x) V7(x) V8(x)];

for i ∈ 1:m
    println("V$i(x) = ", vtemp[i])
end

# Plotting section
for i ∈ 0:(m-1)
    j = i+1;
    print("$i:0.01:$j ")
end

for i ∈ 1:m
    print("w$i(xplot[:, $i]) ")
end

for i ∈ 1:m
    print("θ$i(xplot[:, $i]) ")
end

for i ∈ 1:m
    print("M$i(xplot[:, $i]) ")
end

for i ∈ 1:m
    print("V$i(xplot[:, $i]) ")
end

xplot = [0:0.01:1 1:0.01:2 2:0.01:3 3:0.01:4 4:0.01:5 5:0.01:6 6:0.01:7 7:0.01:8];

wplot = @. [w1(xplot[:, 1]) w2(xplot[:, 2]) w3(xplot[:, 3]) w4(xplot[:, 4]) w5(xplot[:, 5]) w6(xplot[:, 6]) w7(xplot[:, 7]) w8(xplot[:, 8])];

θplot = @. [θ1(xplot[:, 1]) θ2(xplot[:, 2]) θ3(xplot[:, 3]) θ4(xplot[:, 4]) θ5(xplot[:, 5]) θ6(xplot[:, 6]) θ7(xplot[:, 7]) θ8(xplot[:, 8])];

Mplot = @. [M1(xplot[:, 1]) M2(xplot[:, 2]) M3(xplot[:, 3]) M4(xplot[:, 4]) M5(xplot[:, 5]) M6(xplot[:, 6]) M7(xplot[:, 7]) M8(xplot[:, 8])];

Vplot = @. [V1(xplot[:, 1]) V2(xplot[:, 2]) V3(xplot[:, 3]) V4(xplot[:, 4]) V5(xplot[:, 5]) V6(xplot[:, 6]) V7(xplot[:, 7]) V8(xplot[:, 8])];

# plot(xplot,θplot)

println("Done!")