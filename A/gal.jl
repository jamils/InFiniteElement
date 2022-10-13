# Used for math
using LinearAlgebra

# Used for symbolic math
using ModelingToolkit

# Guassian Quadrature package to compute integrals
using QuadGK

@variables c1 c2 c3 c4 c5 c6 c7 c8
@variables r i
@derivatives D'~r

a = 1;
b = 6;
m = 8;

φ₀(r) = 2r - 1;
φᵢ(r,i) = (r-a) * (b-r)^i;
dφ(r,i) = -(b-r)^(i-1) * (i*(r-a)+r+b);

# C(c1, c2, c3, c4, c5, c6, c7, c8) = [c1, c2, c3, c4, c5, c6, c7,c8];

u(r) = φ₀(r) + (c1*φᵢ(r,1) + c2*φᵢ(r,2) + c3*φᵢ(r,3) + c4*φᵢ(r,4) + c5*φᵢ(r,5) + c6*φᵢ(r,6) + c7*φᵢ(r,7) + c8*φᵢ(r,8));
du(r) = expand_derivatives(D(u(r)));

B(r,i) = (r^2 * dφ(r,i) * du(r)) - (r * φᵢ(r,i)*du(r)) + (φᵢ(r,i)*u(r)*(r^4-4));
bf(r,i) = 16*r*cos(π*r)*φᵢ(r,i);
qtemp(r,i) = r^2 * φᵢ(r,i) * du(r);
q(i) = qtemp(b,i) - qtemp(a,i);

fvec = zeros(m); 
#=
for i ∈ 1:m
    btemp, err = quadgk(r -> bf(r,i), a, b);
    println(btemp)
    println(q(i))
    # fvec[i] = q(i) - btemp;
end
=#
K = zeros(m,m);

for i ∈ 1:m
    # B1
    global c1 = 1;
    global c2 = 0;
    global c3 = 0;
    global c4 = 0;
    global c5 = 0;
    global c6 = 0;
    global c7 = 0;
    global c8 = 0;
    temp, err = quadgk(r -> B(r,i), a, b);
    K[1,i] = temp;

    # B2
    global c1 = 0;
    global c2 = 1;
    global c3 = 0;
    global c4 = 0;
    global c5 = 0;
    global c6 = 0;
    global c7 = 0;
    global c8 = 0;
    temp, err = quadgk(r -> B(r,i), a, b);
    K[2,i] = temp;

    # B3
    global c1 = 0;
    global c2 = 0;
    global c3 = 1;
    global c4 = 0;
    global c5 = 0;
    global c6 = 0;
    global c7 = 0;
    global c8 = 0;
    temp, err = quadgk(r -> B(r,i), a, b);
    K[3,i] = temp;

    # B4
    global c1 = 0;
    global c2 = 0;
    global c3 = 0;
    global c4 = 1;
    global c5 = 0;
    global c6 = 0;
    global c7 = 0;
    global c8 = 0;
    temp, err = quadgk(r -> B(r,i), a, b);
    K[4,i] = temp;

    # B5
    global c1 = 0;
    global c2 = 0;
    global c3 = 0;
    global c4 = 0;
    global c5 = 1;
    global c6 = 0;
    global c7 = 0;
    global c8 = 0;
    temp, err = quadgk(r -> B(r,i), a, b);
    K[5,i] = temp;

    # B6
    global c1 = 0;
    global c2 = 0;
    global c3 = 0;
    global c4 = 0;
    global c5 = 0;
    global c6 = 1;
    global c7 = 0;
    global c8 = 0;
    temp, err = quadgk(r -> B(r,i), a, b);
    K[6,i] = temp;

    # B7
    global c1 = 0;
    global c2 = 0;
    global c3 = 0;
    global c4 = 0;
    global c5 = 0;
    global c6 = 0;
    global c7 = 1;
    global c8 = 0;
    temp, err = quadgk(r -> B(r,i), a, b);
    K[7,i] = temp;

    # B8
    global c1 = 0;
    global c2 = 0;
    global c3 = 0;
    global c4 = 0;
    global c5 = 0;
    global c6 = 0;
    global c7 = 0;
    global c8 = 1;
    temp, err = quadgk(r -> B(r,i), a, b);
    K[8,i] = temp;
end

cvec = [c1, c2, c3, c4, c5, c6, c7, c8];

for i ∈ 1:m
    for j ∈ 1:m
        cvec[:] = 0;
        cvec[i] = 1;
        temp, err = quadgk(r -> B(r,i), a, b);
        K[j,i] = temp;
    end
end