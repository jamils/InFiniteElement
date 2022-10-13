# ==========
# This is a program used to check my answer for question 2 of the first part of the exam
# ==========

# Guassian Quadrature package to compute integrals
using QuadGK

# Used for math
using LinearAlgebra

# ==========
# Input parameters
a = 0;  # Lower bound
b = 9;  # Upper bound
λ = 3;  # Number of nodes
m = 4;  # Number of elements

hₑ = 1/m;   # Step size

n = m*(λ-1) + 1;   # The size of the global matrix and vectors is m(λ-1)+1

# Boundary conditions
    BCType = 1; # 1 for only essential,
                # 2 for only natural,
                # 3 for xa and qb

    ## Essential (Dirichlet) boundary conditions
    xa = 1;   # Lower essential boundary condition at u(a)
    xb = 9;     # Upper essential boundary condition at u(b)

    ## Natural (neumann) boundary conditions
    qa = -1;     # Lower natural boundary condition at u'(a)
    qb = 4/3;   # Upper natural boundary condition at u'(b)
    qᵥ = zeros(n);    # Vector for the natural boundary condition

    if BCType == 2
        qᵥ[begin] = -qa;    # First item in vector
        qᵥ[end] = qb;       # Last item in the vector
    elseif BCType == 3
        qᵥ[end] = qb;       # Last item in the vector
    end
# ==========

# ==========
global u = zeros(n);      # The solution of the differential equation
global K = zeros(n, n);   # The global K vector
global β = zeros(n);      # The global b vector
                          # Using β to not conflict with upper limit b

global xᵥ = a:hₑ:b;     # Array of xᵢ node points

# ==========
# Defining N(s) and dN/ds functions for λ nodes
if λ == 2
    function N(ℸ, s)
        if ℸ == 1
            (1 - s) / 2;
        elseif ℸ == 2
            (1 + s)   / 2;
        end
    end

    function dN(ℸ, s)
        if ℸ == 1
            -(1/2);
        elseif ℸ == 2
            (1/2);
        end
    end

elseif λ == 3
    function N(ℸ, s)
        if ℸ == 1
            - (1/2) * (1-s) * s;
        elseif ℸ == 2
            1 - s^2;
        elseif ℸ == 3
            (1/2) * (1+s) * s;
        end
    end

    function dN(ℸ, s)
        if ℸ == 1
            -(1/2) + s;
        elseif ℸ == 2
            -2*s;
        elseif ℸ == 3
            (1/2) + s;
        end
    end
elseif λ == 4
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
end

# ==========
# Solving along each λ-node elements
for ℵ ∈ 0:(m-1)     # Indicie of node, essential the i for xᵢ
    for ℶ ∈ 1:λ     # Looping through each node
        for ℷ ∈ 1:λ     # Looping through each node, again
            k(s) = (hₑ/2) * ( (N(ℶ, s) * N(ℷ, s)) - (4/(hₑ^2) * dN(ℶ, s) * dN(ℷ, s)) );
            ktemp, err = quadgk(s -> -k(s), -1, 1);
            # Solves the integral ∫(hₑ/2)(NᵢNⱼ - (4/hₑ²) dNᵢ/ds dNⱼ/ds) ds from -1 to 1

            if ℵ == 0
                K[ℶ, ℷ] += ktemp; # Sum the above to get Kᵢⱼ
            else
                K[(ℵ*(λ-1))+ℶ, (ℵ*(λ-1))+ℷ] += ktemp;
            end
        end

        β_int(s) = (1/4) * (hₑ*s + (xᵥ[begin+ℵ] + xᵥ[begin+ℵ+1]))^2 * N(ℶ, s) * (hₑ/2);
        βtemp, err = quadgk(s -> -β_int(s), -1, 1);
        # Computes the integral ∫(1/4)(hₑs + (xᵢ+xᵢ₊₁))² Nᵢ (hₑ/2)ds from -1 to 1

        if ℵ == 0
            β[ℶ] += βtemp;
        else
            β[(ℵ*(λ-1))+ℶ] += βtemp;
        end
    end
end

# Setting natural boundary conditions
f = qᵥ .+ β;

# Setting essential boundary conditions
if BCType == 1
    f = (-xa .* K[:,begin]) .+ β;
    K[begin, :] .= 0;
    K[:, begin] .= 0;
    K[begin, begin] = 1;
    f[begin] = xa;
    # ----------
    f = (-xb .* K[:,end]) .+ f;
    K[end, :] .= 0;
    K[:, end] .= 0;
    K[end, end] = 1;
    f[end] = xb;
elseif BCType == 3
    f = (-xa .* K[:,begin]) .+ f;
    K[begin, :] .= 0;
    K[:, begin] .= 0;
    K[begin, begin] = 1;
    f[begin] = xa;
end