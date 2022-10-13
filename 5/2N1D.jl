# Plotting packages
using Makie
using AbstractPlotting
using GLMakie

# Used to parallelize code
using Distributed
using SharedArrays

# Guassian Quadrature package to compute integrals
@everywhere using QuadGK

# Used for math
using LinearAlgebra


#=
This program is written to solve the second-order differential equation:
    u" + u - x² = 0
for x ∈ [0,1]
=#


# ==========
# Input parameters
a = 0;  # Lower bound
b = 1;  # Upper bound
k = 2;  # Number of nodes
m = 5;  # Number of elements

# Boundary conditions
    ## Essential (Dirichlet) boundary conditions
    xa = 1/2;   # Lower essential boundary condition at u(a)
    xb = 1;     # Upper essential boundary condition at u(b)

    essential_a = 1;    # 1 for lower essential true, 0 for false.
    essential_ab = 0;   # Both essential conditions true

    ## Natural (Neumann) boundary conditions
    qa = 0;    # Lower natural boundary condition at u'(a)
    qb = 4/3;   # Upper natural boundary condition at u'(b)
    qᵥ = zeros(m+1);
    qᵥ[begin] = -qa;
    qᵥ[end] = qb;

# ==========

hₑ = 1/m;   # Step size

global u = zeros(m+1);
global K = zeros(m+1, m+1);
global β = zeros(m+1);

global xᵥ = a:hₑ:b; # array of xᵢ node points

# Sets N and dN functions for number of nodes
if k == 2
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

elseif k == 3
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
end
# ==========
# Solving along each 2-node element

for ℵ ∈ 1:m     # Indicie of node, essentially the i for xᵢ
    for ℶ ∈ 1:k
        for ℷ ∈ 1:k
            # k₁ = ( (2/hₑ) * N(ℶ, 1) * dN(ℷ, 1) ) - ( (2/hₑ) * N(ℶ, -1) * dN(ℷ, -1) );
            # Finds (2/hₑ) Nᵢ dNⱼ/ds from -1 to 1

            k₂(s) = ( (N(ℶ, s) * N(ℷ, s)) + (4/(hₑ^2)) * dN(ℶ, s) * dN(ℷ, s) ) * (hₑ/2);
            ktemp, err = quadgk(s -> k₂(s), -1, 1);
            # Solves the integral ∫(he/2)(NᵢNⱼ + (4/hₑ²) dNᵢ/ds dNⱼ/ds) ds from -1 to 1

            K[ℵ+(ℶ-1),ℵ+(ℷ-1)] += ktemp; #(k₁ - ktemp); # Sum the above to get Kᵢⱼ
            # ADD SYMETRY CHECK
        end
    end

    β_integrand_1(s) = (1/4) * (hₑ*s + (xᵥ[ℵ] + xᵥ[ℵ+1]))^2 * N(1, s) * (hₑ/2);
    βtemp_1, err = quadgk(s -> β_integrand_1(s), -1, 1);
    # Computes the integral ∫ (1/4)(hₑs + (xᵢ+xᵢ₊₁))² N₁(hₑ/2)ds from -1 to 1

    β[ℵ] += - βtemp_1;

    β_integrand_2(s) = (1/4) * (hₑ*s + (xᵥ[ℵ] + xᵥ[ℵ+1]))^2 * N(2, s) * (hₑ/2);
    βtemp_2, err = quadgk(s -> β_integrand_2(s), -1, 1);
    # Computes the integral ∫ (1/4)(hₑs + (xᵢ+xᵢ₊₁))² N₂(hₑ/2)ds from -1 to 1

    β[ℵ+1] = - βtemp_2;
end

# ==========
# Time to set essential boundary conditions
if essential_a == 1
    K[begin, :] .= 0;
    K[:, begin] .= 0;
    K[begin, begin] = xa;
    β[begin] = xa;
elseif essential_ab == 1
    K[begin, :] .= 0;
    K[:, begin] .= 0;
    K[begin, begin] = xa;
    β[begin] = xa;
    # ----------
    K[end, :] .= 0;
    K[:, end] .= 0;
    K[end, end] = xb;
    β[end] = xb;
end 

f = qᵥ .+ β;

u = K \ f;

# if essential_a == 1
#     u[begin] = xa;
# elseif essential_ab == 1
#     u[begin] = xa;
#     u[end] = xb;
# end

# Exact solution for 1.
uₑ₁(x) = -2 + x^2 + (5/2) * cos(x) - (1/2) * (5*cot(1) - 4*csc(1))*sin(x);

# Exact solution for 2.
uₑ₂(x) = -2 + x^2 + (1/3) * (2*csc(1) - 3*cot(1))*cos(x) - sin(x);

# Exact solution for 3.
uₑ₃(x) = -2 + x^2 + (5/2) *cos(x) - (1/6)*(4-15*sin(1))*sec(1)*sin(x);

# ==========
# Plotting section

println("Computing plot arrays")

xs = a:0.001:b;
lxs = length(xs);
u_exact_array = zeros(lxs);

for i ∈ 1:lxs
    u_exact_array[i] = uₑ₃(xs[i]);
end

println("Plotting")

scene = lines(xs, u_exact_array, color = :black)
lines!(scene, a:hₑ:b, u, color = 1:20)