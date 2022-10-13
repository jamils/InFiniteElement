# Plotting packages
using Makie
using AbstractPlotting
using GLMakie

# Guassian Quadrature package to compute integrals
@everywhere using QuadGK

# Used for math
using LinearAlgebra

# Used to parallelize code
using Distributed
using SharedArrays


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
m = 10;  # Number of elements

# Boundary conditions
    ## Essential (Dirichlet) boundary conditions
    xa = 1/2;   # Lower essential boundary condition at u(a)
    xb = 1;     # Upper essential boundary condition at u(b)

    essential_a = 0;    # 1 for lower essential true, 0 for false.
    essential_ab = 1;   # Both essential conditions true

    ## Natural (Neumann) boundary conditions
    qa = 0;    # Lower natural boundary condition at u'(a)
    qb = 0;   # Upper natural boundary condition at u'(b)
    qᵥ = zeros(m+1);
    qᵥ[begin] = -qa;
    qᵥ[end] = qb;

# ==========

hₑ = 1/m;   # Step size

global u = zeros(m+1);
global K = zeros(m+1, m+1);
global β = zeros(m+1);

global xᵥ = a:hₑ:b;

# N1(x, xᵢ₊₁) = (xᵢ₊₁ + x) / hₑ;
# N2(x, xᵢ)   = (x - xᵢ)   / hₑ;

function N(ℸ, x, xᵢ, xᵢ₊₁)
    if ℸ == 1
        (xᵢ₊₁ + x) / hₑ;
    elseif ℸ == 2
        (x - xᵢ)   / hₑ; ℸ
    end
end

dN = [-(1/2) * (2/hₑ), (1/2) * (2/hₑ)];

# ==========
# Solving along each 2-node element

for ℵ ∈ 1:m 
    println("ℵ = ", ℵ)
    if isodd(ℵ)
        ℸ = 1;
    elseif iseven(ℵ)
        ℸ = 2;
    end

    for ℶ ∈ 1:k
        println("ℶ = ", ℶ)
        for ℷ ∈ 1:k
            println("ℷ = ", ℷ)
            k₁ = (N(ℷ, xᵥ[ℵ+1], xᵥ[ℵ], xᵥ[ℵ+1]) * dN[ℶ]) - (N(ℷ, xᵥ[ℵ], xᵥ[ℵ], xᵥ[ℵ+1]) * dN[ℶ]);
            # Finds Nᵢ dNⱼ/dx from xᵢ to xᵢ₊₁

            k₂(x) = (N(ℷ, x, xᵥ[ℵ], xᵥ[ℵ+1]) * N(ℶ, x, xᵥ[ℵ], xᵥ[ℵ+1])) + (dN[ℸ] * dN[ℶ]);
            println(k₂(1))
            ktemp, err = quadgk(x -> k₂(x), xᵥ[ℵ], xᵥ[ℵ+1]);
            # Solves the integral ∫(NᵢNⱼ - dNᵢ/dx dNⱼ/dx) dx from xᵢ to xᵢ₊₁

            K[ℵ+(ℶ-1),ℵ+(ℷ-1)] += (k₁ + ktemp); # Sum the above to get Kᵢⱼ
        end
    end

    f_integrand(x) = x^2 * N(ℸ, x, xᵥ[ℵ], xᵥ[ℵ+1]);
    ftemp, err = quadgk(x -> f_integrand(x), xᵥ[ℵ], xᵥ[ℵ+1]);
    # Computes the integral ∫ x²Nᵢ dx from xᵢ to xᵢ₊₁

    β[ℵ] += ftemp;
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

u = K \ (qᵥ .+ β);

if essential_a == 1
    u[begin] = xa;
elseif essential_ab == 1
    u[begin] = xa;
    u[end] = xb;
end

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
    u_exact_array[i] = uₑ₁(xs[i]);
end

println("Plotting")

scene = lines(xs, u_exact_array, color = :black)
lines!(scene, a:hₑ:b, u, color = :blue)