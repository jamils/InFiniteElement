# Plotting packages
using Makie
using AbstractPlotting
using GLMakie

# Guassian Quadrature package to compute integrals
using QuadGK

# Used for math
using LinearAlgebra

#=
This program is written to solve the second-order differential equation:
    u" + u - x² = 0
for x ∈ [0,1]
=#

print("Defining things")
# ==========
# Input parameters
a = 0;  # Lower bound
b = 1;  # Upper bound
λ = 3;  # Number of nodes
m = 10;  # Number of elements

hₑ = 1/m;   # Step size

n = m*(λ-1) + 1;   # The size of the global matrix and vectors is m(λ-1)+1

# Boundary conditions
    ## Essential (Dirichlet) boundary conditions
    xa = 1/2;   # Lower essential boundary condition at u(a)
    xb = 1;     # Upper essential boundary condition at u(b)

    essential_a = 0;    # 1 for lower essential true, 0 for false
    essential_ab = 1;   # Both essential boundary conditions true

    ## Natural (neumann) boundary conditions
    qa = 0;     # Lower natural boundary condition at u'(a)
    qb = 0;   # Upper natural boundary condition at u'(b)
    qᵥ = zeros(n);    # Vector for the natural boundary condition
    qᵥ[begin] = -qa;    # First item in vector
    qᵥ[end] = qb;       # Last item in the vector
# ==========
print(".")

# ==========
global u = zeros(n);       # The solution of the differential equation
global K = zeros(n, n);   # The global K vector
global β = zeros(n);       # The global b vector
                            # Using β to not conflict with upper limit b

global xᵥ = a:hₑ:b;     # Array of xᵢ node points
print(".")

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
end
println(".")

# print("Finding K and β")
# ==========
# Solving along each λ-node elements
for ℵ ∈ 0:(m-1)     # Indicie of node, essentiall the i for xᵢ
    for ℶ ∈ 1:λ     # Looping through each node
        for ℷ ∈ 1:λ     # Looping through each node, again
            k(s) = (hₑ/2) * ( (N(ℶ, s) * N(ℷ, s)) - (4/(hₑ^2) * dN(ℶ, s) * dN(ℷ, s)) );
            ktemp, err = quadgk(s -> k(s), -1, 1);
            # Solves the integral ∫(hₑ/2)(NᵢNⱼ - (4/hₑ²) dNᵢ/ds dNⱼ/ds) ds from -1 to 1

            if λ == 2
                K[ℵ+ℶ, ℵ+ℷ] += ktemp;
            elseif ℵ == 0
                K[ℶ, ℷ] += ktemp; # Sum the above to get Kᵢⱼ
            else
                K[(ℵ*(λ-1))+ℶ, (ℵ*(λ-1))+ℷ] += ktemp;
            end
        end

        β_int(s) = (1/4) * (hₑ*s + (xᵥ[begin+ℵ] + xᵥ[begin+ℵ+1]))^2 * N(ℶ, s) * (hₑ/2);
        βtemp, err = quadgk(s -> β_int(s), -1, 1);
        # Computes the integral ∫(1/4)(hₑs + (xᵢ+xᵢ₊₁))² Nᵢ (hₑ/2)ds from -1 to 1
        # println("xᵢ = ", xᵥ[ℵ], ", ℶ = ", ℶ)
        # println("β = ", βtemp)

        if λ == 2
            β[ℵ+ℶ] += βtemp;
        elseif ℵ == 0
            β[ℶ] += βtemp;
        else
            β[(ℵ*(λ-1))+ℶ] += βtemp;
        end
    end
    # print(".")
end
println(" ")

# println("Global K = ")
# display(K)

# println("Global β = ")
# display(β)

# ====================
K = -K;
β = -β;
# ====================

# Setting natural boundary conditions
if essential_ab == 0
    f = qᵥ .+ β;
end

println("Solving")
# ==========
# Setting essential boundary conditions
if essential_a == 1
    f = (-xa .* K[:,begin]) .+ f;
    K[begin, :] .= 0;
    K[:, begin] .= 0;
    K[begin, begin] = 1;
    f[begin] = xa;
elseif essential_ab == 1
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
end

# println("β with BC = ")
# display(β)

# println("K with BC = ")
# display(K)

# println("f = ")
# display(f)

# Solving the matrix equation Ku=f for u
u = K \ f;

# println("u = ")
# display(u)

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
xus = range(a, b, length=n);
lxs = length(xs);
u_exact_array = zeros(lxs);

for i ∈ 1:lxs
    u_exact_array[i] = uₑ₁(xs[i]);
end

println("Plotting")

scene = lines(xs, u_exact_array, color = :black)
lines!(scene, xus, u, color = 1:20) 