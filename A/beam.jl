# Used for math
using LinearAlgebra

# Guassian Quadrature package to compute integrals
using QuadGK

# Solving the tapered bar problem
# Keeping all units in N and m

λ = 3; # Number of nodes
m = 4; # Number of elements

L = 1; # 1 meter
E = 7.5e10; # 75 GPa, or 7.5×10^10 Pa = 7.5×10^10 N/m^2
P = 1e6; # 1000 kN, or 1×10^6 N
b = 0.03; # 30 mm or 0.03 m
h₀ = 0.1; # 100 mm or 0.1 m
h(x) = h₀ * exp((-4/5) * x);
f(x) = 1e6 * (1 - (3/4)*x^2); # N/m

EA(x) = E * b * h(x) * 1000; # Axial deformation function

hₑ = 1/m; # Step size
n = m*(λ-1) + 1;   # The size of the global matrix and vectors is m(λ-1)+1

# ==========
global u = zeros(n);      # The solution of the differential equation
global K = zeros(n, n);   # The global K vector
global β = zeros(n);      # The global b vector
                          # Using β to not conflict with upper limit b
# global f = zeros(n);
global xᵥ = 0:hₑ:L;     # Array of xᵢ node points

# ==========
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

k3 = zeros(λ, λ);
f3 = zeros(λ);

# ==========
# Solving along each λ-node elements
for ℵ ∈ 0:(m-1)     # Indicie of node, essential the i for xᵢ
    xs(s) = (1/2) * ((2 * xᵥ[begin+ℵ]) + hₑ * (1+s));
    for ℶ ∈ 1:λ     # Looping through each node
        for ℷ ∈ 1:λ     # Looping through each node, again
            k(s) = (2/hₑ) * EA(xs(s)) * (dN(ℶ, s) * dN(ℷ, s)) ;
            ktemp, err = quadgk(s -> -k(s), -1, 1);
            # Solves the integral ∫(2/hₑ) EA(s) (dNᵢ/ds dNⱼ/ds) ds from -1 to 1

            if ℵ == 0
                K[ℶ, ℷ] += ktemp; # Sum the above to get Kᵢⱼ
            else
                K[(ℵ*(λ-1))+ℶ, (ℵ*(λ-1))+ℷ] += ktemp;
            end

            if ℵ == 3
                k3[ℶ, ℷ] = ktemp;
            end
        end

        β_int(s) = f(xs(s)) * N(ℶ, s) * (hₑ/2);
        βtemp, err = quadgk(s -> -β_int(s), -1, 1);
        # Computes the integral ∫f(s) Nᵢ (hₑ/2)ds from -1 to 1

        # q = N(ℶ,1) * EA(xs(1)) - N(ℶ,-1) * EA(xs(-1));

        if ℵ == 3
            f3[ℶ] = βtemp;
        end

        if ℵ == 0
            β[ℶ] += βtemp;
        else
            β[(ℵ*(λ-1))+ℶ] += βtemp;
        end
    end
end

p = zeros(n);
p[end] = P;

fv = p .+ β;

K[begin, :] .= 0;
K[:, begin] .= 0;
K[begin, begin] = 1;
fv[begin] = 0;

u = K \ fv;