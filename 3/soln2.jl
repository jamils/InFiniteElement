# Plotting packages
using Makie
using AbstractPlotting
using GLMakie

# Guassian Quadrature package to compute integrals
@everywhere using QuadGK

# Used to parallelize code
using Distributed
using SharedArrays


# ==========
# Input parameters
k = 10; # Size of C matrix
a = 0; # Lower bound
b = 1; # Upper bound
# ==========

# Boundary Conditions
ua = 0;
ub = 2;

u_exact(x) = (1/sin(b-a)) * ((ub + b^2 - 2)*sin(x-a) + (ua + a^2 - 2)*sin(b-x) - (x^2-2)*sin(b-a));

Cmat = SharedArray{Float64}(k,k);
RH = SharedArray{Float64}(k);

@everywhere B(x,i,j) = sin(i*π*x)*sin(j*π*x);
@everywhere ℓ(x,i) = -x^2 * sin(i*π*x) + 2*x*sin(i*π*x);


println("Computing Matrix")

@time @sync @distributed for i ∈ 1:k
    for j ∈ 1:k
        Cmat[i,j], err = quadgk(x -> B(x,i,j), a, b);
    end

    RH[i], err = quadgk(x -> ℓ(x,i), a,b);
end

C = Cmat \ RH;

φi(x) = sin.(((1:(k)) .* (π *x))); # sin.(((1:(k*m)) .* (π *x)));

u_approx(x) = (2*x) + sum(φi(x) .* C);

# Plotting section

println("Computing plot arrays")

xs = a:0.001:b;
lxs = length(xs);

u_exact_array = zeros(lxs);
u_approx_array = zeros(lxs);

for i ∈ 1:lxs
    u_exact_array[i] = u_exact(xs[i]);
    u_approx_array[i] = u_approx(xs[i]);
end

println("Plotting")

# scene = lines(xs, u_exact_array, color = :black)
# lines!(scene, xs, u_approx_array, color = :blue)