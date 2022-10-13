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


# ==========
# Input parameters
k = 7; # Size of C matrix
m = 2; # Number of subdivisions
a = 0; # Lower bound
b = 1; # Upper bound
# ==========

# Boundary Conditions
ua = 0;
ub = 2;
u_exact(x) = (1/sin(b-a)) * ((ub + b^2 - 2)*sin(x-a) + (ua + a^2 - 2)*sin(b-x) - (x^2-2)*sin(b-a));

sub_vec = SharedArray{Float64}(m);
subdiv = (b-a) / m;

@sync @distributed for i ∈ 1:m
    sub_vec[i] = i * subdiv;
end

xs = a:0.001:b;
lxs = length(xs);

sxs = 1:Int64(floor((lxs-1)/m));

# utot = zeros(m,Int64((lxs-1)/m));
global u_final = zeros(1);

for i ∈ 1:m
    bb = sub_vec[i];
    aa = bb - subdiv;

    println("i = ", i)
    println("aa = ", aa)
    println("bb = ", bb)

    Cmat = SharedArray{Float64}(k,k);
    RH = SharedArray{Float64}(k);

    @everywhere B(x,i,j) = sin(i*π*x)*sin(j*π*x);
    @everywhere ℓ(x,i) = -x^2 * sin(i*π*x) + 2*x*sin(i*π*x);

    println("Computing Matrix")

    @sync @distributed for i ∈ 1:k
        for j ∈ 1:k
            Cmat[i,j], err = quadgk(x -> B(x,i,j), aa, bb);
        end

        RH[i], err = quadgk(x -> ℓ(x,i), aa,bb);
    end

    C = Cmat \ RH;

    φi(x) = sin.(((1:k) .* (π *x)));

    u_approx(x) = (2*x) + sum(φi(x) .* C);

    xxs = aa:0.001:(bb-0.001);

    # u_temp = zeros(length(xxs));

    # push!(u_final, u_approx(xxs[1]));
    # deleteat!(u_final, 1);

    for j ∈ 1:length(xxs)
        # u_temp = u_approx(xxs[j]);
        push!(u_final, u_approx(xxs[j]));
    end

    # uu = append!(u_final, u_temp);
    # u_final = uu;

    # if i == 1
    #     global u_init = u_temp;
    # else
    #     u_final = append!(u_init, u_temp);
    # end

    # u_init = u_final;
end

# u_temp = zeros(k);

# u_temp = xtot[1,:];

# for i ∈ 2:m
#     global u_final = append!(u_temp, xtot[i,:]);
# end

# Plotting section

# for i ∈ 1:(m)
    # push!(u_final, 2)
# end

println("Computing plot arrays")

u_exact_array = zeros(lxs);
# u_approx_array = zeros(lxs);
u_approx_array = u_final;

for i ∈ 1:lxs
    u_exact_array[i] = u_exact(xs[i]);
    # u_approx_array[i] = u_approx(xs[i]);
end

println("Plotting")

scene = lines(xs, u_exact_array, color = :black)
lines!(scene, xs, u_approx_array, color = :blue)