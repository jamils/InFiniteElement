# using Makie
using AbstractPlotting
using DataFrames

# ==========
# Matrix of member length L (m) and angle θ (rad)
# L = √((xⱼ - xᵢ)² + (yⱼ - yᵢ)²)
L = zeros(m);
# θ = tan⁻¹ ( (yⱼ - yᵢ)/(xⱼ - xᵢ) )
θ = zeros(m);

for i ∈ 1:m
    ci = C[i, 1];
    cj = C[i, 2];
    L[i] = sqrt( (xvec[cj]-xvec[ci])^2 + (yvec[cj] - yvec[ci])^2 );
    θ[i] = atan((yvec[cj] - yvec[ci]), (xvec[cj] - xvec[ci]));
end

# ==========
# Compute deformation coordinates using q
x₊ = zeros(λ);
y₊ = zeros(λ);

memvec = [1, 4, 7, 10]; # Starting places of each member

for i ∈ 1:m
    x₊[i] = xvec[i] + q[memvec[i]];
    yi = Int64(memvec[i+1]);
    y₊[i] = yvec[i] + q[yi];
end

x₊[end] = xvec[end] + q[memvec[end]];


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
Fₖ = @. EIₑ * ϵ * 1e-3; # Axial force for every member

# Setup the dataframes for the output values
deform = DataFrame(x=xvec, y=yvec, xp=x₊, yp=y₊)
member = DataFrame(L=L, Lp=L₊, ϵ=ϵ*1e3, F=Fₖ)
values = DataFrame(q=q, Q=Qnew)
#=
println("- Displacement and force for each degree of freedom at each node")
display(values)
println("- Node displacement")
display(deform)
println("- New lengths, axial strain, and axial force for each member")
display(member)
=#