#=
    This takes an input and computes values for a triangle of the following form:


        3
        .
       / \
      / . \
     /     \
    ∙-------∙
  1           2

=#

using LinearAlgebra
using Plots
using ModelingToolkit
using HCubature

# =====
# Input data (to be changed later)

# Node locations

Nodes = [
    [1 2]
    [6 4]
    [3 6]
    [0.0 0.0]
];

# Now to compute node 4, which is in the center of the triangle

# First we will start with computing the halfway points on each element

halfs = zeros(3,2);
# 1
halfs[1,1] = (Nodes[1,1] + Nodes[2,1]) / 2;
halfs[1,2] = (Nodes[1,2] + Nodes[2,2]) / 2;

# 2
halfs[2,1] = (Nodes[2,1] + Nodes[3,1]) / 2;
halfs[2,2] = (Nodes[2,2] + Nodes[3,2]) / 2;

# 3
halfs[3,1] = (Nodes[3,1] + Nodes[1,1]) / 2;
halfs[3,2] = (Nodes[3,2] + Nodes[1,2]) / 2;

# Now we can use this same logic to size down the triangle infinately to find the center
accuracy = eps(); # Machine epsilon
centerarr = zeros(3,2);
temp = halfs;
maxdiff = 1;

@time while maxdiff > accuracy #1e-4
    #1
    centerarr[1,1] = (temp[1,1] + temp[2,1]) / 2;
    centerarr[1,2] = (temp[1,2] + temp[2,2]) / 2;

    # 2
    centerarr[2,1] = (temp[2,1] + temp[3,1]) / 2;
    centerarr[2,2] = (temp[2,2] + temp[3,2]) / 2;

    # 3
    centerarr[3,1] = (temp[3,1] + temp[1,1]) / 2;
    centerarr[3,2] = (temp[3,2] + temp[1,2]) / 2;

    @. temp = centerarr;

    xdiff = [centerarr[1,1] - centerarr[2,1], centerarr[2,1] - centerarr[3,1], centerarr[3,1] - centerarr[1,1]];
    xdiff = @. abs(xdiff);
    ydiff = [centerarr[1,2] - centerarr[2,2], centerarr[2,2] - centerarr[3,2], centerarr[3,2] - centerarr[1,2]];
    ydiff = @. abs(ydiff);

    global maxdiff = maximum(maximum([xdiff, ydiff]));
end

center = [centerarr[1,1], centerarr[1,2]];
@. Nodes[4,:] = center;

#=
│ Row │ x       │ y       │
│     │ Float64 │ Float64 │
├─────┼─────────┼─────────┤
│ 1   │ 1.0     │ 2.0     │
│ 2   │ 6.0     │ 4.0     │
│ 3   │ 3.0     │ 6.0     │
│ 4   │ 3.33333 │ 4.0     │
=#

# =====
# Plotting
linex = [Nodes[1,1], Nodes[2,1], Nodes[3,1], Nodes[1,1]];
liney = [Nodes[1,2], Nodes[2,2], Nodes[3,2], Nodes[1,2]];
# linex = Nodes[:,1];
# liney = Nodes[:,2];
push!(linex, Nodes[1,1]);
push!(liney, Nodes[1,2]);
plot(linex,liney)

mid1x = [Nodes[1,1], Nodes[4,1]];
mid2x = [Nodes[2,1], Nodes[4,1]];
mid3x = [Nodes[3,1], Nodes[4,1]];

mid1y = [Nodes[1,2], Nodes[4,2]];
mid2y = [Nodes[2,2], Nodes[4,2]];
mid3y = [Nodes[3,2], Nodes[4,2]];

plot!(mid1x, mid1y, color = :blue);
plot!(mid2x, mid2y, color = :blue);
plot!(mid3x, mid3y, color = :blue);

scatter!(Nodes[:,1], Nodes[:,2], color=:purple)

# =====
# Next we need to map the elements

#=
Element  Node 1   Node 2   Node 3
1           1       2        4
2           2       3        4
3           1       3        4
=#

N = Nodes;

ElementMap = [
    [1 2 4] ######
    [2 3 4]
    [1 3 4]
];

# Element 1
E1 = [
    [1 N[1,1] N[1,2]]
    [1 N[2,1] N[2,2]]
    [1 N[3,1] N[3,2]]
];
#=
 1.0  1.0  2.0
 1.0  6.0  4.0
 1.0  3.0  6.0
=#

# Element 2
E2 = [
    [1 N[2,1] N[2,2]]
    [1 N[3,1] N[3,2]] 
    [1 N[4,1] N[4,2]]
];
#=
 1.0  6.0      4.0
 1.0  3.0      6.0
 1.0  3.33333  4.0
=#

# Element 3
E3 = [
    [1 N[1,1] N[1,2]]
    [1 N[3,1] N[3,2]] 
    [1 N[4,1] N[4,2]]
];
#=
 1.0  1.0      2.0
 1.0  3.0      6.0
 1.0  3.33333  4.0
=#

iE1 = inv(E1);
iE2 = inv(E2);
iE3 = inv(E3);

@variables x y

# Element 1
N11(x,y) = iE1[1,1] + x*iE1[1,2] + y*iE1[1,3];
N12(x,y) = iE1[2,1] + x*iE1[2,2] + y*iE1[2,3];
N13(x,y) = iE1[3,1] + x*iE1[3,2] + y*iE1[3,3];

# Element 2
N21(x,y) = iE2[1,1] + x*iE2[1,2] + y*iE2[1,3];
N22(x,y) = iE2[2,1] + x*iE2[2,2] + y*iE2[2,3];
N23(x,y) = iE2[3,1] + x*iE2[3,2] + y*iE2[3,3];

# Element 3
N31(x,y) = iE3[1,1] + x*iE3[1,2] + y*iE3[1,3];
N32(x,y) = iE3[2,1] + x*iE3[2,2] + y*iE3[2,3];
N33(x,y) = iE3[3,1] + x*iE3[3,2] + y*iE3[3,3];

# println("N³₁(x,y) = ", N31(x,y))
# println("N³₂(x,y) = ", N32(x,y))
# println("N³₃(x,y) = ", N33(x,y))

# =====
# Time to compute local k

@parameters α β

k1 = Array{Operation}(undef, 3,3);
k2 = Array{Operation}(undef, 3,3);
k3 = Array{Operation}(undef, 3,3);

for i ∈ 1:3
    for j ∈ 1:3
        k1const = α * (iE1[2,i]*iE1[2,j] + iE1[3,i]*iE1[3,j]);
        k1integrand(x,y) = (iE1[1,i] + x*iE1[2,i] + y*iE1[3,i]) * (iE1[1,j] + x*iE1[2,j] + y*iE1[3,j]);

        # hcubature computes the 2-dimensional integral
        temp1, err = hcubature(w -> k1integrand(w[1],w[2]), E1[:,2], E1[:,3]);

        k1[i,j] = k1const + β * temp1;

        # ---

        k2const = α * (iE2[2,i]*iE2[2,j] + iE2[3,i]*iE2[3,j]);
        k2integrand(x,y) = (iE2[1,i] + x*iE2[2,i] + y*iE2[3,i]) * (iE2[1,j] + x*iE2[2,j] + y*iE2[3,j]);

        temp2, err = hcubature(w -> k2integrand(w[1],w[2]), E2[:,2], E2[:,3]);

        k2[i,j] = k2const + β * temp2;

        # ---

        k3const = α * (iE3[2,i]*iE3[2,j] + iE3[3,i]*iE3[3,j]);
        k3integrand(x,y) = (iE3[1,i] + x*iE3[2,i] + y*iE3[3,i]) * (iE3[1,j] + x*iE3[2,j] + y*iE3[3,j]);

        temp3, err = hcubature(w -> k3integrand(w[1],w[2]), E3[:,2], E3[:,3]);

        k3[i,j] = k3const + β * temp3;
    end
end