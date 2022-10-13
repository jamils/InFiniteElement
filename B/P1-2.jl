using DataFrames
using ModelingToolkit
using HCubature

Nodes = [
    [1 1]
    [1 4]
    [5 1]
    [5 4]
];

@variables x y

#=
First, we take the needed pieces of the
  interpolation functions for a quadtratic
  two-dimensional serendipity element:
    ∂N₁/∂x
    ∂N₂/∂x
    ∂N₁/∂y
    ∂N₂/∂y
=#

dNx1(ξ,η) = (1/4) * (1-η) * (2*ξ+η);
dNx2(ξ,η) = (1/4) * (1-η) * (2*ξ-η);

dNy1(ξ,η) = (1/4) * (1-ξ) * (ξ+2*η);
dNy2(ξ,η) = -(1/4) * (1+ξ) * (ξ-2*η);

#=
Now we take the double integral kᵢⱼ = ∫∫(∂N₁/∂x ∂N₂/∂x + ∂N₁/∂y ∂N₂/∂y) dx dy
  with following bounds of integration:
    x ∈ [1,5]
    y ∈ [1,4]
=#

int(x,y) = dNx1(x,y)*dNx2(x,y) + dNy1(x,y)*dNy2(x,y);
k12, err = hcubature(w -> int(w[1],w[2]), [1,1], [5,4]);

#=
Next, we take the needed pieces of the
  interpolation functions for a quadtratic
  two-dimensional serendipity element:
    N₁(x,y)
    N₂(x,y)
=#

N1(ξ,η) = -(1/4) * (1-ξ) * (1-η) * (1+ξ+η);
N2(ξ,η) = -(1/4) * (1+ξ) * (1-η) * (1-ξ+η);

#=
And finally, we take the double integral kᵢⱼ = ∮ N₁N₂ ds = ∫∫N₁N₂ dx dy
  with following bounds of integration:
    x ∈ [1,5]
    y ∈ [1,4]
=#

hatint(x,y) = N1(x,y) * N2(x,y);
k12hat, err = hcubature(w -> hatint(w[1],w[2]), [1,1], [5,4]);

#=
This gives us the resulting values for k₁₂ and k₁₂_hat
    k₁₂     = -18.125 = -145/8
    k₁₂_hat = -58.2   = -291/5
=#

