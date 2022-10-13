using LinearAlgebra
using Plots
using ModelingToolkit
using HCubature
using DataFrames

# Node locations

a = ((4-1)/2) + 1;
b = ((3-1)/2) + 1;

Nodes = [
    [1 1]
    [4 1]
    [4 3]
    [1 3]
    [a 1]
    [4 b]
    [a 3]
    [1 b]
    [a b]
];

# scatter(Nodes[:,1], Nodes[:,2])

nodex = Nodes[:,1];
nodey = Nodes[:,2];

linex1 = [nodex[1]; nodex[2]; nodex[3]; nodex[4]; nodex[1]];
linex2 = [nodex[1]; nodex[3]];
linex3 = [nodex[2]; nodex[4]];
linex4 = [nodex[6]; nodex[8]];
linex5 = [nodex[5]; nodex[7]];

liney1 = [nodey[1]; nodey[2]; nodey[3]; nodey[4]; nodey[1]];
liney2 = [nodey[1]; nodey[3]];
liney3 = [nodey[2]; nodey[4]];
liney4 = [nodey[6]; nodey[8]];
liney5 = [nodey[5]; nodey[7]];

# plot!(linex1, liney1)#, color=:black)
# plot!(linex2, liney2)#, color=:black)
# plot!(linex3, liney3)#, color=:black)
# plot!(linex4, liney4)#, color=:black)
# plot!(linex5, liney5)#, color=:black)

NodeAltx = [
    Nodes[1,1],
    Nodes[5,1],
    Nodes[2,1],
    Nodes[6,1],
    Nodes[3,1],
    Nodes[7,1],
    Nodes[4,1],
    Nodes[8,1]
];

NodeAlty = [
    Nodes[1,2],
    Nodes[5,2],
    Nodes[2,2],
    Nodes[6,2],
    Nodes[3,2],
    Nodes[7,2],
    Nodes[4,2],
    Nodes[8,2]
];

NodeNums = [1 5 2 6 3 7 4 8 9];

# 8 triangular elements
ElementMap = zeros(8,3);

for i ∈ 1:7
    ElementMap[i, 1] = NodeNums[i];
    ElementMap[i, 2] = NodeNums[i+1];
    ElementMap[i,end] = NodeNums[end];
end

ElementMap[8, 1] = NodeNums[end-1]
ElementMap[8, 2] = NodeNums[begin]
ElementMap[8,end]= NodeNums[end];

#=
│ Row │ x1      │ x2      │ x3      │
│     │ Float64 │ Float64 │ Float64 │
├─────┼─────────┼─────────┼─────────┤
│ 1   │ 1.0     │ 5.0     │ 9.0     │
│ 2   │ 5.0     │ 2.0     │ 9.0     │
│ 3   │ 2.0     │ 6.0     │ 9.0     │
│ 4   │ 6.0     │ 3.0     │ 9.0     │
│ 5   │ 3.0     │ 7.0     │ 9.0     │
│ 6   │ 7.0     │ 4.0     │ 9.0     │
│ 7   │ 4.0     │ 8.0     │ 9.0     │
│ 8   │ 8.0     │ 1.0     │ 9.0     │
=#

# Now to compute things at each element

K = Array{Operation}(undef,8,3,3);

for i ∈ 1:8
    if i<8
        E = [
            [1 NodeAltx[i] NodeAlty[i]]
            [1 NodeAltx[i+1] NodeAlty[i+1]]
            [1 Nodes[9,1] Nodes[9,2]]
        ];
    elseif i == 8
        E = [
            [1 NodeAltx[end] NodeAlty[end]]
            [1 NodeAltx[begin] NodeAlty[begin]]
            [1 Nodes[9,1] Nodes[9,2]]
        ];
    end

    iE = inv(E);

    @variables x y

    N = Array{Operation}(undef, 3);

    for j ∈ 1:3
        N[j] = iE[j,1] + x*iE[j,2] + y*iE[j,3];
        # println("N^$i", "_$j = ", N[j])
    end

    @parameters α β

    k = Array{Operation}(undef, 3,3);

    for ℵ ∈ 1:3
        for ℶ ∈ 1:3
            kconst = α * (iE[2,ℵ]*iE[2,ℶ] + iE[3,ℵ]*iE[3,ℶ]);
            kintegrand(x,y) = (iE[1,ℵ] + x*iE[2,ℵ] + y*iE[3,ℵ]) * (iE[1,ℶ] + x*iE[2,ℶ] + y*iE[3,ℶ]);

            # hcubature computes the 2-dimensional integral
            temp, err = hcubature(w -> kintegrand(w[1],w[2]), E[:,2], E[:,3]; maxevals=10000);

            k[ℵ,ℶ] = kconst + β * temp;
        end
    end
    @. K[i,:,:] = k;

    # print("k^$i = ")
    # display(k);
    # println("\n")
end

