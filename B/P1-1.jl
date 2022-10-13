using DataFrames
using ModelingToolkit

#= 
    (a) Obtain the connectivity matrix [C] 
=#

C = [
    [11 13 3 1 12 8 2 7]        # 1
    [15 17 6 4 16 10 5 9]       # 2
    [22 24 13 11 23 19 12 18]   # 3
    [24 26 15 13 25 20 14 19]   # 4
    [26 28 17 15 27 21 16 20]   # 5
    [33 35 24 22 34 30 23 29]   # 6
    [35 37 26 24 36 31 25 30]   # 7
    [37 39 28 26 38 32 27 31]   # 8
];


#=
│ Row │ x1    │ x2    │ x3    │ x4    │ x5    │ x6    │ x7    │ x8    │
│     │ Int64 │ Int64 │ Int64 │ Int64 │ Int64 │ Int64 │ Int64 │ Int64 │
├─────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┼───────┤
│ 1   │ 11    │ 13    │ 3     │ 1     │ 12    │ 8     │ 2     │ 7     │
│ 2   │ 15    │ 17    │ 6     │ 4     │ 16    │ 10    │ 5     │ 9     │
│ 3   │ 22    │ 24    │ 13    │ 11    │ 23    │ 19    │ 12    │ 18    │
│ 4   │ 24    │ 26    │ 15    │ 13    │ 25    │ 20    │ 14    │ 19    │
│ 5   │ 26    │ 28    │ 17    │ 15    │ 27    │ 21    │ 16    │ 20    │
│ 6   │ 33    │ 35    │ 24    │ 22    │ 34    │ 30    │ 23    │ 29    │
│ 7   │ 35    │ 37    │ 26    │ 24    │ 36    │ 31    │ 25    │ 30    │
│ 8   │ 37    │ 39    │ 28    │ 26    │ 38    │ 32    │ 27    │ 31    │
=#

#=
    (b) Give the following components of the assembled global 
     coefficient matrix expressed in terms of element matrix kᵉᵢⱼ:
        K₂₆,₁₃
        K₂₄,₂₆
        K₃₀,₁₉
        K₂₄,₂₄
        K₃₁,₂₄
=#

#=
Say that:
    kᵢⱼ = ∫∫(∂Nᵢ/∂x ∂Nⱼ/∂x + ∂Nᵢ/∂y ∂Nⱼ/∂y) dx dy

where 
    N(x,y) ⇒ serendipity function
=#

K = Array{Operation}(undef, 39,39);

@parameters k11 k12 k13 k14 k15 k16 k17 k18
@parameters k21 k22 k23 k24 k25 k26 k27 k28
@parameters k31 k32 k33 k34 k35 k36 k37 k38
@parameters k41 k42 k43 k44 k45 k46 k47 k48
@parameters k51 k52 k53 k54 k55 k56 k57 k58
@parameters k61 k62 k63 k64 k65 k66 k67 k68
@parameters k71 k72 k73 k74 k75 k76 k77 k78
@parameters k81 k82 k83 k84 k85 k86 k87 k88

k = [
    [k11 k12 k13 k14 k15 k16 k17 k18]
    [k21 k22 k23 k24 k25 k26 k27 k28]
    [k31 k32 k33 k34 k35 k36 k37 k38]
    [k41 k42 k43 k44 k45 k46 k47 k48]
    [k51 k52 k53 k54 k55 k56 k57 k58]
    [k61 k62 k63 k64 k65 k66 k67 k68]
    [k71 k72 k73 k74 k75 k76 k77 k78]
    [k81 k82 k83 k84 k85 k86 k87 k88]
];

# Klocal = Array{Operation}(undef, 8,38,38);

#= 
Sums up the local k into global k follow the below order:
    1:8
    6:13
    10:17
    14:21
    18:25
    22:29
    26:33
    32:39

I had trouble getting the local matricies to index evenly
into the global of size 39×39
=#
u_index = 0;

@. K[:,:] = ModelingToolkit.Constant(0);

for i ∈ 1:8
    if i == 1
        j = i;
    elseif i == 2
        j = u_index - 2;
    elseif i == 8
        j = 39-7;
    else 
        j = u_index - 3;
    end
    l_index = j;
    global u_index = l_index + 7;
    indr = l_index:u_index;
    println(indr)
    @. K[indr, indr] += k;
end

println("K26,13 = ", K[26,13])
println("K24,26 = ", K[24,26])
println("K30,19 = ", K[30,19])
println("K24,24 = ", K[24,24])
println("K31,24 = ", K[31,24])