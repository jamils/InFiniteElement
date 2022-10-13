using LinearAlgebra
using QuadGK

include("func.jl")

C = [[1, 2, 3] [2, 3, 4]];
#=
 1  2
 2  3
 3  4
=#

# =====
# Things for each element
h1 = sqrt(a^2 + h^2);
h2 = b;
h3 = h;

θ1 = atan(h/a);
θ2 = 0;
θ3 = - π/2;

hvec = [h1 h2 h3];
θvec = [θ1 θ2 θ3];
xvec = [0, a, a+b, a+b];
yvec = [0, h, h, 0];

# =====

Q = zeros(dof);
Q[4] = P1;
Q[5] = -P2;
Q[9] = M;

# Tᵉ = [cₑ  sₑ 0 0 0 0
#  -sₑ cₑ 0 0 0 0
#  0 0 1 0 0 0
#  0 0 0 cₑ  sₑ 0
#  0 0 0 -sₑ cₑ 0
#  0 0 0 0 0 1];

Tᵀ(cₑ, sₑ) = 
[cₑ -sₑ 0 0 0 0
 sₑ cₑ 0 0 0 0
 0 0 1 0 0 0
 0 0 0 cₑ -sₑ 0
 0 0 0  sₑ cₑ 0
 0 0 0 0 0 1];

# =====
# Compute some values
K = zeros(dof, dof); # Global stiffness matrix
F = zeros(dof);

# Persistant cache of original values
Korg = zeros(dof, dof);
Forg = zeros(dof);

for i ∈ 1:m
    ftemp = zeros(6);
    for j ∈ 1:6
        ftemp[j] = f(j,i, hvec[i]); # Compute the local f vector
    end
    kmat = k(hvec[i], θvec[i]); # Compute the local k matrix
    # ===
    # This section adds the local k into the elemental Kᵉ matrix
        len = length(kmat[1,:]);
        ib = Int64(  ((dof/4)*(i-1))+1  );
        ie = Int64(  ib + len-1  );
        K[ib:ie,ib:ie] = @. K[ib:ie,ib:ie] + kmat; # Summing Kᵉ
    # ===
    ftemp = Tᵀ(cos(θvec[i]), sin(θvec[i])) * ftemp; # Modifying local f with Tᵀ matrix
    F[ib:ie] = @. F[ib:ie] + ftemp; # Summing Fᵉ
    if i == m
        # Records orginal values in cache
        @. Korg = K;
        @. Forg = F;
    end
end

Kmod = Korg;
Fmod = Forg;

Fmod = @. F + Q; # Include Initial forces into global F

# Boundary Conditions
# !!!
    # I believe the pin shown in problems 2 and 3 are to allow the beams to rotate. In the example (same setup as problem 1) all three DOF's are 0 for each end without any pin
# !!!

# ===
# This determines the boundary conditions determined by the System
if problem == 1
    nodepoints = [1, 2, 3, 10, 11, 12];
elseif problem == 2
    nodepoints = [1, 2, 10, 11, 12];
elseif problem == 3
    nodepoints = [1, 2, 3, 10, 11];
end

# Sets the boundary conditions in the modified K matrix and F vector
for i ∈ nodepoints
    bc(i)
end


q = Kmod \ Fmod; # Computes the displacement vector q
Qnew = Korg*q - Forg; # Computes the Force vector Q from the unmodified K and F values

include("post.jl")