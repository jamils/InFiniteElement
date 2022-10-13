problem = 1;

λ = 4; # Number of nodes in the system
m = 3; # Number of elements in the system
dof = 3*λ; # Total degrees of freedom for the system

a = 5; # (m)
b = 2.5; # (m)
h = 4; # (m)
Aₑ = 11e-3; # (m²)
Iₑ = 198e-6; # (m⁴)
Eₑ = 200e9; # (Pa)
EIₑ = Eₑ*Iₑ; # Bending stiffness
EAₑ = Eₑ*Aₑ; # Axial stiffness
w1 = 10e3; # (N/m) Only in -y direction
w2 = 20e3; # (N/m)
P1 = 100e3; # (N)  x direction
P2 = 120e3; # (N) -y direction
# M = 50e3; # (N⋅m)
M = 0;

include("main.jl")