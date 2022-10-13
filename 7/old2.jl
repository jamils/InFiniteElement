using LinearAlgebra
using QuadGK

# Solving for the second problem

#=
    This problem gives a uniform bar of length a+b
    with a non-uniform downward force f(x). There
    are three marked locations A, B, and C. A is 
    pinned to a wall at x=0. B is on a roller, at
    x=a. C is the very end of the bar with no
    attachment. 

    For this, we will use m=5 elements, splitting
    such at there are three elements in AB and two
    in BC. With the given lengths a=6m and b=4m,
    each element is 2m in length. We will treat this
    as a three node system.
=#

# Keeping all units in N and m

# ==========
# Parameters given in the problem statement

a = 6; # (m) Length of AB
b = 4; # (m) Length of BC
L = a+b; # (m) Total lengh of the bar
EI = 5e6; # (N⋅m²) Axial stiffness is uniform

f(x) = 20*(1+3*((a+b-x)/(a+b))^2) * 1e3; # (N/m)

λ = 3; # Number of nodes per element
m = 5; # Number of elements
hₑ = 1/m; # Step size
n = 2*λ*m; # Degrees of freedom, two for every node in the system

xᵥ = 0:hₑ:L; # Array of xᵢ node points

# ==========
# u = zeros(n);      # The solution of the differential equation
K = zeros(n, n);   # The global K vector
β = zeros(n);      # The global b vector
                   # Using β to not conflict with other symbol b

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

# ==========
# Solving along each λ-node elements
for ℵ ∈ 0:(m-1)     # Indicie of node, essential the i for xᵢ
    xs(s) = (1/2) * ((2 * xᵥ[begin+ℵ]) + hₑ * (1+s));
    for ℶ ∈ 1:λ     # Looping through each node
        for ℷ ∈ 1:λ     # Looping through each node, again
            k(s) = (2/hₑ) * xs(s); #(dN(ℶ, s) * dN(ℷ, s)) ;
            ktemp, err = quadgk(s -> -k(s), -1, 1);
            # Solves the integral ∫(2/hₑ) (dNᵢ/ds dNⱼ/ds) ds from -1 to 1

            if ℵ == 0
                K[ℶ, ℷ] += ktemp * EI; # Sum the above to get Kᵢⱼ
            else
                K[(ℵ*(λ-1))+ℶ, (ℵ*(λ-1))+ℷ] += ktemp;
            end
        end

        β_int(s) = f(xs(s)) * N(ℶ, s) * (hₑ/2);
        βtemp, err = quadgk(s -> -β_int(s), -1, 1);
        # Computes the integral ∫f(s) Nᵢ (hₑ/2)ds from -1 to 1

        if ℵ == 0
            β[ℶ] += βtemp;
        else
            β[(ℵ*(λ-1))+ℶ] += βtemp;
        end
    end
end

error("Stop!")

# Q = zeros(n); # Nodal force vector
q = zeros(n); # Nodal displacement vector

#=
A          B       C
---!---!---!---!---!
Each ! is an element, each - a node
A
    Element 1     x  y
        node 1 -> q1  q2
        node 2 -> q3  q4
        node 3 -> q5  q6
    Element 2 
        node 4 -> q7  q8
        node 5 -> q9  q10
        node 6 -> q11 q12
    Element 3
        node 7 -> q13 114
        node 8 -> q15 q16
        node 9 -> q17 q18
B
    Element 4
        node 10 -> q19 q20
        node 11 -> q21 q22
        node 12 -> q23 q24
    Element 5
        node 13 -> q25 q26
        node 14 -> q27 q28
        node 15 -> q29 q30
C

q1 horizontal, q2 vertical
q1, q2 at A
q19, q20 at B
q29, q30 at C

=#


# Modify K to fit boundary conditions for q 
# q1
K[1,:] .= 0;
K[:,1] .= 0;
K[1,1] = 1;
# q2
K[2,:] .= 0;
K[:,2] .= 0;
K[2,2] = 1;

# q20
K[20,:] .= 0;
K[:,20] .= 0;
K[20,20] = 1;

# Sovling for displacement q
q = K \ β;

q[1] = 0;
q[2] = 0;
q[20] = 0;