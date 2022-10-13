# !!! Need to recompute, don't trust numbers
# Interpolation functions for the plane frame element
function N(ℸ, s, hₑ)
    if ℸ == 1
        (1/2) * (1-s);
    elseif ℸ == 2
        (1/4) * (2+s) * (1-s)^2;
    elseif ℸ == 3
        (1/8) * hₑ * (1+s) * (1-s^2);
    elseif ℸ == 4
        (1/2) * (1+s);
    elseif ℸ == 5
        (1/4) * (2-s) * (1+s)^2;
    elseif ℸ == 6
        -(1/8) * hₑ * (1-s) * (1+s)^2;
    end
end

function dN(ℸ, s)
    if ℸ == 1
        -0.5;
    elseif ℸ == 2
        0.25 * (1 + -s) ^ 2 + -0.5 * (1 + -s) * (2 + s);
    elseif ℸ == 3
        0.125 * (1 + -s^2) * hₑ + -0.25 * (1 + s) * hₑ * s;
    elseif ℸ == 4
        0.5;
    elseif ℸ == 5
        -0.25 * (1 + s)^2 + 0.5 * (2 + -s) * (1 + s);
    elseif ℸ == 6
        0.125 * (1 + s)^2 * hₑ + -0.25 * (1 + -s) * (1 + s) * hₑ;
    end
end

function f(i, j, hₑ)
    if j == 1 && (i == 1 || i == 4)
        w = -w1 * sin(θvec[1]); # wx 
    elseif j == 1 && (i != 1 || i != 4)
        w = -w1 * cos(θvec[1]); # wy
    elseif j == 2 && (i == 1 || i == 4)
        w = 0;
    elseif j == 2 && (i != 1 || i != 4)
        w = -w2;
    elseif j == 3
        w = 0;
    end

    func(s) = w * N(i,s, hₑ) * (hₑ/2);

    temp, err = quadgk(s -> func(s), -1, 1);
        # Computes the integral ∫₋₁¹ w Nᵢ(s) hₑ/2 ds
    return temp
end

# Local frame element stiffness
function k(hₑ, θₑ)
    cₑ = cos(θₑ);
    sₑ = sin(θₑ);
    αₑ = (hₑ^2 * Aₑ) / (2 * Iₑ);
    βₑ = (2 * EIₑ) / (hₑ^3);

    pluscs = αₑ*cₑ^2 + 6*sₑ^2;
    plussc = αₑ*sₑ^2 + 6*cₑ^2;
    par = (αₑ-6)*cₑ*sₑ;
    hs = 3*hₑ*sₑ;
    hc = 3*hₑ*cₑ;
    h2 = hₑ^2;
    
    kᵉ = βₑ * 
    [pluscs    par    -hs     -pluscs   -par    -hs
    par     plussc    hc      -par    -plussc   hc
    -hs     hc      2*h2    hs      -hc     h2
    -pluscs   -par    hs      pluscs    par     hs
    -par    -plussc   -hc     par     plussc    -hc
    -hs     hc      h2      hs      -hc     2*h2];

    return kᵉ
end

function bc(i)
    @. Kmod[i,:] = 0;
    @. Kmod[:,i] = 0;
    Kmod[i,i] = 1;
    Fmod[i] = 0;
end