using LinearAlgebra
using Plots
using ModelingToolkit
using HCubature
using QuadGK

# =====
# Input data (to be changed later)

# Node locations

diff1x = √(1^2 + 6^2)/3;
diff1y = √(3^2 + 5^2)/3;

diff2x = √(6^2 + 3^2)/3;
diff2y = √(5^2 + 8^2)/3;

diff3x = √(1^2 + 3^2)/3;
diff3y = √(3^2 + 8^2)/3;

Nodes = [
    [1 3]
    [6 5]
    [3 8]
    [1+diff1x 3+diff1y]
    [6-diff2x 5+diff2y]
    [3-2*diff3x 8-2*diff3y]
    [1+2*diff1x 3+2*diff1y]
    [6-2*diff2x 5+2*diff2y]
    [3-2*diff3x 8-2*diff3y]
    [0.0 0.0]
];

# Now to compute node 10, which is in the center of the triangle

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
@. Nodes[10,:] = center;

#=
│ Row │ x       │ y       │
│     │ Float64 │ Float64 │
├─────┼─────────┼─────────┤
│ 1   │ 1.0     │ 2.0     │
│ 2   │ 6.0     │ 4.0     │
│ 3   │ 3.0     │ 6.0     │
│ 4   │ 3.33333 │ 4.0     │
=#

N = Nodes;

E = [
    [1 N[1,1] N[1,2] N[1,1]*N[1,2]]
    [1 N[2,1] N[2,2] N[2,1]*N[2,2]]
    [1 N[3,1] N[3,2] N[3,1]*N[3,2]]
    [1 N[4,1] N[4,2] N[4,1]*N[4,2]]
];

iE = inv(E);

@variables x y

N = Array{Operation}(undef, 4);

@. N = iE[1,:] + x*iE[2,:] + y*iE[3,:] + x*y*iE[4,:];

N = @. simplify(N);

# for i ∈ 1:4 println("N$i(x,y) = ", N[i]) end

@derivatives Dx'~x Dy'~y

dNx = Array{Operation}(undef, 4);
dNy = Array{Operation}(undef, 4);

@. dNx = expand_derivatives(Dx(N[:]));
@. dNy = expand_derivatives(Dy(N[:]));

# @. println("dNx = ", dNx)

Γ = Array{Operation}(undef, 3);

xr = E[:,2];
yr = E[:,3];
Γ[1] = ( (xr[2]*yr[1] - xr[1]*yr[2])/(xr[2]-xr[1]) ) + ( (yr[2]-yr[1])/(xr[2]-xr[1]) )*x;
Γ[2] = ( (xr[3]*yr[2] - xr[2]*yr[3])/(xr[3]-xr[2]) ) + ( (yr[3]-yr[2])/(xr[3]-xr[2]) )*x;
Γ[3] = ( (xr[1]*yr[3] - xr[3]*yr[1])/(xr[1]-xr[3]) ) + ( (yr[1]-yr[3])/(xr[1]-xr[3]) )*x;

γ = zeros(3);

γ[1] = √((Nodes[2,1]-Nodes[1,1])^2 + (Nodes[2,2]-Nodes[1,2])^2);
γ[2] = √((Nodes[3,1]-Nodes[2,1])^2 + (Nodes[3,2]-Nodes[2,2])^2);
γ[3] = √((Nodes[3,1]-Nodes[1,1])^2 + (Nodes[3,2]-Nodes[1,2])^2);

# Computing the transformation coordinates for (x,y)→(s,t)
a = γ[1];
s = sum(γ)/2;
Area = √(s*(s-γ[1])*(s-γ[2])*(s-γ[3]));
b = 2*Area/γ[1];
c = sin(acos(b/γ[3]))*γ[3];

function ∫(Nint)
    # evalfunc = eval(build_function(Nint, x,y));
    # func1(x,y) = evalfunc(x,y);
    func1(x,y) = Nint;
    func2(x,y) = func1(y*x,y^2/2);
    val = func2(0,0);
    func3(x,y) = func2(x,y) - val + (val*y);

    Γ1(x) = Γ[1];
    Γ2(x) = Γ[2];
    Γ3(x) = Γ[3];

    out1(x) = func3(x,Γ1(x));
    out2(x) = func3(x,Γ2(x));
    out3(x) = func3(x,Γ3(x));

    return out1(x), out2(x), out3(x)
end

@parameters α β

xmin = minimum(Nodes[:,1]);
xmax = maximum(Nodes[:,1]);
ymin = minimum(Nodes[:,2]);
ymax = maximum(Nodes[:,2]);

K = Array{Operation}(undef, 4,4);
k00 = zeros(4,4);
k10 = zeros(4,4);
k01 = zeros(4,4);

i = 1; j = 2;

# @time for i ∈ 1:4
    # for j ∈ 1:4
        dNxi = eval(build_function(dNx[i], y));
        dNxj = eval(build_function(dNx[j], y));
        dNyi = eval(build_function(dNy[i], x));
        dNyj = eval(build_function(dNy[j], x));

        Ni = eval(build_function(N[i], [x,y]));
        Nj = eval(build_function(N[j], [x,y]));

        ftempx(x) = dNxi(x) * dNxj(x);
        ftempy(y) = dNyi(y) * dNyj(y);
        ftempxy(x,y) = Ni([x,y]) * Nj([x,y]);

        # println(ftempx(x))

##############
# TEST SECTION
##############
TESTfunc1(x,y) = ftempxy(x,y);
TESTfunc2(x,y) = TESTfunc1(y*x,y^2/2);
TESTval = TESTfunc2(0,0);
TESTfunc3(x,y) = TESTfunc2(x,y) - TESTval + (TESTval*y);

TESTΓ1(x) = Γ[1];
TESTΓ2(x) = Γ[2];
TESTΓ3(x) = Γ[3];

TESTout1(x) = TESTfunc3(x,TESTΓ1(x));
TESTout2(x) = TESTfunc3(x,TESTΓ2(x));
TESTout3(x) = TESTfunc3(x,TESTΓ3(x));
##############
##############

        # fint1(x), fint2(x), fint3(x) = ∫(ftempxy(x,y));

        # fint(x) = eval(build_function(tempint(x,Γ1(x))));

        tempx, err = quadgk(x -> ftempx(x), xmin, xmax);
        tempy, err = quadgk(y -> ftempy(y), ymin, ymax);
        # temp1, err = quadgk(x -> fint(x), xmin, xmax);
        # temp1, err = quadgk(x -> fint1(x), xmin, xmax);
        # temp2, err = quadgk(x -> fint(x, Γ[2]), xmin, xmax);
        # temp3, err = quadgk(x -> fint(x, Γ[3]), xmin, xmax);
        # tempxy, err = hcubature(w -> ftempxy(w[1],w[2]), [xmin,ymin], [xmax,ymax]);

        tempxy = temp1 + temp2 + temp3;

        k00[i,j] = tempxy;
        k10[i,j] = tempx;
        k01[i,j] = tempy;

        K[i,j] = α * (tempx + tempy) + β * tempxy;
    # end
# end