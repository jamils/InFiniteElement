# Plotting script
using Plots

title = "Question 2 method comparison";

function infunc()
    include("2.jl")

    return y₊
end

function outfunc()
    include("2.jl")

    return y₊
end

yin = infunc();

yout = outfunc();

xplot = 0:hₑ:(a+b);

yplot = [yin, yout];

labels = ["f(x) in integral" "f(x) outside integral"];

plt = plot(xplot, yplot, label = labels, legend = :bottomleft, title = title, thickness_scaling = 1.3);

display(plt)