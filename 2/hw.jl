n = 1;

c = AbstractVector{Float64};
φ(x,i) = (x^i)*(1-x);

for i=1:n
    global ut(x) = c .* φ(x,i);
end