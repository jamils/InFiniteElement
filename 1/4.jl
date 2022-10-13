u(x) = -(1/(420*x^2))*(38230 - 24667*x + 105*x^4 + 42*x^5 - 7*x^6 + x^7);

a = 2; # Lower boundary
b = 5; # Upper boundary

n = 16; # Number of sub-areas

println("n = ", n)

k = range(a, length=(n+1), stop=b);

uj = zeros(n+1);
hj = zeros(n);

for j in 1:n
    hj[j] = (1/2) * ( u(k[j]) + u(k[j+1]) );
end

Δx = (b-a) / n;

A = sum(hj.*Δx);

println("Approximate area = ", A)

Aact = 5.66748;

e = ((A - Aact) / Aact) * 100;

println("Error = ", e, "%")