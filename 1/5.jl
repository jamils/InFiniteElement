u(x) = -(1/(420*x^2))*(38230 - 24667*x + 105*x^4 + 42*x^5 - 7*x^6 + x^7);

a = 2; # Lower boundary
b = 5; # Upper boundary

n = 5; # Number of quadrature points

n1i = [1,2,3,4,5];
s1i = [0, 1/sqrt(3), sqrt(3/5), sqrt((3/7)+(2/7)*sqrt(6/5)), (1/3)*sqrt(5+2*sqrt(10/7))]; # Gauss-Legendre coordinates

w1i = [2, 1, (5/9), (18-sqrt(30))/36, (322-13*sqrt(70))/900]; # Gauss-Legendre weights

n2i = [2,3,4,5,6];
s2i = [1, 1, 1, sqrt(3/7), sqrt((1/3) + (2 * sqrt(7))/21)]; # Gauss-Lobatto coordinates

w2i = [1, (1/3), (1/6), (49/90), (14-sqrt(7))/30]; # Gauss-Lobatto weights

α = (b-a)/2;



# Computing Gauss-Legendre
γ1 = zeros(n);
x1i = zeros(n);

for j in 1:n
    x1i[j] = (1/2) * ((a+b) + (b-a)*s1i[j]);
    γ1[j] = u(x1i[j]) * w1i[j];
end

A1 = α*sum(γ1);

# Computing Gauss-Lobatto
γ2 = zeros(n);
x2i = zeros(n);

for j in 1:n
    x2i[j] = (1/2) * ((a+b) + ((b-a)*s2i[j]));
    γ2[j] = u(x2i[j]) * w2i[j];
end

A2 = α*sum(γ2);

println(n)
println(A1)
println(A2)

Aact = 5.66748;

e1 = ((A1 - Aact) / Aact) * 100;
e2 = ((A2 - Aact) / Aact) * 100;

println(e1, "%")
println(e2, "%")