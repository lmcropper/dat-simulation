close
clear all

z_range = [-1e-3, 1e-3]; %Axial range in m
r_range = [-500e-6, 500e-6]; %Radial range in m
z_steps = 1000;
r_steps = 1000;

z_lin = linspace(z_range(1), z_range(2), z_steps);
r_lin = linspace(r_range(1), r_range(2), r_steps);

lambda = 633e-9;
n = 1;
u = 3e8 / n;
k = 2*pi/lambda;
t = 0;

a = 12e-3;
f = 25.4e-3;
w_0 = 1.1e-3;
n = 1.5;
P = 1;

u_lin = k * z_lin * power(a/f, 2);      %Normalized axial coordinate
v_lin = k * r_lin * (a/f);         %Normalized radial coordinate

[u, v] = meshgrid(u_lin, v_lin);

J_0 = @(x) integral(@(tau) cos(-x .* sin(tau)), 0, pi);

scale_term = (8*pi*power(a, 4)*P) / (power(lambda, 2)*power(f, 2)*power(w_0, 2));

ab_func = @(x) -(power(a*x, 4)*power(n, 2)) / (8*power(f, 3)*power((n-1), 2));
U=0;
V=0;

I = integral(@(rho) exp(-power(rho, 2)/(w_0/a)) * exp(-1i*(k*ab_func(rho) - ((U*power(rho, 2))/2))) * J_0(V*rho) * rho, 0, 1)




