clear; close all; clc;

N = [30, 50];
a = 0.9; 
b = 1;
alpha = 0; 
beta = 1;

P = @(x) -1./x;
Q = @(x) 0.*x;
R = @(x) -1.6./x.^4;

% Pour N = 30, h=1/30
N_30 = N(1);
x_30 = linspace(a, b, N_30+2)'; 
P_vec30 = P(x_30);
Q_vec30 = Q(x_30);
R_vec30 = R(x_30);

y_30 = problimite(N_30, P_vec30, Q_vec30, R_vec30, a, b, alpha, beta);

%h=1/50
N_50   = N(2);
x_50   = linspace(a,b,N_50+2)';
P_50   = P(x_50);
Q_50 = Q(x_50); 
R_50 = R(x_50);
y_50  = problimite(N_50,P_50,Q_50,R_50,a,b,alpha,beta);

%Solution exacte
syms y(x)
ode  = diff(y,x,2) + diff(y,x)/x + 1.6/x^4 == 0;
cond1 = y(0.9) == 0;
cond2 = y(1)   == 1;
y_exact = dsolve(ode, [cond1 cond2]);
x_fine = 0.9:0.001:1;
y_fine = double(subs(y_exact,x,x_fine));

figure(1);
plot(x_fine,y_fine,'b-', 'LineWidth', 2);
hold on;
plot(x_30,y_30,'-o', x_50,y_50,'-x');
legend('exacte','N=30','N=50');
xlabel('x');
ylabel('y(x)');
title('Comparaison solution exacte et approchées par différences finies');
grid on;
hold off;

% Pour Figure 2 - erreur constante
d = 0.81;
c = 1.4;
y_exact_num = @(x) (c - 0.4./x.^2) - (c - 0.4/d) * log(x) / log(0.9);

h_list = [1e-2, 1e-3, 1e-4, 1e-5];
err_list = [];

for h = h_list
    N = round((b-a)/h) - 1;
    x = linspace(a, b, N+2)';
    
    y_approx = problimite(N, P(x), Q(x), R(x), a, b, alpha, beta);
    
    err = max(abs(y_approx(2:end-1) - y_exact_num(x(2:end-1))));
    err_list = [err_list, err];
end

figure(2);
loglog(h_list, err_list, 'bo-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on; grid on;
loglog(h_list, err_list(1)*(h_list/h_list(1)).^2, 'r--', 'LineWidth', 1.5);
xlabel('h'); ylabel('E(h)');
title('Erreur E(h) en fonction de h (log-log)');
legend('E(h)', 'Pente 2', 'Location', 'best');
hold off;










