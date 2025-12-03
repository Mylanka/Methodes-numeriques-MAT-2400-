function y = problimite(N,P,Q,R,a,b,alpha,beta)

% Entrées:
%   N : nombre de points intérieurs
%   P, Q, R : vecteurs (N+2) contenant p(xi), q(xi), r(xi) aux noeuds
%   i=0,1,...,N+1
%   a, b : bornes de l'intervalle
%   alpha : condition au bord y(a) = alpha
%   beta : condition au bord y(b) = beta
%
% Sortie:
%   y : vecteur (N+2 x 1) solution approchée incluant les bords

% Calcul du pas
h=(b-a)/(N+1);

% Diagonal principal (2+q_i*h^2) pour i=1..N
D = 2 + Q(2:N+1) * h^2;

% Diagonale inférieur (-1-p_i * h/2) pour i=2..N
I = -1 - P(3:N+1) * h/2;

% Diagonale supérieur (-1+p_i*h/2) pour i=1..N-1
S= -1+ P(2:N)*h/2;

%Vecteur b = -r_i*h^2 pour i=1..N
b_vec = -R(2:N+1)*h^2;

%Modification du premier élément -r_1*h^2 + (1+p_1h/2)*alpha
b_vec(1)= b_vec(1)+(1+P(2)*h/2)*alpha;

%Modification du dernier élément -r_Nh^2+(1-p_n*h/2)*beta
b_vec(N)= b_vec(N)+ (1-P(N+1)*h/2)*beta;

% Résolution du système tridiagonal A*y_int = b_vec
y_int = tridiagonal(D, I, S, b_vec);

% Construction du vecteur solution (N+2 points) avec les conditions aux bords
alpha = alpha(:);   
beta  = beta(:);
y_int = y_int(:);   


y = [alpha; y_int; beta];

