function x = tridiagonal(D, I, S, b)
    %taille matrice tridiagonal
    n = length(D);
    %initialisation vecteur y
    y=zeros(1,n);
    %initialisation vecteur x
    x=zeros(1,n);
    y(1) = b(1)/D(1);
    for i=2:n
    % Factorisation de la matrice tridiagonale
    % Ligne de U
        S(i-1) = S(i-1)/D(i-1);
    % Pivot
        D(i) = D(i) - I(i-1)*S(i-1);
    % Descente triangulaire
        y(i) = (b(i)-I(i-1)*y(i-1))/D(i);
    end
    % Remontée triangulaire
    x(n) = y(n);
    for i=n-1:-1:1
        x(i) = y(i) - S(i)*x(i+1);
    end

