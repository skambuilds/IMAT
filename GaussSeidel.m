% Metodo di Gauss-Seidel
% x - soluzione calcolata
% k - numero di iterazioni eseguite
% resvec - vettore che contiene il residuo ad ogni iterazione
function [x,k,resvec] = GaussSeidel(A,b,tau,maxn,x)
% Individuo dimensione del sistema
n=size(A,1);
% Preallocazione risorse
x0=100*ones(n,1);
resvec=zeros(maxn,1);
% Impostazione dello splitting additivo
P = tril(A);
N = P - A;
k=0;
% Ciclo di iterazione con criterio di Cauchy
while(norm(x-x0)>tau*norm(x)) && (k<maxn)
    x0 = x;
    k = k+1;
    % Calcolo del residuo
    r = b-A*x;
    % Memorizzazione del residuo
    resvec(k) = norm(r);
    
    step = N*x0+b;
    x = P\(step);
    
end


