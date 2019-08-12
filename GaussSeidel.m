% Gauss-Seidel method
% x - solution (at step k)
% k - step (number of iterations performed)
% resvec - vector that contains the residual at each iteration

function [x,k,resvec] = GaussSeidel(A,b,tau,maxn,x)
% System size
n = size(A,1);

% Resources pre-allocation
x0 = 100*ones(n,1);
resvec = zeros(maxn,1);

% Setting the additive splitting
P = tril(A);
N = P - A;
k = 0;

% Iteration cycle with Cauchy criterion
while(norm(x-x0) > tau*norm(x)) && (k<maxn)
    x0 = x;
    k = k+1;

    % Calculation of the residual
    r = b-A*x;

    % Storing the residual
    resvec(k) = norm(r);
    
    step = N*x0+b;
    x = P\(step);
end
