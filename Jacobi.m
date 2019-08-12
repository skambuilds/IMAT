% Jacobi method
% x - solution (at step k)
% k - step (number of iterations performed)
% resvec - vector that contains the residual at each iteration

function [x,k,resvec] = Jacobi(A,b,tau,maxn,x)
% System size
n = size(A,1);

% Pre-allocation of resources
x0 = 100*ones(n,1);
resvec = zeros(maxn,1);

% Setting the additive splitting
d = diag(A);
EF = A - diag(d);
k = 0;

% Iteration cycle with Cauchy criterion
while(norm(x-x0) > tau*norm(x)) && (k<maxn)
    x0 = x;
    k = k+1;
    r = b-A*x; % Calculation of the residue
    resvec(k) = norm(r); % Storage of the residual
   
    step = b - EF*x0;
    x = (step)./d;    
end
