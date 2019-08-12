% Classic gradient method
% x - solution (at step k)
% k - step (number of iterations performed)
% resvec - vector that contains the residual at each iteration

function [x,k,resvec] = SelfGradient(A,b,tau,maxn,x)
% System size
n = size(A,1);

% Initialize the vector x0
x0 = 100*ones(n,1);

% Initialize the residue
r = b - A*x;

% Iterations counter
k = 0;

% Starting the algorithm cycle with control condition. You can choose between:
% - residual control: norm(r) > tau*norm(b)
% - Cauchy condition: norm(x-x0) >t au*norm(x)
% I set however a maximum number of iterations
% Preallocation of resources for the residual vector
resvec=zeros(maxn,1);

while(norm(x-x0) > tau*norm(x)) && (k<maxn)
    x0 = x;
    k = k+1;

    % % Storing the residue
    resvec(k) = norm(r);

    % Optimize by calculating the matrix-vector product only once
    s = A*r; 

    % Step calculation
    alpha = (r'*r)/(r'*s);

    % New solution
    x = x0+alpha*r;

    % Residual update
    r = r-alpha*s; 
end
