% Classic preconditioned gradient method
% x - solution (at step k)
% k - step (number of iterations performed)
% resvec - vector that contains the residual at each iteration

function [x,k,resvec] = SelfPreGradient(A,b,tau,maxn,Rt,R,x)
% System size
n = size(A,1);

% Initialize the vector x0
x0 = 100*ones(n,1);

% Initialize the residual
r = b - A*x;

% Resolution of a linear system
y = Rt\r;
z = R\y;

% Iterations counter
k = 0;

% Starting the algorithm cycle with control condition. You can choose between:
% - residue control: norm(r) > tau*norm(b)
% - condizione di Cauchy: norm(x-x0) > tau*norm(x)
% I set however a maximum number of iterations
% Preallocation of resources for the residual vector
resvec=zeros(maxn,1);

while(norm(x-x0) > tau*norm(x)) && (k<maxn)
    x0 = x;
    k = k+1;

    % Memorize the norm of the residual in the vector
    resvec(k) = norm(r);

    % Optimize by calculating the matrix-vector product only once
    s = A*z; 

    % Step calculation
    alpha = (z'*r)/(z'*s);

    % New solution
    x = x0+alpha*z;

    % Residual update
    r = r-alpha*s;

    % I solve the system Pz = r per k+1
    y = R'\r;
    z = R\y;
end
