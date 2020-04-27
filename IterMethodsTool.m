fprintf('---IMAT (Iterative methods analysis tool)---\n');
fprintf('Press space to start\n');
pause
clear

% List of matrices
mat1 = '[1] - Diagonally dominant matrix sparse';
mat2 = '[2] - Positive symmetric definite matrix';
mat3 = '[3] - Positive and sparse symmetric definite matrix';
mat4 = '[4] - Poisson''s matrix';

% Choice of matrix type
matrix_str = sprintf('%s\n%s\n%s\n%s\n%s\n',mat1,mat2,mat3,mat4);
fprintf(matrix_str);
choice_str = 'Choose the matrix to use ';
mtype = input(choice_str);

steps_str = input('Enter the dimension steps you want to test separated by spaces or commas: ', 's');
dimens = str2num(steps_str);
steps = size(dimens,2);

% Parameters for performing the analysis
tol = 1e-12;
nmax = 200;
omega = 0.5;

fprintf('List of settings used for performing the test:\n')
fprintf('%d - Tolerance\n', tol);
fprintf('%d - Maximum number of iterations\n', nmax);
fprintf('%d - Relaxation parameter\n', omega);
prompt_3 = sprintf('Do you want to change these values? [Y][N]');

r = input(prompt_3, 's');
if isempty(r)
    r = 'N';
end
if strcmp(r,'y') || strcmp(r,'Y')
    prompt_n = 'Enter the tolerance value: ';
    tol = input(prompt_n);
    prompt_n = 'Enter the maximum number of iterations: ';
    nmax = input(prompt_n);
    prompt_n = 'Enter the omega value for relaxation operations: ';
    omega = input(prompt_n);
end

% We use an incomplete preconditioner - Cholesky
% Required parameters:
opts.type = 'ict'; % Select Cholesky with threshold dropping (ICT)
opts.droptol = 1e-3; % Below which threshold values must be set to zero (small numbers become zero)
% If you increase this value becomes more sparse preconditioner
opts.shape = 'upper'; % choose the upper (RT)(R) format
% opts is actually a struct that contains the parameters

% Initialization of variables
time1 = zeros(steps,1);
time1p = zeros(steps,1);
time2 = zeros(steps,1);
time2p = zeros(steps,1);
time3 = zeros(steps,1);
time3p = zeros(steps,1);
time4 = zeros(steps,1);
time5 = zeros(steps,1);
time6 = zeros(steps,1);

x1 = cell(1,steps);
iter1 = zeros(steps,1);
resvec1 = cell(1,steps);
error1 = zeros(steps,1);

x1p = cell(1,steps);
iter1p = zeros(steps,1);
resvec1p = cell(1,steps);
error1p = zeros(steps,1);

x2 = cell(1,steps);
iter2 = zeros(steps,1);
resvec2 = cell(1,steps);
error2 = zeros(steps,1);

x2p = cell(1,steps);
iter2p = zeros(steps,1);
resvec2p = cell(1,steps);
error2p = zeros(steps,1);

x3 = cell(1,steps);
iter3 = zeros(steps,1);
resvec3 = cell(1,steps);
error3 = zeros(steps,1);

x3p = cell(1,steps);
iter3p = zeros(steps,1);
resvec3p = cell(1,steps);
error3p = zeros(steps,1);

x4 = cell(1,steps);
iter4 = zeros(steps,1);
resvec4 = cell(1,steps);
error4 = zeros(steps,1);

x5 = cell(1,steps);
iter5 = zeros(steps,1);
resvec5 = cell(1,steps);
error5 = zeros(steps,1);

x6 = cell(1,steps);
iter6 = zeros(steps,1);
resvec6 = cell(1,steps);
error6 = zeros(steps,1);

for i=1:steps
    
    fprintf ('Iteration [%d] with dimensioning [%d]\n', i, dimens(i));
    A = MatrixCreator(mtype,dimens(i));
    figure(i);
    spy(A); % In blue the non-zero elements of the selected matrix
    titleTimeBar = sprintf('Sparsity pattern of the matrix with dim= %d', dimens(i));
    title(titleTimeBar)   
    
    % Construction of the linear system to be solved
    n = size (A, 1);
    % System setting, with solution equal to a vector of 1:
    sol = ones (n, 1);
    b = A * sol;
    % Alternatively we choose the known term with a random vector of size n:
    % b = rand (n, 1);
    x0 = zeros (n, 1);
    
    % Matlab gradient conjugate
    tic
    [x, flag, relres, iter, resvec] = pcg (A, b, tol, nmax, [], [], x0); % we can also not insert tol and nmax because they are set by default
    time1(i) = toc;
    error1(i) = norm(sol-x)/norm(sol); % relative error
    x1(i) = {x};
    iter1(i) = iter;
    resvec1(i) = {resvec};
    
    % Evaluate whether to perform the incomplete LU according to the matrices chosen
    % setup.type = 'crout';
    % setup.milu = 'row';
    % setup.droptol = 1e-3;
    % setup.udiag = 1;
    % [L,U] = ilu(A,setup);
    
    % Matlab conjugate gradient preconditioned
    R = ichol(A, opts);
    tic
    [x, flag, relres, iter, resvec] = pcg(A,b,tol,nmax,R',R,x0);
    time1p(i) = toc;
    error1p(i) = norm(sol-x)/norm(sol); % relative error
    x1p(i) = {x};
    iter1p(i) = iter;
    resvec1p(i) = {resvec};
    
    % Method of steepest descent or stationary-phase method or saddle-point method
    tic
    [x, iter, resvec] = SelfGradient(A,b,tol,nmax,x0);
    time2(i) = toc;
    error2(i) = norm(sol-x)/norm(sol); % relative error
    x2(i) = {x};
    iter2(i) = iter;
    resvec2(i) = {resvec};
    
    % Classic preconditioned gradient
    tic
    [x, iter, resvec] = SelfPreGradient(A,b,tol,nmax,R',R,x0);
    time2p(i) = toc;
    error2p(i) = norm(sol-x)/norm(sol); % relative error
    x2p(i) = {x};
    iter2p(i) = iter;
    resvec2p(i) = {resvec};
    
    % Conjugate gradient
    tic
    [x, iter, resvec] = SelfConiugGradient(A,b,tol,nmax,x0); % we can also not insert tol and nmax because they are set by default
    time3(i) = toc;
    error3(i) = norm(sol-x)/norm(sol); % relative error
    x3(i) = {x};
    iter3(i) = iter;
    resvec3(i) = {resvec};
    
    % Preconditioned conjugate gradient
    tic
    [x, iter, resvec] = SelfPreConiugGradient(A,b,tol,nmax,R',R,x0);
    time3p(i) = toc;
    error3p(i) = norm(sol-x)/norm(sol); % relative error
    x3p(i) = {x};
    iter3p(i) = iter;
    resvec3p(i) = {resvec};
    
    % Jacobi method
    tic
    [x, iter, resvec] = Jacobi(A,b,tol,nmax,x0);
    time4(i) = toc;
    error4(i) = norm(sol-x)/norm(sol); % relative error
    x4(i) = {x};
    iter4(i) = iter;
    resvec4(i) = {resvec};
    
    % Gauss Seidel method
    tic
    [x, iter, resvec] = GaussSeidel(A,b,tol,nmax,x0);
    time5(i) = toc;
    error5(i) = norm(sol-x)/norm(sol); % relative error
    x5(i) = {x};
    iter5(i) = iter;
    resvec5(i) = {resvec};   

end

timeMatrix = [time1 time1p time2 time2p time3 time3p time4 time5];
minTime = min(min(timeMatrix));
maxTime = max(max(timeMatrix));

T.Dimensione = dimens';
T.MCG = time1;
T.MPCG = time1p;
T.SG = time2;
T.SPG = time2p;
T.SCG = time3;
T.SPCG = time3p;
T.Jacobi = time4;
T.GaussSeidel = time5;

% Time table
timeTable = struct2table(T)

I.Dimensione = dimens';
I.MCG = iter1;
I.MPCG = iter1p;
I.SG = iter2;
I.SPG = iter2p;
I.SCG = iter3;
I.SPCG = iter3p;
I.Jacobi = iter4;
I.GaussSeidel = iter5;

% Iteration table
iterTable = struct2table(I)

E.Dimensione = (dimens)';
E.MCG = error1;
E.MPCG = error1p;
E.SG = error2;
E.SPG = error2p;
E.SCG = error3;
E.SPCG = error3p;
E.Jacobi = error4;
E.GaussSeidel = error5;

% Errors table
errorTable = struct2table(E)

% Residual charts
figure(steps+1)
for i=1:steps

subplot(1,steps,i)
semilogy(resvec1{i}, '--g', 'LineWidth',5)
hold on
semilogy(resvec1p{i}, '-g', 'LineWidth',5)
hold on
semilogy(resvec2{i}, '--r', 'LineWidth',3)
hold on
semilogy(resvec2p{i}, '-r', 'LineWidth',3)
hold on
semilogy(resvec3{i}, '--b', 'LineWidth',2)
hold on
semilogy(resvec3p{i}, '-b', 'LineWidth',2)
hold on
semilogy(resvec4{i}, '-y', 'LineWidth',2)
hold on
semilogy(resvec5{i}, '-m', 'LineWidth',2)
hold off
legend('MCG','MPCG','SG','SPG','SCG','SPCG','Jacobi','GaussSeidel')
titleGradPlot = sprintf('Residual Iterative Methods, DIM = %d', dimens(i));
title(titleGradPlot)
end

% Bar charts (iterations and times)
figure(steps+2)
for i=1:steps
    iterData = [iter1(i) iter1p(i) iter2(i) iter2p(i) iter3(i) iter3p(i) iter4(i) iter5(i)];
    ax1 = subplot(2,steps,i);
    c = categorical({'MCG','MPCG','SG','SPG','SCG','SPCG','Jacobi','GaussSeidel'});
    c = reordercats(c,{'MCG','MPCG','SG','SPG','SCG','SPCG','Jacobi','GaussSeidel'});
    b = bar(ax1,c,iterData);
    b.FaceColor = 'flat';
    b.CData(1,:) = [0,1,0];
    b.CData(2,:) = [0,1,0];
    b.CData(3,:) = [1,0,0];
    b.CData(4,:) = [1,0,0];
    b.CData(5,:) = [0,0,1];
    b.CData(6,:) = [0,0,1];
    b.CData(7,:) = [1,1,0];
    b.CData(8,:) = [1,0,1];
    ylim([0 nmax])
    titleIterBar = sprintf('Iteractions, DIM = %d', dimens(i));
    title(titleIterBar)
    
    timeData = [time1(i) time1p(i) time2(i) time2p(i) time3(i) time3p(i) time4(i) time5(i)];
    ax1 = subplot(2,steps,i+steps);
    c = categorical({'MCG','MPCG','SG','SPG','SCG','SPCG','Jacobi','GaussSeidel'});
    c = reordercats(c,{'MCG','MPCG','SG','SPG','SCG','SPCG','Jacobi','GaussSeidel'});
    b = bar(ax1,c,timeData);
    % ylim([minTime maxTime])
    ylim([0.0001 0.4])
    b.FaceColor = 'flat';
    b.CData(1,:) = [0,1,0];
    b.CData(2,:) = [0,1,0];
    b.CData(3,:) = [1,0,0];
    b.CData(4,:) = [1,0,0];
    b.CData(5,:) = [0,0,1];
    b.CData(6,:) = [0,0,1];
    b.CData(7,:) = [1,1,0];
    b.CData(8,:) = [1,0,1];
    titleTimeBar = sprintf('Time, DIM = %d', dimens(i));
    title(titleTimeBar)   
end