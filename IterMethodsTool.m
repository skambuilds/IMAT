fprintf('---Strumento di analisi metodi iterativi---\n');
fprintf('Premi Invio per iniziare\n');
pause
clear
% Elenco descrizione matrici
mat1 = '[1] - Matrice diagonalmente dominante sparsa';
mat2 = '[2] - Matrice simmetrica definita positiva';
mat3 = '[3] - Matrice simmetrica definita positiva e sparsa';
mat4 = '[4] - Matrice di Poisson';

% Scelta del tipo di matrice
matrix_str = sprintf('%s\n%s\n%s\n%s\n%s\n',mat1,mat2,mat3,mat4);
fprintf(matrix_str);
choice_str = 'Quale matrice desideri utilizzare? ';
mtype = input(choice_str);

% Scelta degli step di dimensione
steps_str = input('Inserisci gli step di dimensione che desideri testare separati da spazi o virgole: ', 's');
dimens = str2num(steps_str);
steps = size(dimens,2);

% Parametri per l'esecuzione dell'analisi
tol = 1e-12;
nmax = 200;
omega = 0.5;

fprintf('Di seguito le impostazioni per l''esecuzione del test:\n')
fprintf('%d - Tolleranza\n', tol);
fprintf('%d   - Numero massimo di iterazioni\n', nmax);
fprintf('%d   - Parametro di rilassamento\n', omega);
prompt_3 = sprintf('Vuoi modificarli? [S][N]');

r = input(prompt_3, 's');
if isempty(r)
    r = 'N';
end
if strcmp(r,'s') || strcmp(r,'S')
    prompt_n = 'Inserisci il valore di tolleranza: ';
    tol = input(prompt_n);
    prompt_n = 'Inserisci il numero massimo di iterazioni: ';
    nmax = input(prompt_n);
    prompt_n = 'Inserisci il valore omega per le operazioni di rilassamento: ';
    omega = input(prompt_n);
end


% Utilizziamo un pre-condizionatore - cholesky incompleta
% Paramentri necessari:
opts.type = 'ict'; % Seleziono la Cholesky incompleta con sogliatura - Cholesky with threshold dropping (ICT)
opts.droptol = 1e-3; % Sotto quale soglia i valori devono essre posti a zero (i numeri piccoli diventano zero)
% se aumento questo valore il precondizionatore diviene più sparso
opts.shape = 'upper'; % scelgo il formato upper (RT)(R)
% opts è di fatto una struct che contiene i parametri

% Inizializzazione delle variabili
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
    
    fprintf('Iterazione [%d] con dimensionamento [%d]\n', i, dimens(i));
    A = MatrixCreator(mtype,dimens(i));
    figure(i);
    spy(A); % In blu gli elementi non nulli della matrice selezionata
    
    % Costruzione del sistema lineare da risolvere
    n = size(A,1);
    % Impostazione del sistema, con soluzione pari ad un vettore di 1:
    sol = ones(n,1);
    b = A*sol;
    % In alternativa scegliamo il termine noto con un vettore casuale di dimensione n:
    % b = rand(n,1);
    x0= zeros(n,1);
    
    % Matlab gradiente coniugato
    tic
    [x, flag, relres, iter, resvec] = pcg(A,b,tol,nmax,[],[],x0); %potremo anche non inserire tol e nmax perchè vengono impostati di default
    time1(i) = toc;
    error1(i) = norm(sol-x)/norm(sol); % calcolo errore relativo
    x1(i) = {x};
    iter1(i) = iter;
    resvec1(i) = {resvec};
    
    % Valutare se eseguire la LU incompleta a seconda delle matrici scelte
%     setup.type = 'crout';
%     setup.milu = 'row';
%     setup.droptol = 1e-3;
%     setup.udiag = 1;
%     [L,U] = ilu(A,setup);
    
    % Matlab gradiente coniugato precondizionato
    R = ichol(A, opts);
    tic
    [x, flag, relres, iter, resvec] = pcg(A,b,tol,nmax,R',R,x0);
    time1p(i) = toc;
    error1p(i) = norm(sol-x)/norm(sol); % calcolo errore relativo
    x1p(i) = {x};
    iter1p(i) = iter;
    resvec1p(i) = {resvec};
    
    % Gradiente classico o di massima discesa
    tic
    [x, iter, resvec] = SelfGradient(A,b,tol,nmax,x0);
    time2(i) = toc;
    error2(i) = norm(sol-x)/norm(sol); % calcolo errore relativo
    x2(i) = {x};
    iter2(i) = iter;
    resvec2(i) = {resvec};
    
    % Gradiente classico precondizionato
    tic
    [x, iter, resvec] = SelfPreGradient(A,b,tol,nmax,R',R,x0);
    time2p(i) = toc;
    error2p(i) = norm(sol-x)/norm(sol); % calcolo errore relativo
    x2p(i) = {x};
    iter2p(i) = iter;
    resvec2p(i) = {resvec};
    
    % Gradiente coniugato
    tic
    [x, iter, resvec] = SelfConiugGradient(A,b,tol,nmax,x0); %potremo anche non inserire tol e nmax perchè vengono impostati di default
    time3(i) = toc;
    error3(i) = norm(sol-x)/norm(sol); % calcolo errore relativo
    x3(i) = {x};
    iter3(i) = iter;
    resvec3(i) = {resvec};
    
    % Gradiente coniugato precondizionato
    tic
    [x, iter, resvec] = SelfPreConiugGradient(A,b,tol,nmax,R',R,x0);
    time3p(i) = toc;
    error3p(i) = norm(sol-x)/norm(sol); % calcolo errore relativo
    x3p(i) = {x};
    iter3p(i) = iter;
    resvec3p(i) = {resvec};
    
    % Metodo di Jacobi
    tic
    [x, iter, resvec] = Jacobi(A,b,tol,nmax,x0);
    time4(i) = toc;
    error4(i) = norm(sol-x)/norm(sol); % calcolo errore relativo
    x4(i) = {x};
    iter4(i) = iter;
    resvec4(i) = {resvec};
    
    % Metodo di Gauss Seidel
    tic
    [x, iter, resvec] = GaussSeidel(A,b,tol,nmax,x0);
    time5(i) = toc;
    error5(i) = norm(sol-x)/norm(sol); % calcolo errore relativo
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

% Tabella tempi
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

% Tabella iterazioni
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

% Tabella errori
errorTable = struct2table(E)


% Grafici residui
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
titleGradPlot = sprintf('Residui Metodi Iterativi, DIM = %d', dimens(i));
title(titleGradPlot)
end

% Grafici a barre per iterazioni e tempi
figure(steps+2)
for i=1:steps
    iterData = [iter1(i) iter1p(i) iter2(i) iter2p(i) iter3(i) iter3p(i) iter4(i) iter5(i)];
    ax1 = subplot(2,steps,i);
    c = categorical({'MCG','MPCG','SG','SPG','SCG','SPCG','Jacobi','GaussSeidel'});
    c = reordercats(c,{'MCG','MPCG','SG','SPG','SCG','SPCG','Jacobi','GaussSeidel'});
    b=bar(ax1,c,iterData);
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
    titleIterBar = sprintf('Iterazioni, DIM = %d', dimens(i));
    title(titleIterBar)
    
    timeData = [time1(i) time1p(i) time2(i) time2p(i) time3(i) time3p(i) time4(i) time5(i)];
    ax1 = subplot(2,steps,i+steps);
    c = categorical({'MCG','MPCG','SG','SPG','SCG','SPCG','Jacobi','GaussSeidel'});
    c = reordercats(c,{'MCG','MPCG','SG','SPG','SCG','SPCG','Jacobi','GaussSeidel'});
    b=bar(ax1,c,timeData);
    %     ylim([minTime maxTime])
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
    titleTimeBar = sprintf('Tempi, DIM = %d', dimens(i));
    title(titleTimeBar)   
   
end
