% Metodo del gradiente coniugato
% x - soluzione calcolata
% k - numero di iterazioni eseguite
% resvec - vettore che contiene il residuo ad ogni iterazione
function [x,k,resvec] = SelfConiugGradient(A,b,tau,maxn,x)
% Determino la dimensione della matrice
n=size(A,1);
% Inizializzo il vettore x0
x0=100*ones(n,1);
% Inizializzazione residuo
r = b - A*x;
% Inizializzazione direzione iniziale
p = r;
% Inizializzazione contatore iterazioni
k = 0;
% Avvio ciclo dell'algoritmo con condizione di controllo, posso scegliere
% tra:
% - controllo del residuo: norm(r)>tau*norm(b)
% - condizione di Cauchy: norm(x-x0)>tau*norm(x)
% Imposto comunque un numero di iterazioni massimo
% Preallocazione risorse per il vettore residuo
resvec=zeros(maxn,1);
while(norm(x-x0)>tau*norm(x)) && (k<maxn)
    x0=x;
    k=k+1;
    % Memorizzo la norma del residuo nel vettore
    resvec(k)=norm(r);
    % Ottimizzo calcolando solo una volta il prodotto matrice vettore:
    s = A*p;
    delta = p'*s;    
    % Calcolo del passo:
    alpha = (p'*r)/delta;
    % Calcolo nuovo soluzione
    x=x0+alpha*p;
    % Aggiornamento residuo
    r=r-alpha*s;
    beta = s'*r/delta;
    p = r-beta*p;    
end

