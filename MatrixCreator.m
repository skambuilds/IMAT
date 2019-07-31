function [A] = MatrixCreator(mtype, n)

switch mtype
    case 1
        m = 2;
        prompt_dominante = sprintf('Il fattore dominante è al momento settato a %d, vuoi modificarlo? [S][N]', m);
        r = input(prompt_dominante, 's');
        if isempty(r)
            r = 'N';
        end
        if strcmp(r,'s') || strcmp(r,'S')
            prompt_dominante = 'Inserisci il fattore dominante: ';
            m = input(prompt_dominante);
        end
        A = rand(n);
        A = A - diag(diag(A));
        s = abs(A)*(ones(n,1));
        dA = A + m*diag(s);
        A = sparse(dA); % rende sparsa una matrice piena
                
    case 2
        rcond = 0.01; % reciproco del numero di condizionamento (in questo caso sarà circa 100)
        prompt_tau_rcond = 'Il reciproco del numero di condizionamento è fissati a 0.01, vuoi modificarlo? [S][N]';
        r = input(prompt_tau_rcond, 's');
        if isempty(r)
            r = 'N';
        end
        if strcmp(r,'s') || strcmp(r,'S')
            prompt_rcond = 'Inserisci il reciproco del numero di condizionamento (es. 0.01 per 100): ';
            rcond = input(prompt_rcond);
        end
        A = sprandsym(n, 1, rcond, 1); % random simmetrica (con l'ultimo parametro specifico che sia definita positiva)
            
    case 3
        tau = 0.05; % 5% di sparsità per la matrice
        rcond = 0.01; % reciproco del numero di condizionamento (in questo caso sarà circa 100)
        prompt_tau_rcond = 'La percentuale di sparsità e il reciproco del numero di condizionamento sono fissati a 0.05 e 0.01, vuoi modificarli? [S][N]';
        r = input(prompt_tau_rcond, 's');
        if isempty(r)
            r = 'N';
        end
        if strcmp(r,'s') || strcmp(r,'S')
            prompt_tau = 'Inserisci la percentuale di sparsità (0.1 per ottenere il 10% di elementi non nulli): ';
            tau = input(prompt_tau);
            prompt_rcond = 'Inserisci il reciproco del numero di condizionamento (es. 0.01 per 100): ';
            rcond = input(prompt_rcond);
        end
        A = sprandsym(n, tau, rcond, 1); % sparsa random simmetrica (con l'ultimo parametro specifico che sia definita positiva)
            
    case 4
        A = gallery('poisson', n); % avremo una matrice molto ristretta vicino alla diagonale (Block tridiagonal matrix from Poisson's equation (sparse))
    otherwise
        disp('Selezione non valida');
end

fprintf('Matrice tipologia [%d] con dimensione [%d] creata con successo\n', mtype, n);

