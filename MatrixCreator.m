function [A] = MatrixCreator(mtype, n)

switch mtype
    case 1
        m = 2;
        prompt_dominant = sprintf('The dominant factor is currently set to %d, do you want to change it? [Y][N]', m);
        r = input(prompt_dominant, 's');
        if isempty(r)
            r = 'N';
        end
        if strcmp(r,'s') || strcmp(r,'Y')
            prompt_dominante = 'Enter the dominant factor: ';
            m = input(prompt_dominant);
        end
        A = rand(n);
        A = A - diag(diag(A));
        s = abs(A)*(ones(n,1));
        dA = A + m*diag(s);
        A = sparse(dA); % makes a full matrix sparse
                
    case 2
        rcond = 0.01; % reciprocal of the condition number (in this case it will be about 100)
        prompt_tau_rcond = 'The reciprocal of the condition number is set at 0.01, do you want to change it? [Y][N]';
        r = input(prompt_tau_rcond, 's');
        if isempty(r)
            r = 'N';
        end
        if strcmp(r,'s') || strcmp(r,'Y')
            prompt_rcond = 'Enter the reciprocal of the condition number (es. 0.01 per 100): ';
            rcond = input(prompt_rcond);
        end
        A = sprandsym(n, 1, rcond, 1); % random simmetrica (con l'ultimo parametro specifico che sia definita positiva)
            
    case 3
        tau = 0.05; % 5% of sparsity for the matrix
        rcond = 0.01; % reciprocal of the condition number (in this case it will be about 100)
        prompt_tau_rcond = 'The percentage of sparsity and the reciprocal of the condition number are set at 0.05 and 0.01, do you want to change them? [Y][N]';
        r = input(prompt_tau_rcond, 's');
        if isempty(r)
            r = 'N';
        end
        if strcmp(r,'s') || strcmp(r,'Y')
            prompt_tau = 'Enter the sparse percentage (0.1 to get 10% of non-null elements): ';
            tau = input(prompt_tau);
            prompt_rcond = 'Enter the reciprocal of the condition number (eg 0.01 for 100): ';
            rcond = input(prompt_rcond);
        end
        A = sprandsym(n, tau, rcond, 1); % random symmetric sparse (with the last specific parameter being positive definite)
            
    case 4
        A = gallery('poisson', n); % Block tridiagonal matrix from Poisson's equation (sparse)
    otherwise
        disp('Invalid selection');
end

fprintf('Matrix of type [%d] with dimension [%d] successfully created\n', mtype, n);

