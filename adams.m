function f = adams(x,y,A)
% function variables
    A = readmatrix("Adjacency_matrix.csv")
    x = zeros(91, 1)
    y = zeros(91, 1)

% ADAM optimization constraints
    gnorm = zeros(length(A), 1); 
    beta1 = 0.9;
    beta2 = 0.999; 
    e = 10e-8;
    alpha = 0.001;
    m = x; % first moment vector
    v = y; % second moment vector

    for k = 1:kmax 
        indices = randperm(n,bsz);
        
        g = gfun(indices, w);
        f(k) = fun(indices,w);
        gnorm(k) = norm(g);

        m = beta1*m + (1-beta1).*g; % Update biased first moment estimate
        v = beta2*v + (1-beta2).*g.^2; % Update biased second moment estimate
        mt = m ./ (1-beta1^k); % bias corrected first moment estimate
        vt = v ./ (1-beta2^k); % bias corrected second raw moment estimate
        
        w = w -  alpha*mt ./ (sqrt(vt)  + e);

        
        % Check for convergence
        if gnorm(k) < tol
            break;
        end
    end
end