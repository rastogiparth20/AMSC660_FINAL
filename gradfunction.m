function gradfunction()
    A = readmatrix("Adjacency_matrix.csv");

    kmax = 1e5; % Iterations
    N = 91; % Dimention
    tol = 1e-6; % tolerance value
    gnorm = zeros(N, 1);
    beta1 = 0.9;
    beta2 = 0.999; 
    e = 1e-8; % Adjust the small constant
    alpha = 0.001;
    w = randn(2*N, 1); % randomly generated vector
    m = zeros(2*N, 1); % Initialize first moment vector (Hint - column vector with 2*N components) from force.m code file
    v = zeros(2*N, 1); % Initialize second moment vector (Hint - column vector with 2*N components) from force.m code file
    residuals = 1000;

    for k = 1:kmax
        
        g = -forces(w(1:N),w(N+1:end),A); % force(x,y,A) 
        gnorm(k) = norm(g);

        m = beta1 * m + (1 - beta1) * g; % Update biased first moment estimate
        v = beta2 * v + (1 - beta2) * g.^2; % Update biased second moment estimate
        mt = m / (1 - beta1^k); % Bias-corrected first moment estimate
        vt = v / (1 - beta2^k); % Bias-corrected second raw moment estimate
        w = w - alpha * mt ./ (sqrt(vt) + e);

        if rem(k, residuals) == 0
            fprintf('Iteration %d, Norm of Residual: %e\n', k, gnorm(k));
        end

        if gnorm(k) < tol % Check for convergence
            break;
        end
    end
    
    plot_graph(w(1:N), w(N+1:end),A); % plot graph

    figure;
    plot(1:k, gnorm(1:k), 'o');
    set(gca, 'YScale', 'log');
    title("Norm of Force vs Iteration Number");
    xlabel("Iteration Number");
    ylabel("Log Norm of Force");
    grid on;
end