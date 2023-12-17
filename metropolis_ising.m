N = 30; % lattice size
beta = 0.2:0.01:1;
beta_plot = 0.2:0.0001:1
kmax = 1e8;
magvec = zeros(1, length(beta));
beta_vec = zeros(size(beta_plot));
Tc = 0.4408;

for m = 1:length(beta)
    spins = sign(ones(N,N)); % up spin as mentioned in question
    run_mean = 0;
    run_var = 0;

    for k = 1:kmax
        % Choosing a random point on the N by N lattice
        i = randi(N);
        j = randi(N);
        % delta H from part 1 of question 3
        deltaH = 2*spins(i,j)*(spins(mod(i-2,N)+1,j)+spins(mod(i,N)+1,j)+spins(i,mod(j-2,N)+1)+spins(i,mod(j,N)+1));
        
        if or(deltaH <= 0, rand() < exp(-beta(m) * deltaH))
            spins(i, j) = -spins(i, j);
        end

        % Running mean and variance
        mag = sum(spins,'all')/N^2;
        run_mean = (k*run_mean+mag)/(k+1);
        run_var = ((k-1)*run_var+(mag-run_mean)^2)/k;

    end
    magvec(m) = run_mean; % input run_mean into vector to plot graph
end

% for loop for the analytical function
for p = 1:length(beta_plot)
    if beta_plot(p) > Tc
        beta_vec(p) = (1-(sinh(2*beta_plot(p)))^(-4))^(1/8);
    else 
        beta_vec(p) = 0;
    end
end

%plot graph
figure;
plot(beta_plot, beta_vec, '-', 'LineWidth',2);
hold on;
plot(beta, magvec, 'o');
hold on;
plot(beta, magvec + sqrt(run_var), '-.', 'LineWidth',2);
hold on;
plot(beta, magvec - sqrt(run_var), '-.', 'LineWidth',2)
xlabel('Beta');
ylabel('Mean Magnetization');
title('2D Ising Model Plot')

