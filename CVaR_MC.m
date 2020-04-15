function  scenarios = CVaR_MC(So, mu, Q, npaths,T, sigma)
rng(1);
    sigma = sqrt(diag(Q));
    steps = 1;
    T = 26;
    dt = T / steps;
    n = size(Q,1);
    for i = 1:n   
        for j = 1:n
            rho(i,j) = Q(i,j)/(sigma(i)*sigma(j));
        end
    end
    
    L = chol(rho, 'lower');
    scenarios = zeros(npaths, n);

    for i = 1:npaths
        for j = 1:n
            xi = L*randn(n,1);
            S = So(j)*exp((mu(j) - 0.5 * sigma(j)^2 ) * dt ...
               + sigma(j) * sqrt(dt) * xi(j));
            scenarios(i,j) = (S - So(j))/(26*So(j));
        end
    end
end 