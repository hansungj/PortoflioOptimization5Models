function  x_optimal = MVO_MD(mu, Q, rho, targetRet, k)
    
    %Most Divere MVO : a combination of the index tracking problem and the
    %nominal MVO problem, solve the index tracking problem to find a subset
    %of 12 assets ( = k) then use MVO to construct out portfolio. Allowed to short
    %sell
    
    %in our case index = original 20 assets with mu and Q
    
    n = size(rho,1);
    
    % y = n x 1, z = n x n
    rho = reshape(rho,[1,n*n]);
   
    %portfolio size constraint for y
    A_size = [zeros(1,n*n) ones(1,n)];
    b_size = k;
    
    
    %each asset has exactly one representative
    A_rep = zeros(n,n*n+n);
    for i = 1:n
        A_rep(i,((i-1)*10+1):i*10) = ones(1,10);
    end 
    b_rep = ones(n,1);
    
    %zij <= yj constraint
    A_ineq = [eye(n*n) -repmat(eye(n),n,1)];
    b_ineq = zeros(n*n,1);
   
    
    %construct the model
    clear model
    
    %set objective: z
    model.obj = [rho zeros(1,n)];
    model.modelsense = 'max';
 
    %Add constraints
    model.A = sparse([A_size;A_rep;A_ineq]);
    model.rhs = [b_size; b_rep;b_ineq];
    model.sense = [repmat('=', n+1,1); repmat('<', n*n,1)];
    model.vtype = 'B'; %[repmat('B', n*n+n, 1)];
  
    results = gurobi(model);
    
    y = results.x(n*n+1:end,1);
   
    Q = Q.*(y*y.')
    
    Aeq = ones(1,n);
    x_optimal = quadprog(Q,[],(-mu.*y).',-targetRet,Aeq,1);
    
    


    %----------------------------------------------------------------------
    
end