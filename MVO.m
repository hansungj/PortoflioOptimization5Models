function  x_optimal = MVO(mu, Q, targetRet)
    
    % Find the total number of assets
    n = size(Q,1); 

    %budget constraint
    Aeq = ones(1,n); 
    % A = -mu
    % b = -targetRet
    
    x_optimal = quadprog(Q,[],-mu.',-targetRet,Aeq,1);
    
    
    
    
    
end