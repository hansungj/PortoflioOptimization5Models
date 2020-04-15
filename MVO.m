function  x_optimal = MVO(mu, Q, targetRet)
    
    % Use this function to construct your MVO portfolio subject to the
    % target return, with short-selling allowed. 
    %
    % You may use quadprog, Gurobi, or any other optimizer you are familiar
    % with. Just be sure to comment on your code to (briefly) explain your
    % procedure.

    % Find the total number of assets
    n = size(Q,1); 
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    %budget constraint
    Aeq = ones(1,n); 
    % A = -mu
    % b = -targetRet
    
    x_optimal = quadprog(Q,[],-mu.',-targetRet,Aeq,1);
    
    
    
    
    
    
    
    

    %----------------------------------------------------------------------
    
end