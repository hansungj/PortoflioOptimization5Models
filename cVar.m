function  x = cVar(ScenRets, beta, targetRet)

   [npaths n] = size(ScenRets);
  
   A = [-ScenRets -eye(npaths) -ones(npaths,1)];
   A = [-ones(1,n) zeros(1,npaths + 1); A];
  
   b = [-targetRet; zeros(npaths,1)];
   
   %f = [x z gamma]
   f = [zeros(1, n) (1/((1-beta)* npaths))*ones(1, npaths) 1];
   
   Aeq = [ones(1,n) zeros(1,npaths) 0];
   beq = 1;
   
   lb = [zeros(n,1); zeros(npaths,1); -Inf];
   ub = [ones(n,1); Inf*ones(npaths,1); Inf];
   
   x_optimal = linprog(f, A, b, Aeq, beq, lb, ub);
   x = x_optimal(1:n);
end
