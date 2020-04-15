function  x_optimal = MVO_Res(mu, Q, T, reopt, targetRet)
rng(1); 
%MVO_Resampling Code
%allowed to short sell . Draw T = 100 randomly generated observations to
%re-estimate mu' and Q', re-optimize 60 times to construct your average
%optimal portfolio
   n = size(Q,1); 
   %budget constraint
   Aeq = ones(1,n); 
   
   x_star = [];
   
   %mu should be a row vector 
   %collect T observations with mean = mu and variance = Q
   %mvnrnd computes R(i,:) using MU(i,:) and SIGMA(:,:,i)
   for i = 1:reopt
   R = mvnrnd(mu,Q,T);
   
   %each row corresponds to a randomly generate generated values
   
   mu_star = mean(R,1);
   Q_star = cov(R);
   
   x_star(:,i) = quadprog(Q_star,[],-mu_star,-targetRet,Aeq,1);
   
   end
   
   x_optimal = mean(x_star,2);
  
    %----------------------------------------------------------------------
    
end