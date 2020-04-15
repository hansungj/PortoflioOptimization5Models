function  x_optimal = MVO_robust(mu, Q, lambda, alpha)
%MVO_Robust code

n = size(Q,1);

%Robust model
%needed parameters
epsilon = sqrt(chi2inv(1 - alpha,n));
theta = [diag(diag(Q))./n zeros(n,1); zeros(1,n) -1];

%Gurobi 
clear model
model.obj = [-mu.' epsilon];
model.Q = sparse(lambda*2*[Q zeros(n,1);  zeros(1,n+1)]);
model.modelsense = 'min';

%constraints
model.A = sparse([ones(n,1);0]');
model.rhs = 1;
model.sense = '=';
model.lb = zeros(n+1,1);

%Quadratic Constraints
model.quadcon(1).Qc = sparse(theta);
model.quadcon(1).q = zeros(n+1,1);
model.quadcon(1).rhs = 0;
model.quadcon(1).sense = '=' ;

%params
clear params;
params.TimeLimit = 200;
params.OutputFlag = 0;

robust = gurobi(model,params);
x_optimal = robust.x(1:n);

    
end