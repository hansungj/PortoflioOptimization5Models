clc
clear all
format long
gurobi_setup
rng(1)
% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the stock weekly prices and factors weekly returns
adjClose = readtable('Project2_Data_adjClose.csv', 'ReadRowNames', true);
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Properties.RowNames));

factorRet = readtable('Project2_Data_FF_factors.csv', 'ReadRowNames', true);
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

riskFree = factorRet(:,4);
factorRet = factorRet(:,1:3);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of shares outstanding per company
mktNoShares = [36.40 9.86 15.05 10.91 45.44 15.43 43.45 33.26 51.57 45.25 ...
                69.38 8.58 64.56 30.80 36.92 7.70 3.32 54.68 38.99 5.78]' * 1e8;

% Initial budget to invest
initialVal = 100;

% Start of in-sample calibration period 
calStart = datetime('2012-01-01');
calEnd   = calStart + calmonths(12) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2013-01-01');
testEnd   = testStart + calmonths(6) - days(1);

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Investment strategies
% Note: You must populate the functios MVO.m, MVO_card.m and BL.m with your
% own code to construct the optimal portfolios. 
funNames  = {'MVO' 'Robust_MVO' 'Resampling_MVO' 'Most_D_MVO' 'cVar' };
NoMethods = length(funNames);

funList = {'MVO' 'MVO_robust' 'MVO_Res', 'MVO_MD' 'cVar'};
funList = cellfun(@str2func, funList, 'UniformOutput', false);


%--------------------------------------------------------------------------

    %parameters for Robust_MVO
    lambda = 50;
    alpha = 0.1; %confidence interval of 1 - alpha = 0.9
    
    %parameters for Resampling MVO
    T = 100;
    reopt = 60; % number of times you need to re-optimize
    
    %paramerest for Index traking (Most-Diverse)
    z = 12;
    %rho = correlation coefficient matrix will be calculated in the for
    %loop

%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Construct and rebalance your portfolios
%
% Here you will estimate your input parameters (exp. returns, cov. matrix,
% etc) from the Fama-French factor models. You will have to re-estimate 
% your parameters at the start of each rebalance period, and then 
% re-optimize and rebalance your portfolios accordingly. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate counter for the number of observations per investment period
toDay = 0;
npaths_S = [100 500 1000 1500 2500 3000];
NoS = size(npaths_S',1);

Beta_S = [0.7 0.75 0.8 0.85 0.9];

for t = 1 : NoPeriods
    
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    periodReturns = table2array( returns( calStart <= dates & dates <= calEnd, :) );
    periodFactRet = table2array( factorRet( calStart <= dates & dates <= calEnd, :) );
    currentPrices = table2array( adjClose( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :) )';
    periodriskFree = table2array( riskFree( calStart <= dates & dates <= calEnd, :) );
    
    % Subset the prices corresponding to the current out-of-sample test 
    % period.
    periodPrices = table2array( adjClose( testStart <= dates & dates <= testEnd,:) );
    
    % Set the initial value of the portfolio or update the portfolio value
    
        for i = 1 : NoMethods + NoS + size(Beta_S',1)
            if t == 1
        
            currentVal(t,i) = initialVal;
        
            else
            currentVal(t,i) = currentPrices.' * NoShares{i};
            
            end
            %Ex Post Sharpe ratio              
            if t ~= 1
               SR_ep(t,i) = mean(periodReturns*x{i}(:,t-1) - periodriskFree)./var(periodReturns*x{i}(:,t-1));
            else
                SR_ep(t,i) = 0;
            end
 
        end
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    %----------------------------------------------------------------------
    
    % Calculate your initial exp. returns and covariance matrix using the
    % Fama-French factor model. You may find these values as you prefer,
    % the code given below is just an example. 
    
   
    %design matrix B
    B = [ ones(size(periodFactRet,1),1) periodFactRet ];
    
    %coefficient matrix calculation
    Coefficients = [];
    epsilon = [];
    
    for j = 1:size(periodReturns,2)
        
       Coefficients = [Coefficients inv(B.'*B)*B.'*periodReturns(:,j)];
       epsilon = [epsilon periodReturns(:,j) - B*Coefficients(:,j)];
       
    end
    
    A = Coefficients(1,:).';
    V = Coefficients(2:4,:);    % m x n matrix of betas
    
    f_bar = (geomean(periodFactRet + 1) - 1).'; % m x 1 vector of factor expected returns
    F = cov(periodFactRet);       % m x m factor covariance matrix
    
    epsilon = epsilon; % Regression residuals
    D = diag(var(epsilon));  % Diagonal n x n matrix of residual variance

    mu = A + V.'*f_bar; % n x 1 vector of asset exp. returns
    Q  = V.'*F*V + D;   % n x n asset covariance matrix
    
    %correlation coefficient matrix
    rho = [];
    for i = 1:size(Q)
        for j = 1:size(Q)
            rho(i,j) = Q(i,j)/(sqrt(Q(i,i))*sqrt(Q(j,j)));
            rho(j,i) = Q(i,j)/(sqrt(Q(i,i))*sqrt(Q(j,j)));
        end;
    end;
    
    %parameters for Monte-Carlo
    So = periodPrices(end,:);
    sigma = sqrt(diag(Q));
    Beta = 0.95;
    T = 26;
    npaths = 2000;
   
    %----------------------------------------------------------------------
    
    % Define the target return for the 2 MVO portfolios
    targetRet = mean(mu);
    
    % Optimize your portfolios to get the weights 'x'
    
    %MVO
    x{1}(:,t) = funList{1}(mu, Q, targetRet); 
    %MVO_Robust
    x{2}(:,t) = funList{2}(mu, Q, lambda, alpha); 
    %MVO_Resampling
    x{3}(:,t) = funList{3}(mu, Q, T, reopt, targetRet);
    %Most-divere MVO
    x{4}(:,t) = funList{4}(mu, Q, rho, targetRet, z);
    %CVaR Optimization with Monte Carlo simulation
    Scenario = CVaR_MC(So, mu, Q, npaths,T, sigma);
    x{5}(:,t) = funList{5}(Scenario,Beta, targetRet);
   
    NoS = size(npaths_S',1);
    
    for k = 1:NoS
        Scenario = CVaR_MC(So, mu, Q, npaths_S(k),T, sigma);
        x{5+k}(:,t) = funList{5}(Scenario, Beta, targetRet);
    end
    %Effect of Beta on Scenarios
    
    for k = 1:size(Beta_S',1)
        Scenario = CVaR_MC(So, mu, Q, npaths,T, sigma);
        x{5+ NoS + k}(:,t) = funList{5}(Scenario,Beta_S(k), targetRet);
    end
    
    %risk free rate for the period
    Rf = mean(periodriskFree); 
    % Calculate the optimal number of shares of each stock you should hold
    for i = 1:(NoMethods + NoS + size(Beta_S',1))
        

        % Number of shares your portfolio holds per stock
        NoShares{i} = x{i}(:,t) .* currentVal(t,i) ./ currentPrices;
        
        
        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay,i) = periodPrices * NoShares{i};

        %------------------------------------------------------------------
        
        % Calculate your transaction costs for the current rebalance
        % period. The first period does not have any cost since you are
        % constructing the portfolios for the 1st time.
        % transaction cost is 0.5% of the traded volume
        
        %%%%Portfolio returns using all scenarios%%%%
        
        monte_port = Scenario*x{i}(:,t);
        A = sort(monte_port);
        VaR(i,t) = 26*A(round(npaths*(1-Beta)),:);
        CVaR_port(i,t) = 26*mean(A(1:round(size(Scenario,1)*(1-Beta)),1));
        eReturns(t,i) = mu'*x{i}(:,t);
        portVar(t,i) = (x{i}(:,t))'*Q*(x{i}(:,t));
        ePeriod(t,i) = (1+eReturns(t,i))^(toDay - fromDay);
        
        %Ex Ante Sharpe Ratio Calculation
        SR_ea(t,i) = (mu.'*x{i}(:,t) - Rf)/sqrt(x{i}(:,t).'*Q*x{i}(:,t));
        
        if t ~= 1
            
           tCost(t-1, i) = 0.005 * (dot(abs(NoShares{i} - NoSharesOld{i}),currentPrices));
           
        end
        
        
        
        NoSharesOld{i} = NoShares{i};
        %------------------------------------------------------------------
        
    end
     
    % Update your calibration and out-of-sample test periods
    calStart = calStart + calmonths(6);
    calEnd   = calStart + calmonths(12) - days(1);
    
    testStart = testStart + calmonths(6);
    testEnd   = testStart + calmonths(6) - days(1);
    
    %----------------------------------------------------------------------
%%%%
    %----------------------------------------------------------------------

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------

% Calculate the portfolio average return, variance (or standard deviation),
% or any other performance and/or risk metric you wish to include in your
% report.

%FINAL CVAR
%Weekly Returns
weeklyReturn = zeros(156, NoMethods+NoS + size(Beta_S',1));
for t = 1:156
    for i = 1:NoMethods + NoS + size(Beta_S',1)
        weeklyReturn(t,i) =  (portfValue(t+1,i) - portfValue(t,i))/portfValue(t,i);
    end
end

%Statistics
weeklyReturnVar = std(weeklyReturn);
avgWeeklyReturns = geomean(weeklyReturn+1)-1;
avgEReturns = geomean(eReturns+1)-1;

for i = 1:NoMethods
    overallReturn(i) = (portfValue(toDay,i) - initialVal)/initialVal;
end

for i = 1:NoPeriods
    periodWeeklyReturns(i,:) = geomean(weeklyReturn(i:i+25,:)+1);
end

investHorizon = 3;
for i = 1:size(funNames,2)
    %average return per year
    avgReturns(1,i) = (portfValue(end,i) - 100)/investHorizon;
    %variance of returns
    varReturns(1,i) = var(portfValue(:,i));
end 

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 4.1 Plot the portfolio values 
%--------------------------------------------------------------------------
testStart = datetime('2013-01-01');
testEnd = testStart + 6*calmonths(6) - days(1);
plotDates = dates(testStart <= dates); 


fig1 = figure(1);
plot(plotDates, portfValue(:,1))
hold on
plot(plotDates, portfValue(:,2))
hold on
plot(plotDates, portfValue(:,3))
hold on
plot(plotDates, portfValue(:,4))
hold on
plot(plotDates, portfValue(:,5))
legend(funNames, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio Value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig1,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig1,'fileName','-dpng','-r0');

%--------------------------------------------------------------------------
% 4.2 Plot the portfolio weights 
%--------------------------------------------------------------------------

% MVO Plot
fig2 = figure(2);
area(x{1}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('MVO Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig2,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig2,'fileName2','-dpng','-r0');


% MVO_Robust
fig3 = figure(3);
area(x{2}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Robust MVO Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig3,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig3,'Position');
set(fig3,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig3,'fileName3','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig3,'fileName3','-dpng','-r0');


% Resampling MVO
fig4 = figure(4);
area(x{3}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Resampling MVO Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig4,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig4,'Position');
set(fig4,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);
print(fig4,'fileName4','-dpng','-r0');

% Most-diverse MVO
fig5 = figure(5);
area(x{4}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Most Diverse MVO Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig5,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig5,'Position');
set(fig5,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

print(fig5,'fileName5','-dpng','-r0');

%CVaR Optimization with Monte Carlo
fig6 = figure(6);
area(x{5}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('CVaR Portfolio Weights', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig6,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig6,'Position');
set(fig6,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

print(fig6,'fileName6','-dpng','-r0');

t = 1:NoPeriods;

fig7 = figure(7);
for i = 1:NoMethods
    plot(t,VaR(i,:))
    %plot(t, CVaR_port(i,:))
    hold on
end

legend(funNames, 'Location', 'eastoutside','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);
set(gca,'XTickLabelRotation',30);
title('Value at Risk (Beta = 0.95, Scenarios = 2000)', 'FontSize', 14)
ylabel('Value at Risk','interpreter','latex','FontSize',12);
print(fig7,'fileName7','-dpng','-r0');



fig8 = figure(8);
for i = 1:NoMethods
    plot(t,CVaR_port(i,:))
   
    hold on
end

legend(funNames, 'Location', 'eastoutside','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);
set(gca,'XTickLabelRotation',30);
title('Conditional Value at Risk (Beta = 0.95, Scenarios = 2000', 'FontSize', 14)
ylabel('CVaR','interpreter','latex','FontSize',12);
print(fig8,'fileName8','-dpng','-r0');
numpaths_s = {'2000 Scenarios' '100 Scenarios' '500 Scenarios' '1000 Scenarios' '1500 Scenarios' '2500 Scenarios' '3000 Scenarios'};

fig9 = figure(9);
for i = 1:NoS+1
    plot(plotDates, portfValue(:,4+i))
    hold on
end

legend(numpaths_s, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Effect of Scenarios on CVaR Portfolio Value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);
% Define the plot size in inches
set(fig9,'Units','Inches', 'Position', [0 0 8, 5]);
pos9 = get(fig9,'Position');
set(fig9,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos9(3), pos9(4)]);
print(fig9,'fileName9','-dpng','-r0');
Beta_s = {'0.95' '0.7' '0.75' '0.8' '0.85' '0.90'}
fig10 = figure(10);
plot(plotDates, portfValue(:,5))
hold on
for i = 1:size(Beta_S',1)
    plot(plotDates, portfValue(:,5+NoS+i))
    hold on
end

legend(Beta_s, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Effect of Beta on CVaR Portfolio Value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);
% Define the plot size in inches
set(fig10,'Units','Inches', 'Position', [0 0 8, 5]);
pos10 = get(fig10,'Position');
set(fig10,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos10(3), pos10(4)]);
print(fig10,'fileName10','-dpng','-r0');
P = [1;2;3;4;5]
fig11 = figure(11)

plot(P, tCost(:,1))
hold on
plot(P, tCost(:,2))
hold on
plot(P, tCost(:,3))
hold on
plot(P, tCost(:,4))
hold on
plot(P, tCost(:,5))
legend(funNames, 'Location', 'eastoutside','FontSize',12);
title('Transaction Costs', 'FontSize', 14)
xlabel('Rebalance Period','interpreter','latex','FontSize',12);
ylabel('Value','interpreter','latex','FontSize',12);
set(fig11,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig11,'Position');
set(fig11,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

print(fig11,'fileName11','-dpng','-r0');

a = [zeros(1,5); weeklyReturn(:,1:5)];
fig12 = figure(12);
plot(plotDates, a(:,1))
hold on
plot(plotDates, a(:,2))
hold on
plot(plotDates, a(:,3))
hold on
plot(plotDates, a(:,4))
hold on
plot(plotDates, a(:,5))
legend(funNames, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Weekly Portfolio Returns', 'FontSize', 14)
ylabel('Percentage','interpreter','latex','FontSize',12);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End