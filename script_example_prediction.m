%% Load data
clear
load('Data/cvs');
load('Data/data');

%% Prepare results
Result = [];
Result.MSE = [];
Result.Correlation = [];
Result.RSquare = [];
Result.Yhat = {};
% 
%% 
%%% Run HCCA
% [Us,~] = HCCA({Xexp;Xmut;Xmet;Climate}, 4, 0.85);
% X = Us{end};

%%% Run HCCA with PPI network
% T1 = load('Data/STRING/PPI_Exp');
% T2 = load('Data/STRING/PPI_Mut');
% T3 = load('Data/STRING/PPI_Met');
% NW_F_Exp = T1.PPI;
% NW_F_Mut = T2.PPI;
% NW_F_Met = T3.PPI;
% LexpF = fct_graphLaplacian(NW_F_Exp);
% LmutF = fct_graphLaplacian(NW_F_Mut);
% LmetF = fct_graphLaplacian(NW_F_Met);
% [Us,~] = HCCA_PPI({Xexp;Xmut;Xmet;Climate}, 4, 0.85, {LexpF;LmutF;LmetF;[]});
% X = Us{end};

%%%%%%%%% Baselines:

%%% Run basic CCA with Xexp and Climate
[U1,V1,~,~,~] = MyCCA(Xexp, Climate, 4, 0.85);
X = [U1 V1];

%%% Run basic pairwise CCA with Xexp, Xmut and Climate
% [U1,V1,~,~,~] = MyCCA3(Xexp, Xmut, Climate, 4, 0.85);
% X = [U1 V1];

%%% Run basic pairwise CCA with Xexp, Xmut, Xmet and Climate
% [U1,V1,~,~,~] = MyCCA4(Xexp, Xmut, Xmet, Climate, 4, 0.85);
% X = [U1 V1];

%%% Run SNF for Xexp and Climate
% Xexp = zscore(Xexp);
% Climate = zscore(Climate);
% Dist1 = dist2(Xexp,Xexp);
% Dist2 = dist2(Climate,Climate);
% W1 = affinityMatrix(Dist1, 10, 0.7);
% W2 = affinityMatrix(Dist2, 10, 0.7);
% W = SNF({W1,W2}, 10, 20);
% [U,S,V] = svd(W);
% X = U*sqrtm(S);
% X = X(:,1:40);

%%% Run Stacked Xexp and Climate
% X = [Xexp Climate];

%% Run SVR
for j=1:200 % 200 partitions
    j
    cv = cvs{j};

    yhatall = [];

    for i=1:10 % 10-fold CV
        rng(1);
        Xtrn = X(cv.training(i),:);
        ytrn = y(cv.training(i),:);
        model = fitrsvm(Xtrn,ytrn,'KernelFunction', 'gauss','KernelScale','auto','Standardize',true);

        Xtst = X(cv.test(i),:);
        ytst = y(cv.test(i));
        yhat = predict(model,Xtst);
        yhatall = [yhatall; yhat];

        Result.MSE(j,i) = immse(ytst, yhat);
        corr = corrcoef(ytst, yhat);
        Result.Correlation(j,i) = corr(1,2);
        Result.RSquare(j,i) = MyRSquare(ytst, yhat, mean(ytrn));
    end
    Result.Yhat{j} = yhatall;
end

%% Report results
disp(num2str(mean(mean(Result.RSquare,2))));
