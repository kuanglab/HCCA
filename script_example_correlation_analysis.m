clear
% load ResHCCAPPI % load results

% assuming that we are evaluating correlation betweeb gene expression and
% climate, and the HCCA order is (Xmut,(Xmet,(Xexp,Cli)))

Ui1 = Uis{1};
Uj1 = Ujs{1};
Uj2 = Ujs{2};
Uj3 = Ujs{3};

cli = 138; % index of climate feature
n_genes = 500;

genes = cell(500,4);

% baseline
[~,idx] = sort(normc(Climate(:,cli))'*normc(Xexp),'descend');
genes(:,1) = Gexp(idx(1:n_genes));

% level 2
cliCorr = normc(Climate)'*normc(Uj1(:,1:2));
expCorr = normc(Xexp)'*normc(Ui1(:,1:2));
dist = pdist2(cliCorr,expCorr);
[~,idx] = sort(dist(cli,:));
genes(:,2) = Gexp(idx(1:n_genes));

% level 3
cliCorr = normc(Climate)'*normc(Uj2(:,1:2));
expCorr = normc(Xexp)'*normc(Uj2(:,1:2));
dist = pdist2(cliCorr,expCorr);
[~,idx] = sort(dist(cli,:));
genes(:,3) = Gexp(idx(1:n_genes));

% level 4
cliCorr = normc(Climate)'*normc(Uj3(:,1:2));
expCorr = normc(Xexp)'*normc(Uj3(:,1:2));
dist = pdist2(cliCorr,expCorr);
[~,idx] = sort(dist(cli,:));
genes(:,4) = Gexp(idx(1:n_genes));
