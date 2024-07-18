function experiment_cancer

% cancer datasets

load('Data_PCAProcessed','NUM','PAM50_extend','CENSOR','TIME');
X = NUM;
PAM50 = PAM50_extend;
clear NUM PAM50_extend

% centeralize data 
% X = compute_kernel(X);
X = X - repmat( mean(X,2), [1 size(X,2)] );

% get 95% dimensionality
fea = X';
[coeff,score,latent] = pca(fea);
cum = cumsum(latent);
petag = cum ./ cum(end);
idx = find(petag > 0.9);
rdim = idx(1);

% parameter settings
params.maxIter = 20;
params.eps = 1e-3;

params.dim = 80;%rdim;
params.lambda = 5 * size(X,2); 
params.sigma = 0.001;
params.gamma = 10;


% run DDRTree algorithm
[W, Z, stree, Y, history] = DDRTree(X, params);
plot_figure(Z, Y, stree, PAM50);

% % run DRTree algorithm
% [W, Z, stree, Y, history] = DDRTree_1(X, params);
% plot_figure(Y, Y, stree, PAM50);

function plot_figure(X, mu, stree, PAM50)

% draw cancer data
line_config ={'>','s','o','d','^','p'};
% names = {'1','2','3','4','5','6'};
names = {'Basal','Her2+','luminal A', 'luminal B', 'normal-like','normal'};
mpdc6 = distinguishable_colors(6);
numL = max(PAM50);

figure;
for i=1:numL
    idx = find(PAM50==i);
    hp(i)=plot3(X(1,idx), X(2,idx), X(3,idx),line_config{i},...
        'MarkerSize',5, 'MarkerEdgeColor',mpdc6(i,:),'MarkerFaceColor',mpdc6(i,:));
    hold on;
end
grid on;

for n = 1:size(stree)
    index = find(stree(n,:)>0);
    for m = 1:length(index)
        hline=plot3([mu(1,n), mu(1,index(m))],[mu(2,n), mu(2,index(m))],...
            [mu(3,n), mu(3,index(m))],'-k','LineWidth',3);
    end
end

set(gca, 'FontSize',16);
legend([hp hline],[ names(1:numL), 'tree structure']);