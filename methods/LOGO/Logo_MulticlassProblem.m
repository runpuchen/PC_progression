function Weight = Logo_MulticlassProblem(patterns, targets, Para)

%Logo_MulticlassProblem: logo algorithm for feature selction for multiclassication problem
%Y. Sun, S. Todorovic, and S. Goodison,
%Local Learning Based Feature Selection for High Dimensional Data Analysis
%IEEE Trans. on Pattern Analysis and Machine Intelligence, vol. 32, no. 9, pp. 1610-1626, 2010.
%--------------------------------------------------------------------------
%INPUT:
%     patterns:  training data: [x1,x2,...xn] Each column is an observation
%      targets:  class label = {1,2,...,C}
%         Para:  parameters. More specifically,
%   Para.sigma:  kernel width
%  Para.lambda:  regulariztion parameter
%    Para.plot:  1: plot of the learning process; 0: do not plot
%OUTPUT:
%       Weight:  weight of features
%--------------------------------------------------------------------------
%by Yijun Sun @University at Buffalo
%update history: Feb. 10/March 20, 2007/OCT 10, 2010/July 18, 2013
%% ==========================================================================
sigma = Para.sigma;                 % kernel width
lambda = Para.lambda;               % regulariztion parameter
plotfigure = Para.plotfigure;       % whether the progress and feature weights are plotted

Uc = unique(targets);
if min(Uc)~=1 | max(Uc)~=length(Uc);
    error('Targets should run from 1 to C !')
end
[dim,N_patterns] = size(patterns);

for n=1:length(Uc)
    temp = find(targets==n);
    index{n} =temp;
    N(n) = length(temp);
end

Original_dim = dim;
Original_index = 1:dim;
History = [];
Weight =  1/sqrt(dim)*ones(dim,1); %initial guess
History(:,1) = Weight;
P.lambda = lambda;
Difference =1;t=0;theta =[];

while  Difference>0.01 & t<=10;
    t = t+1;
    NM = zeros(dim,N_patterns);
    NH = zeros(dim,N_patterns);
    V = (Weight(:).^2)';
    
    for i = 1:N_patterns
        index_SameClass = index{targets(i)};
        index_DiffClass = [];
        for c = 1:length(Uc)
            if c~=targets(i)
                temp = index{c};
                index_DiffClass = [index_DiffClass;temp(:)];
            end
        end
        
        Temp_SameClass            = abs(patterns(:,index_SameClass) - patterns(:,i)*ones(1,length(index_SameClass)));
        Temp_DiffClass            = abs(patterns(:,index_DiffClass) - patterns(:,i)*ones(1,length(index_DiffClass)));
        
        if t==1
            dist_SameClass    = sum(Temp_SameClass,1)/sqrt(dim);
            dist_DiffClass    = sum(Temp_DiffClass,1)/sqrt(dim);
        else
            dist_SameClass    = (V)*Temp_SameClass;
            dist_DiffClass    = (V)*Temp_DiffClass;
        end
        temp_index_SameClass = find(dist_SameClass==0);
        
        prob_SameClass = exp(-dist_SameClass/sigma);prob_SameClass(temp_index_SameClass(1)) = 0;
        if sum(prob_SameClass)~=0;prob_S = prob_SameClass/sum(prob_SameClass);else;[dum,I] = sort(dist_SameClass);prob_S(I(1))=1;end
        prob_DiffClass = exp(-dist_DiffClass/sigma);
        if sum(prob_DiffClass)~=0;prob_D = prob_DiffClass/sum(prob_DiffClass);else;[dum,I] = sort(dist_DiffClass);prob_D(I(1))=1;end
        
        NH(:,i) = Temp_SameClass*prob_S(:);
        NM(:,i) = Temp_DiffClass*prob_D(:);
    end
    
    Z = NM-NH;
    CostDiff = 1000; Cost(1) = 10000;
    j=1;
    
    while CostDiff>0.0001*Cost(j) & j <= 50
        j= j+1;
        a = (Weight.^2)'*Z; % Margin
        Result = 1./(1+exp(a));
        descent = lambda*Weight-(Z*Result(:)).*Weight;
        
        P.Weight = Weight;
        P.descent = descent;
        [alpha, Cost(j)] = fminbnd(@(p) logocost_multiclass_cost(p, P, Z), 0, 1);
        
        Weight = Weight-alpha*descent;
        CostDiff = abs(Cost(j)-Cost(j-1));
    end
    
    Weight = abs(Weight);
    Difference = norm(abs(Weight/max(Weight)-History(:,t)/max(History(:,t))));%max(abs(Weight/max(Weight)-History(:,t)/max(History(:,t))));
    theta(t) = Difference;
    History(:,t+1) = Weight;
    
    if t==1;index_zeros = find(Weight/max(Weight)<=10^(-5));end
    if t>=2;index_zeros = find(Weight/max(Weight)<=10^(-5));end
    patterns(index_zeros,:)=[];
    dim = size(patterns,1);
    Weight(index_zeros)=[];
    History(index_zeros,:)=[];
    Original_index(index_zeros)=[];
end
temp = zeros(1,Original_dim);
temp(Original_index) = Weight.^2;
Weight = temp;

%Monitoring the feature weights
if plotfigure ==1
    figure;
    semilogy(theta,'-o','LineWidth',1,'MarkerFaceColor','w','MarkerSize',10)
    title('Theta');
    xlabel('Number of Iterations');
    ylabel('Difference')
    grid on
    boldify1
    drawnow
end

return
%% ==================End of the code===================================





