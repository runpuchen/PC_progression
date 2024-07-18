function TestMe_Multiclass(Num_RandomFeature)

%% Conduct experiment on a toy data to test the code
%Input:                     
%      Num_RandomFeature: number of randomly distributed features. You can even set No_RandomFeature = 10^6, if your computer has a large RAM and you are patient.
%% ---------------------------------------------------------------
close all
disp(['>>> The number of irrelevant features is ' num2str(Num_RandomFeature) '...'])

%% three-cluster data
N = 100;
cluster_1 = [randn(1,N)+3;randn(1,N)-3];
cluster_2 = [randn(1,N)+3;randn(1,N)+3];
cluster_3 = randn(2,N)-3;
targets = [ones(1,N),ones(1,N)*2,ones(1,N)*3];     %class label

Data = [cluster_1,cluster_2,cluster_3];            %pool them together    
Original_dim = size(Data,1);
Data = [Data;randn(Num_RandomFeature,size(Data,2))]; %add gaussian noise
[MIN,I] = min(Data,[],2);
[MAX,I] = max(Data,[],2);
for n=1:size(Data,1)
    Data(n,:) = (Data(n,:)-MIN(n))/(MAX(n)-MIN(n)); % normalize data so that each feature is comparable
end

%plot the first two features
index_1 = find(targets==1);
index_2 = find(targets==2);
index_3 = find(targets==3);

figure(1);hold on
plot(Data(1,index_1),Data(2,index_1),'o', 'LineWidth',1,'MarkerSize',8)
plot(Data(1,index_2),Data(2,index_2),'r*','LineWidth',1,'MarkerSize',8)
plot(Data(1,index_3),Data(2,index_3),'gs','LineWidth',1,'MarkerSize',8)
title('The first two features')
axis square;axis tight
boldify1
drawnow

%parameters
Para.plotfigure = 1;     % 1: plot of the result of each iteration; 0: do not plot
Para.sigma= 1;           % kernel width; If the algorithm does not converge, use a larger kernel width.
Para.lambda = 1;         % regularization parameter
%I arbitarily set sigma= 2 and lambda = 1. The proposed algorithm is not sensitive to parameters. 
%The algorithm can used for classification. The parameters can be learning via cross-validation (see the paper).

s = cputime;
Weight_2 = Logo_MulticlassProblem(Data, targets, Para);
CPUTime = cputime-s;
disp(['>>> The total CPU time is ' num2str(CPUTime) ' seconds.'])

figure;
semilogx(Weight_2/max(Weight_2),'-o','LineWidth',1,'MarkerFaceColor','w','MarkerSize',10)
hold on
plot([Original_dim,Original_dim],[0,1],'r--', 'LineWidth',2);
xlabel('Features')
ylabel('Feature Scores')
title('Feature Weights');
axis tight
boldify1

return
%% =====================End of The Code=====================================



