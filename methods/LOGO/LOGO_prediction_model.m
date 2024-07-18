function prediction = LOGO_prediction_model(data_test, data_train, label_train, Para)

if isfield(Para, 'soft' )
    soft = Para.soft;
else
    soft = 1;
    Para.soft = soft;
end

if isfield(Para, 'sigma' )
    sigma = Para.sigma;                 % kernel width/ the number of the samples in a kernel
else
    if soft ==1; sigma = mean(std(patterns'));end
    if soft ==0; sigma = 10; end
    Para.sigma = sigma;
end

if isfield(Para, 'distance' )
    distance = lower(Para.distance);
    if ~any(strcmp({'sqeuclidean', 'cityblock'},distance))
        error('Distance should be sqeuclidean or cityblock')
    end
else
    distance = 'sqeuclidean';
    Para.distance = distance;
end

Weight = Para.Weight;

[dim, N_patterns] = size(data_test);

Weight = reshape(Weight, [dim, 1]);


Num_sample = size(data_test, 2);
% Multi classes
Uc = unique(label_train);
N_class = length(Uc);
for n=1:length(Uc)
    temp = find(label_train==n);
    index{n} =temp;
    N(n) = length(temp);
end
prediction = zeros(1, Num_sample);

for n = 1:size(data_test,2)
    Dist = zeros(1, N_class);
    test = data_test(:,n);
    for nc = 1:N_class
        data_nc = data_train(:, index{nc});
        if strcmp(lower(Para.distance), 'cityblock')
            
            temp = abs(data_nc - test*ones(1, N(nc)));
        elseif strcmp(lower(Para.distance), 'euclidean')
            temp = (data_nc - test*ones(1, N(nc))).^2;
        else
            disp('Wrong distance parameter')
        end
        
        dist    = (Weight)'*temp;
        prob = exp(-dist/sigma);
        if sum(prob)~=0
            prob = prob/sum(prob);
        else
            [~,I] = sort(dist);
            prob=zeros(size(I));
            prob(I(1))=1;
        end
        
        
        Dist(nc) = sum(dist.*prob);
        
        
    end
    [m, id] = min(Dist);
    prediction(n) = Uc(id);
    
end


end

function kernel_coef = KernelFunction(distance_sorted, K, kernel)
if K <= length(distance_sorted)
    u = distance_sorted(1:K)/distance_sorted(K);
else u = distance_sorted(1:end)/distance_sorted(end);
end
switch lower(kernel)
    case {'epanechnikov'}
        kernel_coef = (1-u.^2);
    case {'cosine'}
        kernel_coef = cos(u*pi/2);
    case {'quartic'}
        kernel_coef = (1-u.^2).^2;
    case {'tricube'}
        kernel_coef = (1-u.^3).^3;
    case {'triweight'}
        kernel_coef = (1-u.^2).^3;
end

end
