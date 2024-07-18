function p = estimate_PValue(dist, target, hit, max_value)
N_permute = 10000;
N_larger = 0;
N_samples = length(dist);
for i = 1:N_permute
    id = randsample(1:N_samples, N_samples);
    target_permute = target(id);
    temp = calculate_ES(dist, target_permute, hit);
    m = max(temp);
    if m>max_value
        N_larger = N_larger+1;
    end



end

p = N_larger/N_permute;

end