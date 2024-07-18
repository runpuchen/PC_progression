function ES_curve = calculate_ES(distance, target, hit)
N = length(distance);
[distance,id]=sort(distance);
target = target(id);
ES_curve = zeros(N,1);
h_hit = 1/sum(target == hit);
h_miss = 1/sum(target ~= hit);
for i = 2:N
    if target(i) == hit
        ES_curve(i) = ES_curve(i-1)+h_hit;
    else ES_curve(i) = ES_curve(i-1)-h_miss;
    end


end


end