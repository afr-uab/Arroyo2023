function esr = gsea(score, hits)


% hits is rank of hits in ordered list

n = length(score);
nh = length(hits);

p = 1;

nr = sum(abs(score(hits)) .^ p);

phit = zeros(n, 1);
pmiss = zeros(n, 1);

for i = 1:n
    
    k = find(hits <= i);
    
    if ~isempty(k)
        phit(i) = sum((abs(score(hits(k))) .^ p) / nr);
    end
    
    pmiss(i) = (1 / (n - nh)) * (i - length(k));
    
end

esr = phit - pmiss;

return


