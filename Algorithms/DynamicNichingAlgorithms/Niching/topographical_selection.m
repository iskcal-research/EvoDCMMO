function species = topographical_selection(matdis, N, D)
    k = round(0.333*D+0.777*sqrt(N));
    
	[~, nn] = mink(matdis, k+1, 2);
    nn = nn(:, 2:end);
    
    outdegree = repmat([1:N]', 1, k) > nn;
    outdegree = sum(outdegree, 2);

    seeds = find(outdegree==0);
    
    seed_dis = matdis(:, seeds);
    [~, min_idx] = min(seed_dis, [], 2);

    species = struct();
    flag = false(1, length(seeds)); % sometimes two individuals are identical which makes the species empty
    for i=1:length(seeds)
        species(i).idx = sort(find(min_idx==i));
        species(i).seed = min(species(i).idx);
        species(i).len = length(species(i).idx);
        flag(i) = species(i).len == 0;
    end
    species(flag)=[];
end