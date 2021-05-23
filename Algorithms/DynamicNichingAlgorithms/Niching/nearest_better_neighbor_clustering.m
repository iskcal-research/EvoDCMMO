function species = nearest_better_neighbor_clustering(decs)
    alpha = 5.0;
    N = size(decs, 1);
    slbest = zeros(N, 3);
    slbest(:, 1) = 1:N;
    
    matdis = pdist2(decs, decs);
    matdis(logical(eye(N))) = inf;
    [slbest(:, 3), slbest(:, 2)] = min(matdis, [], 2);
    
    meandis = alpha * mean(slbest(:, 3));
    
    matdis(matdis >= meandis) = inf;
    matdis(triu(true(N,N))) = inf;
    [slbest(:, 3), slbest(:, 2)] = min(matdis, [], 2);
    slbest(slbest(:, 3) == inf, 2) = slbest(slbest(:, 3) == inf, 1);
    
    sgbest = zeros(N, 2);
    sgbest(:, 1) = 1:N;
    
    for i=1:N
        j = slbest(i, 2);
        k = slbest(i, 1);
        while j ~= k
            k = j;
            j = slbest(k, 2);
        end
        sgbest(i, 2) = k;
    end
    
    seed_index = sort(unique(sgbest(:,2))); 
    seed_len = length(seed_index);
    species  = struct();
    for i = 1:seed_len
        species(i).seed = seed_index(i);
        species(i).idx= find(sgbest(:,2)==seed_index(i));
        species(i).len = length(species(i).idx);
    end
    
    seed_x = decs(cat(2, species.seed), :);
    seed_dis = pdist2(seed_x, seed_x);
    seed_dis(logical(eye(seed_len))) = inf;
    
    mark = zeros(seed_len, 2);
    mark(:, 1) = 1:seed_len;
    mark(:, 2) = mark(:, 1);
    
    for i = 1:seed_len
        [~,midx] = min(seed_dis(i,:));
        if all(species(i).seed > species(midx).idx) 
            mark(i,2) = midx;
        end
    end
    
    for i = 1 : seed_len 
        j = mark(i,2);
        k = mark(i,1);
        while j~= k 
             k = j;
             j = mark(k,2); 
        end
        mark(i,2) = k;
    end
    
    flag = zeros(1,seed_len);
    for i = 1:seed_len
        if mark(i,1)~=mark(i,2)
            flag(i) = 1;
            sgbest(species(i).idx,2) = species(mark(i,2)).seed;
            species(mark(i,2)).idx = sort([species(i).idx;species(mark(i,2)).idx]);
            species(mark(i,2)).len = length(species(mark(i,2)).idx);
        end   
    end
    species(flag == 1) = [];
end