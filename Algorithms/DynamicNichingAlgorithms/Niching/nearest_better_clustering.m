function species = nearest_better_clustering(decs, phi)
% find the nearest better neighour and the distance of each individual
    matdis = pdist2(decs, decs);
    NP = size(matdis, 1);
    nbc = zeros(NP, 3);
    nbc(1, :) = [1 -1 0]; % the best individual do not have the nearest better neighbour
    for i = 2:NP
        nbc(i, 1) = i;
        [nbc(i, 3), nbc(i, 2)] = min(matdis(i, 1:i-1));
    end

    meandis=phi*mean(nbc(2:NP,3));
    nbc(nbc(:,3)>meandis,2)=-1;
    nbc(nbc(:,3)>meandis,3)=0;
    seeds=nbc(nbc(:,2)==-1,1);
    
    % set the root of subtree which contain the current individual
    m = zeros(NP, 2);
    m(1:NP, 1) = 1:NP;
    for i = 1:NP
       j = nbc(i, 2);
       k = j;
       while j ~= -1
           k =j;
           j = nbc(j, 2);
       end
       if k == -1
           m(i, 2) = i;
       else
           m(i, 2) = k;
       end
    end
    
    % construct the result
    species = struct();
    for i = 1:length(seeds)
        species(i).seed = seeds(i);
        species(i).idx = sort(m(m(:, 2) == seeds(i), 1));
        species(i).len = length(species(i).idx);
    end
end

