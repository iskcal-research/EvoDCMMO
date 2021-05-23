function species = jarvis_partrick(matdis)
    N = size(matdis, 1);
    params = round(lhsdesign(1, 2) * 25);
    J = max(params);
    K = min(params);
    
    [~, nn] = mink(matdis, J+1, 2);
    nn = nn(:, 2:end);
    
    m = zeros(N, 2);
    m(:, 1) = 1:N;    
    
    while any(m(:,2) == 0)
        idx = find(m(:, 2) == 0);
        idx = idx(1);
        queue = idx;
        
        while ~isempty(queue)
            cur = queue(1);
            if m(cur, 2) == 0 %进行处理
                cur_nn = nn(cur, :);
                cur_nn_nn = nn(cur_nn, :);
                for i = 1:J
                    if m(cur_nn(i), 2) == 0 && ismember(cur, cur_nn_nn(i, :))
                        % cur_nn(i) does not cluster and cur is an element
                        % of cur_nn_nn(i, :)
                        temp = intersect(cur_nn, cur_nn_nn(i, :));
                        if length(temp) > K
                            %m(cur_nn(i), 2) = idx;
                            queue = [queue cur_nn(i)];
                        end
                    end
                end
            end
            m(cur, 2) = idx;
            queue = queue(2:end);
        end
    end
    
    % construct the result
    seeds = unique(m(:, 2));
    species = struct();
    for i = 1:length(seeds)
%         species(i).seed = seeds(i);
        species(i).idx = sort(m(m(:, 2) == seeds(i), 1));
        species(i).seed = min(species(i).idx);
        species(i).len = length(species(i).idx);
    end
    
%     if length(species) > 1
%         disp(1);
%     end
end