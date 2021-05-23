classdef DynamicNichingBSO < DynamicNichingBase
    methods
        function obj = DynamicNichingBSO(nt, rt, ct)
            obj = obj@DynamicNichingBase(nt, rt, ct);
        end
        
        function GetNextPop(obj)
            niching = obj.GetNiching();
            
            offspring = obj.pop;
            for i = 1:length(niching)
                parent = obj.pop(niching(i).idx);
                off_decs = obj.GetOffspring(parent);
                
                off = obj.hub.GetIndis(off_decs);
                if ~isempty(off)
                    offspring(niching(i).idx(1:length(off))) = off;
                end
            end
            
            cmp = obj.CompareWithConstraint(obj.pop, offspring);
            obj.pop(cmp) = offspring(cmp);
            
            [indis, index] = obj.MaintainDiversity(obj.pop);
            if ~isempty(index)
%                 disp(obj.hub.evaluated);
                obj.pop(index) = indis;
            end
        end
        
        function off_decs = GetOffspring(obj, parent)
            cluster_num = min(5, length(parent));
            idx = kmeans(parent.decs, cluster_num);
            cluster = obj.TransformToClusters(idx, parent);
            
            if rand(obj.algRand) < 0.2
                cluster_idx = randi(obj.algRand, length(cluster));
                rand_ind = obj.hub.GetIndis(rand(obj.algRand, 1, obj.hub.D) .* (obj.hub.upper-obj.hub.lower) + obj.hub.lower);
                if ~isempty(rand_ind)
                    parent(cluster(cluster_idx).seed) = rand_ind;
                end
            end
            
            selected_decs = obj.GetSelectedIndividuals(parent.decs, cluster);
            
            xi = logsig((0.5*obj.maxG-obj.curG)/20) .* rand(obj.algRand, length(parent), obj.hub.D);
            off_decs = selected_decs + xi .* randn(obj.algRand, length(parent), obj.hub.D);
            off_decs = boundary_check(off_decs, obj.hub.lower, obj.hub.upper);
        end
        
        function cluster = TransformToClusters(obj, idx, parent)
            num = length(unique(idx));
            cluster = struct();
            for i = 1:num
                cluster(i).idx = find(idx==i);
                cluster(i).len = length(cluster(i).idx);
%                 sort_idx = obj.SortWithContraint(parent(cluster(i).idx));
%                 cluster(i).seed = cluster(i).idx(sort_idx(1));
                cluster(i).seed = cluster(i).idx(1);
            end
        end
        
        % obtain the selected individuals for each individuals
        function decs = GetSelectedIndividuals(obj, parent, cluster)
            [n, d] = size(parent);
            decs = zeros(n, d);
            
            for i = 1:length(cluster)
                pro_one = rand(obj.algRand, cluster(i).len, 1) < 0.8;
                num_one = sum(pro_one);
                
                % construct the individuals from one cluster
                rand_idx = rand(obj.algRand, num_one, 1) * n;
                selected_cluster = ones(num_one, 1);
                cluster_num = 0;
                for j = 1:length(cluster)
                    cluster_num = cluster_num + cluster(j).len;
                    selected_cluster(rand_idx >= cluster_num) = selected_cluster(rand_idx >= cluster_num) + 1;
                end
                indi_one = zeros(num_one, d);
                for j = 1:num_one
                    if rand(obj.algRand) < 0.4
                        indi_one(j, :) = parent(cluster(selected_cluster(j)).seed, :);
                    else
                        selected_idx = randi(obj.algRand, cluster(selected_cluster(j)).len);
                        indi_one(j, :) = parent(cluster(selected_cluster(j)).idx(selected_idx), :);
                    end
                end
                
                % construct the individuals from two species
                indi_two = zeros(cluster(i).len - num_one, d);
                for j = 1:(cluster(i).len - num_one)
                    if length(cluster) > 1
                        cluster_idx = randperm(obj.algRand, length(cluster), 2);
                    else
                        cluster_idx = [1 1];
                    end
                    if rand(obj.algRand) < 0.5
                        selected_idx = [cluster(cluster_idx).seed];
                    else
                        r1 = cluster(cluster_idx(1)).idx(randi(obj.algRand, cluster(cluster_idx(1)).len));
                        r2 = cluster(cluster_idx(2)).idx(randi(obj.algRand, cluster(cluster_idx(2)).len));
                        selected_idx = [r1 r2];
                    end
                    par = rand(obj.algRand);
                    indi_two(j, :) = par .* parent(selected_idx(1), :) + (1-par) .* parent(selected_idx(2), :);
                end
                
                sub_decs = zeros(cluster(i).len, d);
                sub_decs(pro_one, :) = indi_one;
                sub_decs(~pro_one, :) = indi_two;
                decs(cluster(i).idx, :) = sub_decs;
            end
        end
    end
end