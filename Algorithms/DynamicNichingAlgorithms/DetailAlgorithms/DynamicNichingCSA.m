classdef DynamicNichingCSA < DynamicNichingBase
    properties(SetAccess=protected)
        copy = 2;
%         p_max = 10;
        min_delete = 10;
        % nc = 0.1;
    end
    
    methods(Access=public)   
        function obj = DynamicNichingCSA(nt, rt, ct)
            obj = obj@DynamicNichingBase(nt, rt, ct);
            obj.maxG = round(obj.hub.freq / (obj.hub.N*obj.copy));
        end
        
        function GetNextPop(obj)
%             disp(obj.hub.evaluated);
            niching = obj.GetNiching();
            
            [niching, indis, dlist] = obj.GetRedundancy(niching);
            dlist = dlist(1:length(indis));
            obj.pop(dlist) = indis;
            
            for i = 1:length(niching)
                sub_index = niching(i).idx;
                next = obj.GetNextSubPop(sub_index);
                if ~isempty(next)
                    obj.pop(sub_index(1:length(next))) = next;
                end
            end
            
            [indis, index] = obj.MaintainDiversity(obj.pop);
            if ~isempty(index)
%                 disp(obj.hub.evaluated);
                obj.pop(index) = indis;
            end
        end
        
        function [niching, indis, dlist] = GetRedundancy(obj, niching)
            cur_d = 0;  % current delete number
            dlist = []; % the index of delete
            % the size of the species can not exceed 'p_max'
%             disp(obj.hub.evaluated);
%             for i = 1:length(niching)
%                if niching(i).len > obj.p_max 
%                    cur_d = cur_d + niching(i).len - obj.p_max;
%                    dlist = [dlist; niching(i).idx(obj.p_max+1:end)];
%                    niching(i).idx = niching(i).idx(1:obj.p_max);
%                    niching(i).len = obj.p_max;
%                end
%             end
%             disp(obj.hub.evaluated);

            % ensure minimum number of deletions 
            i = length(niching);
            while cur_d < obj.min_delete
                cur_d = cur_d + 1;
                dlist = [dlist; niching(i).idx(end)];
                niching(i).idx(end) = [];
                niching(i).len = niching(i).len - 1;
                if niching(i).len == 0
                    niching(i) = [];
                end
                i = i - 1;
                if i == 0
                    i = length(niching);
                end
            end
            indis = obj.hub.GetIndis(rand_indi(obj.algRand, length(dlist), obj.hub.D, obj.hub.lower, obj.hub.upper));
        end
        
         %% Generate the offspring individuals with the selected parents
        function next = GetNextSubPop(obj, select_idx)
            n = length(select_idx);
            D = obj.hub.D;
            pop_decs = obj.pop(select_idx).decs;
%             pop_fits = obj.pop(select_idx).fits;
%             pop_cons = obj.pop(select_idx).cons;
            matdis = pdist2(pop_decs, pop_decs);
            if n == 1
                c = 1;
            else
                c = max(matdis(:));
            end
            alpha = 0.5 * c / sqrt(D);                  % alpha
            next = [];
            
            for i=1:n
                anti = pop_decs(i, :);
                pool_decs = ones(obj.copy+1, 1) * anti;          % clone pool
%                 pool_fit = ones(obj.copy+1, 1) .* pop_fits(i);   % clone fitnees
%                 pool_con = ones(obj.copy+1, 1) .* pop_cons(i);
                if i == 1 % the seed uses Gaussian mutation
                    pool_decs = pool_decs + [zeros(1, D); alpha .* randn(obj.algRand, obj.copy, D)];
                else
                    if rand(obj.algRand) < 0.5
                        step = [zeros(1, D); 0.5 .* (pop_decs(1, :) - anti); alpha .* randn(obj.algRand, obj.copy-1, D)];
                        pool_decs = pool_decs + step;
                    else
                        pool_decs = pool_decs + [zeros(1, D); alpha .* randn(obj.algRand, obj.copy, D)];
                    end
                end
                pool_decs = boundary_check(pool_decs, obj.hub.lower, obj.hub.upper);      % check the boundary
                pool = [obj.pop(select_idx(i)); obj.hub.GetIndis(pool_decs(2:end, :))];
%                 pool = obj.hub.GetIndis(pool_decs);
                if ~isempty(pool)
                    sort_idx = obj.SortWithContraint(pool);
                    next = [next; pool(sort_idx(1))];
                end
            end
        end

    end
end