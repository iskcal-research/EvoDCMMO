classdef DynamicNichingDE < DynamicNichingBase
    methods
        function obj = DynamicNichingDE(nt, rt, ct)
            obj = obj@DynamicNichingBase(nt, rt, ct);
        end
        
        function GetNextPop(obj)
            niching = obj.GetNiching();
            
            for i = 1:length(niching)
                sub_index = niching(i).idx;
                
                if length(sub_index) >= 4
                    n = length(sub_index);
                    idx = DE_rand_idx(obj.algRand, n, n, 3);     % DE/rand/1
                    pop_decs = obj.pop(sub_index).decs;
                   
                    off_decs = pop_decs(idx(:, 1), :) + 0.5 * (pop_decs(idx(:, 2), :) - pop_decs(idx(:, 3), :));
                    off_decs = boundary_check(off_decs, obj.hub.lower, obj.hub.upper);

                    off_decs = crossover(obj.algRand, off_decs, pop_decs, 0.9);
                    offspring = obj.hub.GetIndis(off_decs);
                    
                    if ~isempty(offspring)
                        sub_index = sub_index(1:length(offspring));
                        parent = obj.pop(sub_index);
                        select = obj.CompareWithConstraint(parent, offspring);

                        obj.pop(sub_index(select)) = offspring(select);
                    end
                end
            end
            
            [indis, index] = obj.MaintainDiversity(obj.pop);
            if ~isempty(index)
%                 disp(obj.hub.evaluated);
                obj.pop(index) = indis;
            end
        end
    end
end