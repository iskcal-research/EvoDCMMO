classdef DynamicNichingPSO < DynamicNichingBase
    properties(SetAccess=protected)
        pbest;
        w = 0.729;
        c1 = 2.05;
        c2 = 2.05;
    end
    
    methods
        function obj = DynamicNichingPSO(nt, rt, ct)
            obj = obj@DynamicNichingBase(nt, rt, ct);
        end
        
        function pop = Record(obj)
            pop = obj.pbest;
        end
        
        %% Initialize the population
        function Initialize(obj)
            obj.pop = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper), struct("v", zeros(obj.hub.N, obj.hub.D)));
            obj.pbest = obj.pop;
        end
        
        function GetNextPop(obj)
%             disp(obj.hub.evaluated);
            niching = obj.GetNiching();
            
%             select = obj.CompareWithConstraint(obj.pbest, obj.pop);
%             assert(all(select==0));
            
            for i = 1:length(niching)
                sub_index = niching(i).idx;
                off = obj.GetOffspring(sub_index, niching(i).seed);
                if ~isempty(off)
                    select_index = sub_index(1:length(off));
                    obj.pop(select_index) = off;
                    select = obj.CompareWithConstraint(obj.pbest(select_index), off);
                    obj.pbest(select_index(select)) = off(select);
                end
            end

%             select = obj.CompareWithConstraint(obj.pbest, obj.pop);
%             assert(all(select==0));
%             obj.pbest(select) = obj.pop(select);
            
            [indis, index] = obj.MaintainDiversity(obj.pbest);
            if ~isempty(index)
%                 disp(obj.hub.evaluated);
                indis.SetAdds(struct("v", obj.pop(index).adds("v")));
                obj.pop(index) = indis;
            end
            
            select = obj.CompareWithConstraint(obj.pbest, obj.pop);
            obj.pbest(select) = obj.pop(select);
        end
        
         %% Generate the offspring individuals
        function offspring = GetOffspring(obj, select_idx, seed)
%             disp(obj.hub.evaluated);
            N = length(select_idx);
            pbest_decs = obj.pbest(select_idx).decs;
            pop_decs = obj.pop(select_idx).decs;
            v = obj.pop.adds("v");
            v = v(select_idx, :);
            gbest_decs = repmat(obj.pbest(seed).dec, N, 1);
            v = obj.w .* (v + obj.c1 .* rand(obj.algRand, N, obj.hub.D) .* (pbest_decs - pop_decs) + obj.c2 .* rand(obj.algRand, N, obj.hub.D) .* (gbest_decs - pop_decs));
            off_decs = pop_decs + v;
            off_decs_new = boundary_check(off_decs, obj.hub.lower, obj.hub.upper);
            v = off_decs_new - pop_decs;            
            offspring = obj.hub.GetIndis(off_decs_new, struct("v", v));
        end
        
        %% dynamic response 
        function RespondChange(obj)
            obj.pop = obj.pbest;
            RespondChange@DynamicNichingBase(obj);
            obj.pop.SetAdds(struct("v", zeros(obj.hub.N, obj.hub.D)));
            obj.pbest = obj.pop;
        end
    end
end