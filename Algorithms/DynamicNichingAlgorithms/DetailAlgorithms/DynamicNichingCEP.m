classdef DynamicNichingCEP < DynamicNichingBase
   properties(SetAccess=protected)
        tau;
        tau1;
        q = 10;
    end
    methods
        function obj = DynamicNichingCEP(nt, rt, ct)
            obj = obj@DynamicNichingBase(nt, rt, ct);
        end
        
        function Initialize(obj)
            obj.tau = 1/(sqrt(2*sqrt(obj.hub.D)));
            obj.tau1 = 1/(sqrt(2*obj.hub.D));
            obj.pop = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper), struct("eta", 3 * ones(obj.hub.N, obj.hub.D)));
        end
        
        function GetNextPop(obj)
%             disp(obj.hub.evaluated);
            niching = obj.GetNiching();
            
            for i = 1:length(niching)
                sub_index = niching(i).idx;
                off = obj.GetOffspring(sub_index);
                if ~isempty(off)
                    par = obj.pop(sub_index);
                    obj.pop(sub_index(1:length(par))) = obj.Select(par, off);                    
                end
            end
            
            [indis, index] = obj.MaintainDiversity(obj.pop);
            if ~isempty(index)
%                 disp(obj.hub.evaluated);
                indis.SetAdds(struct("eta", obj.pop(index).adds("eta")));
                obj.pop(index) = indis;
            end
        end
        
        function offspring = GetOffspring(obj, select_idx)
            pop_decs = obj.pop(select_idx).decs;
            eta = obj.pop(select_idx).adds("eta");
            
            off_decs = pop_decs + eta .* randn(obj.algRand, size(pop_decs));
            off_decs = boundary_check(off_decs, obj.hub.lower, obj.hub.upper);
            eta_new = eta .* exp(obj.tau1 .* repmat(randn(obj.algRand, size(pop_decs, 1), 1), 1, obj.hub.D) +  obj.tau .* randn(obj.algRand, size(pop_decs)));
            offspring = obj.hub.GetIndis(off_decs, struct("eta", eta_new));
        end
        
        function next_pop = Select(obj, parent, offspring)
            indis = [parent; offspring];
            index = obj.SortWithContraint(indis);
            indis = indis(index);
            s = length(indis);
            n = length(parent);
            k = min(s, obj.q);
            
            win = zeros(s, 2);
            win(:, 2) = s:-1:1;

            for i = 1:s
               select_idx = randperm(obj.algRand, s, k);
               compared = indis(select_idx);
               original = repmat(indis(i), k, 1);
               select = obj.CompareWithConstraint(compared, original);
               win(i, 1) = sum(select);
            end
            [~, idx] = sortrows(win, 'descend');

            next_pop = indis(idx(1:n));
        end
        
         %% dynamic response 
        function RespondChange(obj)
            RespondChange@DynamicNichingBase(obj);
            obj.pop.SetAdds(struct("eta", 3 * ones(obj.hub.N, obj.hub.D)));
        end
    end
end