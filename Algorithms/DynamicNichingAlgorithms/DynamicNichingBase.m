classdef DynamicNichingBase < DynamicAlgorithm
    properties(SetAccess=protected)
        nct;
        dyt;
        cht;

        curG;
        maxG;
    end
    methods(Access = public)
        function obj = DynamicNichingBase(nt, rt, ct)
            obj.nct = struct();
            obj.dyt = struct();
            obj.cht = struct();
            
            obj.nct.type = nt;
            obj.dyt.type = rt;
            obj.cht.type = ct;
            
            obj.curG = 0;
            obj.maxG = round(obj.hub.freq / obj.hub.N);
        end

        function AfterInitial(obj)

            % constraint parameter
            if obj.cht.type == 2
                % ec
                con = obj.pop.cons;
                obj.cht.ec.p = 0.5;
                sort_con = sort(con, 'ascend');
                theta = round(0.9 * obj.hub.N);
                obj.cht.ec.epsilon_0 = sort_con(theta);
                obj.cht.ec.epsilon_g = 0;
                obj.cht.ec.cp = -(log(obj.cht.ec.epsilon_0)+5) / log(1-obj.cht.ec.p);
            elseif obj.cht.type == 3
                obj.cht.pf.r = 10;
            elseif obj.cht.type == 4
                obj.cht.apf.r = 10;
            elseif obj.cht.type == 5
                % sr
                obj.cht.sr.pf = 0.45;
            elseif obj.cht.type ~= 1
                error("invalid constraint number");
            end

            % dynamic parameter
            if obj.dyt.type == 1
                obj.dyt.random.per = 0.1;
            elseif obj.dyt.type == 2
                obj.dyt.elite.per = 0.1;
            elseif obj.dyt.type == 3
                obj.dyt.memory.size = 50;
                obj.dyt.memory.maxNum = 5;
                obj.dyt.memory.pop = [];
            else
                error("invalid dynamic number");
            end

            % niching parameter
            if obj.nct.type == 3
                obj.nct.nbc.phi = 2;
            end
        end
        
        function Evolve(obj)
            obj.curG = obj.curG + 1;
            obj.EvolveParams();

            obj.GetNextPop(); %% used for optimizor
            %% obj.MaintainDiversity(); %% used for the dynamic strategy
        end

        function GetNextPop(obj)
            
        end

        function index = SortWithContraint(obj, pop)
            fit = pop.fits;
            con = pop.cons;
            if obj.cht.type == 1
                data = [con -fit];
                [~, index] = sortrows(data);
            elseif obj.cht.type == 2
                new_con = con;
                new_con(new_con <= obj.cht.ec.epsilon_g) = 0;
                data = [new_con -fit];
                [~, index] = sortrows(data);
            elseif obj.cht.type == 3
                fit_con = fit - con .^2 .* obj.cht.pf.r;
                [~, index] = sort(fit_con, 'descend');
            elseif obj.cht.type == 4
                fit_con = fit - con .^2 .* obj.cht.apf.v.* obj.cht.apf.r;
                [~, index] = sort(fit_con, 'descend');        
            elseif obj.cht.type == 5
                index = stochastic_sort(fit, con, obj.cht.sr.pf, obj.algRand);
            else
                error("the constraint techeque number should be between 1 and 5.");
            end
        end
        
        function niching = GetNiching(obj)
            pop = obj.Record();
            index = obj.SortWithContraint(pop);
            pop = pop(index);
            switch(obj.nct.type)
                case 1 % jp
                    matdis = pdist2(pop.decs, pop.decs);
                    niching = jarvis_partrick(matdis);
                case 2 % ts
                    matdis = pdist2(pop.decs, pop.decs);
                    niching = topographical_selection(matdis, obj.hub.N, obj.hub.D);
                case 3 % nearest better clutering
                    niching = nearest_better_clustering(pop.decs, obj.nct.nbc.phi);
                case 4 % nearest better neighbor clutering
                    niching = nearest_better_neighbor_clustering(pop.decs);
                otherwise
                    error("wrong niching type in dynamic niching algorithms");
            end
%             idxes = ones(2, length(pop));
%             idxes(1, :) = 1: length(pop);
%             idxes(2, :) = index;
            for i = 1:length(niching)
                niching(i).seed = index(niching(i).seed);
                niching(i).idx = index(niching(i).idx);
            end
        end

        function select = CompareWithConstraint(obj, pop, off)
            pop_fit = pop.fits;
            off_fit = off.fits;
            pop_con = pop.cons;
            off_con = off.cons;


            cmp_fit = off_fit > pop_fit; 
            cmp_con = off_con < pop_con; 
        
            if obj.cht.type == 1 % SFS
                cmp_fit_index = (pop_con == 0 & off_con == 0) | (pop_con == off_con); 
                cmp_con_index = ~cmp_fit_index; % ~t
                select = logical(cmp_fit_index .* cmp_fit + cmp_con_index .* cmp_con);
            elseif obj.cht.type == 2 % EC
                cmp_fit_index = (pop_con <= obj.cht.ec.epsilon_g & off_con <= obj.cht.ec.epsilon_g) | (pop_con == off_con);
                cmp_con_index = ~cmp_fit_index;
                select = logical(cmp_fit_index .* cmp_fit + cmp_con_index .* cmp_con);
            elseif obj.cht.type == 3 % PF
                off_fit_con = off_fit - off_con .^2 .* obj.cht.pf.r;
                pop_fit_con = pop_fit - pop_con .^2 .* obj.cht.pf.r;
                select = off_fit_con > pop_fit_con;
            elseif obj.cht.type == 4 %DPF
                % Varying Fitness Functions in Genetic Algorithms: Studying the Rate of Increase of the Dynamic Penalty Terms
                
                off_fit_con = off_fit - off_con .^2 .* obj.cht.apf.v.* obj.cht.apf.r;
                pop_fit_con = pop_fit - pop_con .^2 .* obj.cht.apf.v.* obj.cht.apf.r;
                select = off_fit_con > pop_fit_con;     
            elseif obj.cht.type == 5 % SS
                cmp_fit_index = pop_con == 0 & off_con == 0;
                cmp_con_count = sum(~cmp_fit_index);
                cmp_fit_index(~cmp_fit_index) = rand(obj.algRand, cmp_con_count, 1) < obj.cht.sr.pf;
                cmp_con_index = ~cmp_fit_index;
                select = logical(cmp_fit_index .* cmp_fit + cmp_con_index .* cmp_con);
            else
                error("invalid constraint type number");
            end
        end

        %% used for random and elite immigrants
        function [indis, index] = MaintainDiversity(obj, pop)
            rest = obj.hub.freq - rem(obj.hub.evaluated, obj.hub.freq);
            if rest == obj.hub.freq
                indis = [];
                index = [];
                return;
            end
%             disp(obj.hub.evaluated);
            
            if obj.dyt.type == 1 % random immigrants
                sort_index = obj.SortWithContraint(pop);
                index = sort_index(obj.hub.N-min(obj.hub.N*obj.dyt.random.per, rest)+1:end);
                indis = obj.hub.GetIndis(rand_indi(obj.algRand, length(index), obj.hub.D, obj.hub.lower, obj.hub.upper));
            elseif obj.dyt.type == 2 % elite immigrants
                sort_index = obj.SortWithContraint(pop);
                niching = obj.GetNiching();
                nums = min([length(niching), round(obj.hub.N*obj.dyt.elite.per), rest]);
                index = sort_index(obj.hub.N-nums+1:end);
                
                radii = zeros(nums, 1);
                for i = 1:nums
                    sub_index = niching(i);
                    radii(i) = max(pdist2(pop(sub_index.seed).dec, pop(sub_index.idx).decs))/(2*sqrt(obj.hub.D));
                end
                elite = pop(cat(1,niching.seed));
                elite = elite(1:nums);
                indis = obj.hub.GetIndis(gauss_indi(obj.algRand, elite.decs, radii, obj.hub.lower, obj.hub.upper));
            else
                indis = [];
                index = [];
            end
            assert(length(indis) == length(index));
        end

        function EvolveParams(obj)
            if obj.cht.type == 2
                if obj.curG / obj.maxG < obj.cht.ec.p 
                    obj.cht.ec.epsilon_g = obj.cht.ec.epsilon_0 * (1 - obj.curG/obj.maxG) ^ obj.cht.ec.cp;
                else
                    obj.cht.ec.epsilon_g = 0;
                end
            elseif obj.cht.type == 4
                obj.cht.apf.v = 1 - exp(-10*obj.curG / obj.maxG);
            end
        end
        
        %% used for memory immigrants
        function RespondChange(obj)
            if obj.dyt.type == 3
                niching = obj.GetNiching();

                seeds_idx = [];
                for i = 1:min(length(niching), obj.dyt.memory.maxNum)
                    seeds_idx = [seeds_idx niching(i).seed];
                end
                %[obj.pop, obj.dyt.memory.pop] = use_memory(obj.algRand, obj.Record(), obj.dyt.memory.pop, obj.dyt.memory.size, seeds_idx, obj.hub.N, obj.hub.D, obj.hub.lower, obj.hub.upper);
                
                pop = obj.Record();
                mem = obj.dyt.memory.pop;
                
                pop = obj.hub.GetIndis(pop.decs);
                if ~isempty(mem)
                    mem = obj.hub.GetIndis(mem.decs);
                    mem_sort_idx = obj.SortWithContraint(mem);
                    mem = mem(mem_sort_idx);
                end
                
                temp = pop(seeds_idx);
                sort_idx = obj.SortWithContraint(pop);
                pop = pop(sort_idx);

                i = round(obj.hub.N / 2);
                select_num = min([length(seeds_idx), 5, length(mem)]);
                if ~isempty(mem)
                    pop(i+1:i+select_num) = mem(1:select_num);
                    i = i + select_num;
                end
                pop(i+1:end) = obj.hub.GetIndis(rand_indi(obj.algRand, obj.hub.N - i, obj.hub.D, obj.hub.lower, obj.hub.upper));
                obj.pop = pop;
                
                for i = length(temp):-1:1
                    if length(mem) < obj.dyt.memory.size
                        mem = [mem; temp(i)];
                    else
                        dis = pdist2(temp(i).dec, mem.decs);
                        [~, nearest_idx] = min(dis);
                        if obj.CompareWithConstraint(mem(nearest_idx), temp(i))%temp(i).fit > mem(nearest_idx).fit
                            mem(nearest_idx) = temp(i);
                        end
                    end
                end
                obj.dyt.memory.pop = mem;
            end

            obj.curG = 0;
            
            if obj.cht.type == 2
                sort_con = sort(obj.pop.cons, 'ascend');
                theta = round(0.9 * obj.hub.N);
                obj.cht.ec.epsilon_0 = sort_con(theta);
                obj.cht.ec.epsilon_g = 0;
                obj.cht.ec.cp = -(log(obj.cht.ec.epsilon_0)+5) / log(1-obj.cht.ec.p);
            end
        end
    end
end