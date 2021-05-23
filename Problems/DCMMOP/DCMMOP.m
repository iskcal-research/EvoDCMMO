classdef DCMMOP < DynamicProblem
    properties(SetAccess=private)
        %proRand;                        % rand stream
        %D;                              % dimension (D)
        G;                              % number of the global peaks (G)
        GArr;                           % range of the global peaks (GArr)
        L;                              % number of the local peaks (L)
        num_peaks;                     % total number of the peaks (G+L)
        s = 1.0;                        % shift length (s)
        lb = -100;                        % lower boundary (lb)
        ub = 100;                         % upper boundary (ub)
        %num_change = 60;                % the number of the environment change (Dy)
        %freq = 3000 * D;                % the maximum evoluations of an environment (freq)
        %maxFEs = num_change * freq;     % the total evolutions (MaxFEs)
        dpeaks = 0.1;                   % the minimum distance between the peaks (dpeaks)
        
        h;
        w;
        height = 50;                    % initial height (height)
        lh_min = 5;                     % minimal height difference (dif_{min})
        lh_max = 20;                    % maximal height difference (dif_{max})
        h_min = 30;                     % minimal height (h_{min})
        h_max = 70;                     % maximal height (h_{max})
        h_s = 7.0;                      % change severity of height (h_s)
        w_min = 1;                      % minimal width (w_{min})
        w_max = 12;                     % maximal width (w_{max})
        w_s = 1.0;                      % change severity of width (w_s)
        
        ch;
        cw;
        phi;
        phi_f;
        cr;
        P;                              % the number of first type of feasible regions (P)
        PArr;                           % the range of first type of feasible regions (PArr)
        R;                              % the number of second type of feasible regions (R)
        cnum_peaks;                     % the total number of feasible regions (P+R)
        cheight = 50;                   % initial constraint height (height^c)
        ch_min = 30;                    % minimal constraint height (h_{min}^c)
        ch_max = 70;                    % maximal constraint height (h_{max}^c)
        ch_s = 7.0;                     % change severity of constraint height (h_s^c)
        cw_min = 1;                     % minimal constraint width (w_{min}^c)
        cw_max = 12;                    % maximal constraint width (w_{max}^c)
        cw_s = 1.0;                     % change severity of constraint width (w_s^c)
        cv = 100;                         % constraint control value (cv)
        cfv = 100;
    
        rest;                           % the rest evaluations in the current envirnment
        X;                              % positions of fitness peaks
        cX;                             % positions of constraint peaks
    end 
    methods
        function obj = DCMMOP(D, G, L, P, R)
            %obj.proRand = RandStream('mt19937ar','Seed', 0);
            %obj.hub.D = D;
            if ~isscalar(G)
                obj.GArr = G;
                G = round((G(1)+G(2))/2);
            end
            if ~isscalar(P)
                obj.PArr = P;
                P = round((P(1)+P(2))/2);
            end
            
            if P > G
                R = R + (P - G);
                P = G;
            end
            
            obj.G = G;
            obj.L = L;
            obj.num_peaks = G + L;
            obj.P = P;
            obj.R = R;
            obj.cnum_peaks = P + R;
            
            obj.lb = ones(1, D) .* obj.lb;
            obj.ub = ones(1, D) .* obj.ub;
            
            obj.hub.D = D;
            obj.hub.lower = obj.lb;
            obj.hub.upper = obj.ub;
            obj.hub.freq = 3000 * D;
            obj.hub.evaluation = 60 * obj.hub.freq;
            obj.rest = obj.hub.freq;
        end
        
        function Initialize(obj)
            obj.X = initialize_peaks_with_gap(obj.proRand, obj.num_peaks, obj.hub.D, obj.lb, obj.ub, obj.dpeaks);
            
            obj.h = ones(1, obj.num_peaks) .* obj.height - [zeros(1, obj.G), rand(obj.proRand, 1, obj.L) .* (obj.lh_max-obj.lh_min) + obj.lh_min];
            obj.w = rand(obj.proRand, 1, obj.num_peaks) * (obj.w_max-obj.w_min)+obj.w_min;

            obj.ch = ones(1, obj.cnum_peaks) .* obj.cheight;
            obj.cw = rand(obj.proRand, 1, obj.cnum_peaks) .* (obj.cw_max-obj.cw_min)+obj.cw_min;
            obj.phi = min(obj.ch) / obj.cv;                  % the parameter \delta(t) in the constraint landscape
            obj.cr = sqrt((obj.ch ./ obj.phi - 1) ./ obj.cw);  % the radius of feasible regions
            obj.phi_f = min(obj.h) - obj.cfv;
            
            % generate the constraint peaks location
            select_peaks = randperm(obj.proRand, obj.G, obj.P);
            select_lb = obj.X(select_peaks, :) - repmat(obj.cr(1:obj.P)', 1, obj.hub.D) ./ sqrt(obj.hub.D);
            select_ub = obj.X(select_peaks, :) + repmat(obj.cr(1:obj.P)', 1, obj.hub.D) ./ sqrt(obj.hub.D);

            obj.cX = rand(obj.proRand, obj.P, obj.hub.D) .* (select_ub - select_lb) + select_lb;
            obj.cX = [obj.cX; rand(obj.proRand, obj.R, obj.hub.D) .* (obj.ub - obj.lb) + obj.lb];
        end
        
        function fits = GetFits(obj, decs)
            dist = pdist2(decs, obj.X);
            fits = max(obj.h - obj.w .* dist, [], 2);
        end
        
        function cons = GetCons(obj, decs)
            cdist = pdist2(decs, obj.cX);
            cons1 = obj.phi - max(obj.ch ./ (1 + obj.cw .* (cdist .^ 2)), [], 2);
            cons2 = obj.phi_f - obj.GetFits(decs);
            cons = max(cons1, 0) + max(cons2, 0);
%             cons = max(cons1, 0);
        end
        
        function bestfit = GetBestFit(obj)
            bestfit = max(obj.h);
        end
        
        function bestpos = GetBestPos(obj)
            cons = obj.GetCons(obj.X(1:obj.G, :));
            
            pos_idx = cons == 0;
            bestpos = obj.X(pos_idx, :);
        end
        
        function found = GetFoundPeak(obj, pop, epsilons)
            best_fit = obj.GetBestFit();
            best_pos = obj.GetBestPos();
            dist = pdist2(pop.decs, best_pos);
            
            cons = obj.GetCons(obj.X(1:obj.G, :));
            pos_idx = cons == 0;
            
            found = zeros(1, length(epsilons));
            
            for i = 1:length(epsilons)
                [~, peak_idx] = max(obj.h(pos_idx) - obj.w(pos_idx) .* dist, [], 2);
                select_idx = (abs(best_fit - pop.fits) <= epsilons(i)) & (pop.cons == 0);
                peak_idx = peak_idx(select_idx);
    
                found(i) = length(unique(peak_idx));
            end
        end
        
        function ChangeDynamic(obj)
            obj.ChangeNumsOfPeaks();            
            
            obj.X = move_peaks_with_gap(obj.proRand, obj.X, obj.num_peaks, obj.hub.D, obj.lb, obj.ub, obj.dpeaks, obj.s);
            
            % update the height
            obj.height = obj.height + randn(obj.proRand, 1, 1) * obj.h_s;
            obj.height = max(obj.height, obj.h_min+obj.lh_min);
            obj.height = min(obj.height, obj.h_max);

            obj.h = ones(1, obj.num_peaks) * obj.height - [zeros(1, obj.G), rand(obj.proRand, 1, obj.L) * (obj.lh_max-obj.lh_min) + obj.lh_min];
            obj.h = max(obj.h, obj.h_min);
            obj.h = min(obj.h, obj.h_max);

            % update the width
            obj.w = obj.w + randn(obj.proRand, 1, obj.num_peaks) * obj.w_s;
            obj.w = max(obj.w, obj.w_min);
            obj.w = min(obj.w, obj.w_max);

            % update the constraints
            obj.ch = obj.ch + randn(obj.proRand, 1, obj.cnum_peaks) * obj.ch_s;
            obj.ch = max(obj.ch, obj.ch_min);
            obj.ch = min(obj.ch, obj.ch_max);
            obj.phi = min(obj.ch)/obj.cv;
            obj.phi_f = min(obj.h) - obj.cfv;
            obj.cw = obj.cw + randn(obj.proRand, 1, obj.cnum_peaks) * obj.cw_s;
            obj.cw = max(obj.cw, obj.cw_min);
            obj.cw = min(obj.cw, obj.cw_max);
            obj.cr = sqrt((obj.ch ./ obj.phi - 1) ./ obj.cw); 

            % generate the constraint peaks location
            select_peaks = randperm(obj.proRand, obj.G, obj.P);
            select_lb = obj.X(select_peaks, :) - repmat(obj.cr(1:obj.P)', 1, obj.hub.D) ./ sqrt(obj.hub.D);
            select_ub = obj.X(select_peaks, :) + repmat(obj.cr(1:obj.P)', 1, obj.hub.D) ./ sqrt(obj.hub.D);

            obj.cX = rand(obj.proRand, obj.P, obj.hub.D) .* (select_ub - select_lb) + select_lb;
            obj.cX = [obj.cX; rand(obj.proRand, obj.R, obj.hub.D) .* (obj.ub - obj.lb) + obj.lb];
        end
        
        function ChangeNumsOfPeaks(obj)
            % change of the number of peaks and constrained peaks
            if ~isempty(obj.GArr)
                r_dir = rand(obj.proRand);
                if r_dir < 0.4
                    g = obj.G + 1;
                elseif r_dir < 0.8
                    g = obj.G - 1;
                else
                    g = obj.G;
                end
                g = max(obj.GArr(1), g);
                g = min(obj.GArr(2), g);
                if g > obj.G
                    rand_w = rand(obj.proRand) * (obj.w_max-obj.w_min)+obj.w_min;
                    rand_x = rand(obj.proRand, 1, obj.hub.D) .* (obj.hub.upper - obj.hub.lower) + obj.hub.lower;
                    obj.w = [obj.w(1:obj.G) rand_w obj.w(obj.G+1:end)];
                    obj.X = [obj.X(1:obj.G, :); rand_x; obj.X(obj.G+1:end, :)];
                elseif g < obj.G
                    obj.w = [obj.w(1:obj.G-1) obj.w(obj.G+1:end)];
                    obj.X = [obj.X(1:obj.G-1,:); obj.X(obj.G+1:end, :)];
                end
                obj.G = g;
                obj.num_peaks = obj.G + obj.L;
            end
            
            if ~isempty(obj.PArr)
                r_dir = rand(obj.proRand);
                if r_dir < 0.4
                    p = obj.P + 1;
                elseif r_dir < 0.8
                    p = obj.P - 1;
                else
                    p = obj.P;
                end
                p = max(obj.PArr(1), p);
                p = min([obj.PArr(2), p, obj.G]);
                if p > obj.P
                    obj.ch = [obj.ch(1:obj.P) obj.cheight obj.ch(obj.P+1:end)];
                    rand_cw = rand(obj.proRand) .* (obj.cw_max-obj.cw_min)+obj.cw_min;
                    obj.cw = [obj.cw(1:obj.P) rand_cw obj.cw(obj.P+1:end)];
                    
                    obj.cX = [obj.cX(1:obj.P, :); zeros(1, obj.hub.D); obj.cX(obj.P+1:end, :)];
                elseif p < obj.P
                    obj.ch = [obj.ch(1:obj.P-1) obj.ch(obj.P+1:end)];
                    obj.cw = [obj.cw(1:obj.P-1) obj.cw(obj.P+1:end)];
                    obj.cX = [obj.cX(1:obj.P-1, :); obj.cX(obj.P+1:end, :)];
                end
                obj.P = p;
                obj.cnum_peaks = obj.P + obj.R;
            end
        end
    end
end