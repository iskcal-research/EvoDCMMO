function X = move_peaks_with_gap(proRand, X, num_peaks, D, lb, ub, dpeaks, s)
    for i = 1:num_peaks
        accept = false;
        while ~accept
            v = rand(proRand, 1, D) - 0.5;
            len = sqrt(sum(v.^2, 2));
            if len ~= 0
                len = s / len;
            end
            v = v * len;
            X(i,:) = X(i, :) + v;
            X(i,:) = boundary_check(X(i,:), lb, ub);
            
            % check the distance
            if i == 1 
                accept = true;
            else
                dis = pdist2(X(i, :), X(1:i-1, :));
                if all(dis(:) > dpeaks) % accept the position
                    accept = true;
                end
            end
        end
    end
end