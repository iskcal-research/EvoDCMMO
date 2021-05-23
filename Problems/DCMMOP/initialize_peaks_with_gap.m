function X = initialize_peaks_with_gap(proRand, num_peaks, D, lb, ub, dpeaks)
    X = zeros(num_peaks, D);              
    for i = 1:num_peaks
       accept = false;
       while ~accept
           X(i, :) = rand(proRand, 1, D) .* (ub - lb) + lb;
           if i == 1
               accept = true;
           else
               dis = pdist2(X(i, :), X(1:i-1, :));
               if all(dis(:) > dpeaks)
                   accept = true;
               end
           end
       end
    end
end