function index = stochastic_sort(fit, con, pf, algRand)
    N = size(fit, 1);
    index = 1:N;
    for i = 1:N-1
        for j = i+1:N
            if (((con(i) == 0) && (con(j) == 0)) || (rand(algRand) < pf))
                if fit(i) < fit(j)
                    index([i j]) = index([j i]);
                    fit([i j]) = fit([j i]);
                    con([i j]) = con([j i]);
                end
            else
                if con(i) > con(j)
                    index([i j]) = index([j i]);
                    fit([i j]) = fit([j i]);
                    con([i j]) = con([j i]);
                end
            end
        end
    end
end