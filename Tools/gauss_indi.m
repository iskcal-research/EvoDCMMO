function pop = gauss_indi(algRand, decs, radii, lb, ub)
    pop = decs + repmat(radii, 1, size(decs, 2)) .* randn(algRand, size(decs));

    pop = boundary_check(pop, lb, ub);
end