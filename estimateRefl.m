function [refl] = estimateRefl(ims, mask)
    ims = max(ims, 0.004);
    log_ims = log(ims);
    
    [i_y_all, i_x_all] = Poisson.get_gradients(log_ims);
    r_y = median(i_y_all, 3);
    r_x = median(i_x_all, 3);
    log_refl = Poisson.solve(r_y, r_x, mask);

    refl = mask .* exp(log_refl);
    refl = refl ./ max(refl(mask));
end

