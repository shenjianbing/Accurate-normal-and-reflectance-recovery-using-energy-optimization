function [n_new, l_new] = updateNormCF(albedo, n_init, ims, mask, lit_mask, light, shad)

F = size(ims,4);
P = sum(mask(:));
[h, w] = size(mask);
ims_g = squeeze(mean(ims,3));
ims_g = reshape(ims_g, [h*w F]);
ims_g = ims_g ./ repmat(albedo(:), [1 F]);

shad  = reshape(shad, [h*w F]);
weight  = shad .^ (1/2);
ims_g   = ims_g .* weight;

n   = zeros(P,3);
ind = find(mask);

warning off;
for i = 1:P
    p = ind(i);
    I = ims_g(p,:);
    w = weight(p,:);
    
    L = light * diag(w);
    n(i,:) = I/L;
end
warning on;

n = normr(n);
n_init(mask,:) = n;
n_new = n_init;

n_new(isnan(n_new)) = 0;
n_new(~mask,:) = 0;

% update pixels which always stay in shadow
all_dark = true(size(mask));
for i = 1:F
    all_dark = all_dark & (~lit_mask(:,:,i));
end
n_tmp = n_new;
for i = 2:h-1
    for j = 2:w-1
        if all_dark(i,j)
            n_tmp(i*w+j, :) = 0.25*( ...
                n_new((i-1)*w+j,:) + n_new((i+1)*w+j,:) + ...
                n_new(i*w+j-1,:) + n_new(i*w+j+1,:) ...
            );
        end
    end
end
n_new(mask,:) = normr(n_tmp(mask,:));

% update lighting direction by robust pixels
all_mask = true(size(mask));
for i = 1:size(lit_mask,3)
    all_mask = all_mask & lit_mask(:,:,i);
end
l_new = n_new(all_mask,:) \ ims_g(all_mask,:);
l_new = l_new * 0.05 + light * 0.95;

end
