function [n_res l_res r_res] = estimateNR(ims, mask, n_init, r_init, l_init, tr, ts)

ims_gray    = squeeze(mean(ims,3));
refl        = estimateRefl(ims_gray, mask);
shad        = zeros(size(ims_gray));
lit_mask    = logical(false(size(ims_gray)));

% shadow field
shading = zeros(size(mask));
for i = 1:size(ims_gray,3)
    im = ims_gray(:,:,i);
    shading(mask) = im(mask) ./ refl(mask);
    
    shading(mask) = shading(mask)/max(shading(mask));
    theta = getThreshold(shading, mask);    
    
    lit_mask(:,:,i) = shading > theta;
    shad(:,:,i) = shading;
end

lighting = l_init;
albedo = r_init;
normal = n_init;

display('start iteration ...')
niter = 2;
for iter = 1:niter
    iter
    albedo = updateReflCF(albedo, normal, ims, mask, lit_mask, tr, ts);
    [normal, lighting] = updateNormCF(albedo, normal, ims, mask, lit_mask, lighting, shad);
end

n_res = normal;
l_res = lighting;
r_res = albedo;
