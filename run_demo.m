
clear;
addpath('MinimizeEntropy', 'DiffuseMaxima', 'GMM-GMR-v2.0');

tag = 'panther';
maskname = 'mask.png';

fig_path = ['data/' tag '/'];
mask = logical(imread([fig_path maskname]));
mask = mask(:,:,1);
[h, w] = size(mask);

% 去除边界非mask值
mask(1,:) = false;
mask(:,1) = false;
mask(h,:) = false;
mask(:,w) = false;

% 使用11幅图像，包括original.png
files = dir([fig_path 'light*.png']);
files = cat(1, files, dir([fig_path 'original.png']));
F = length(files);

ims = zeros([h w 3 F]);
I = zeros(h*w, F);
for i = 1:F
    fname = [fig_path files(i).name];
    im = im2double(imread(fname));
    ims(:,:,:,i) = im;
    
    im = mean(im,3);
    I(:,i) = im(:);
end

sign = -1;
[n_dm, l_dm, r_dm] = photometricStereo(I, mask, 'dm', sign);
% Auxiliary.showN(n_dm, mask, 'dm');
% Auxiliary.showR(r_dm, mask, ims, 'dm');

tr = 0.006;
ts = 0.008;
[n_our, l_our, r_our] = estimateNR(ims, mask, n_dm, r_dm, l_dm, tr, ts);
Auxiliary.showN(n_our*sign, mask, 'normal', 1);
Auxiliary.showR(r_our, mask, ims, 'reflectance', 1);

res_path = [tag '_result\'];
if ~exist(res_path, 'dir')
    mkdir(res_path);
end
Auxiliary.saveR(r_our, mask, ims, [res_path 'reflectance.png'], 1, 0);
Auxiliary.saveN(n_our*sign, mask, [res_path 'normal.png'], 1, 0);

shad04 = mean(ims(:,:,:,4), 3) ./ r_our;
shad04(~mask) = 1;
shad04(mask) = mat2gray(shad04(mask));
imwrite(shad04, [res_path 'shading04.png']);
figure, imshow(shad04);
