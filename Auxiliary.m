classdef Auxiliary < handle
    
    methods(Static)
        function showR(r, mask, ims, ti, padding)
            if nargin < 4
                ti = 'reflectance';
            end
            if nargin < 5
                padding = 0;
            end
            
            show(rMap(r, mask, ims, padding), ti);
        end
        
        function showN(n, mask, ti, padding)
            if nargin < 3
                ti = 'normal';
            end
            if nargin < 4
                padding = 0;
            end
            
            show(nMap(n, mask, padding), ti);
        end
        
        function saveR(r, mask, ims, fname, padding, margin)
            if nargin < 5
                padding = 0;
                margin = 0;
            end
            
            refl = rMap(r, mask, ims, padding);
            if margin > 0
                refl = Auxiliary.clip(refl, mask, margin);
            end
            imwrite(refl, fname);
        end
        
        function saveN(n, mask, fname, padding, margin)
            if nargin < 4
                padding = 0;
                margin = 0;
            end
            
            normal = nMap(n, mask, padding);
            if margin > 0
                normal = Auxiliary.clip(normal, mask, margin);
            end
            imwrite(normal, fname);
        end
        
        function [clipped] = clip(im, mask, margin)
            if nargin < 3
                margin = 3;
            end
            [row, col] = ind2sub(size(mask), find(mask));
            [h, w] = size(mask);
            
            row_start = max(1, min(row)-margin);
            row_end = min(h, max(row)+margin);
            col_start = max(1, min(col)-margin);
            col_end = min(w, max(col)+margin);
            clipped = im(row_start:row_end, col_start:col_end, :);
        end
        
        function saveSR(est_s, est_r, tag, prefix, estimatorName)
            if ~exist(prefix, 'dir')
                mkdir(prefix);
            end
            mask = load_object(tag, 'mask');   
            [ims] = load_multiple(tag, mask);
            
            est_s(mask) = mat2gray(est_s(mask));
            est_s(~mask) = 1;
            shad = Auxiliary.clip(est_s, mask);
            s_name = [prefix tag '_' estimatorName 's.png'];
            imwrite(shad, s_name);
            
            est_r(mask) = mat2gray(est_r(mask));
            img_chroma = sum(ims, 4);
            img_chroma = img_chroma ./ repmat(max(eps,mean(img_chroma,3)),[1,1,size(img_chroma,3)]);
            albedo = repmat(est_r,[1 1 3]).*img_chroma;
            albedo( repmat(~mask,[1 1 3]) ) = 1;
            
            refl = Auxiliary.clip(albedo, mask);
            r_name = [prefix tag '_' estimatorName 'r.png'];
            imwrite(refl, r_name);
        end
        
        function [score] = score_image(gt_s, gt_r, est_s, est_r, mask)
            ncorr = Metric.ncorr(gt_s, gt_r, est_s, est_r, mask);
            score = 0.5*(ncorr(1)+ncorr(2));
        end
        
    end
end

function show(im, ti)
    figure, imshow(im), title(ti);
end

function [n_map] = nMap(n, mask, padding)
    [h w] = size(mask);
    n = normr( reshape(n, [h*w 3]) );
    n_map = uint8((n+1)*128);
    n_map(~mask,:) = 255*padding;
    n_map = reshape(n_map, [h w 3]);
end

function [r_map] = rMap(r, mask, ims, padding)
    chroma = sum(ims,4);
    chroma = chroma ./ repmat(max(eps, mean(chroma,3)), [1 1 size(chroma,3)]);
    r_map = repmat(r, [1 1 3]) .* chroma;
    r_map( repmat(~mask,[1 1 3]) ) = padding;
end


