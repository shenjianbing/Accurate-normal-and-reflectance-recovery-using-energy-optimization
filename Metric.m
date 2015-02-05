classdef Metric < handle
    methods(Static)
        function [mean_angle, angle] = angularN(gt_n, est_n, mask)
            est_n = normr(est_n(mask,:));
            gt_n  = normr(gt_n(mask,:));
            
            inner = abs(sum(est_n.*gt_n,2));
            inner(inner > 1) = 1;
            angle = acosd(inner);
            angle(isnan(angle)) = 0;
            mean_angle = mean(angle);
            
            tmp = zeros(size(mask));
            tmp(mask) = angle;
            angle = tmp;
        end
        
        function [mean_angle, angle] = angularL(gt_l, est_l)
            gt_l = normc(gt_l);
            est_l = normc(est_l);
            inner = abs(sum(gt_l.*est_l, 1));
            inner(inner > 1) = 1;
            angle = acosd(inner);
            mean_angle = mean(angle);
        end
        
        function err = mse(gt, est, mask)
            ssq = ssq_error(gt, est, mask);
            total = sum(mask(:) .* gt(:) .^ 2);
            assert(total > eps);
            err = ssq / total;
        end
        
        function [err] = lmse(gt_s, gt_r, est_s, est_r, mask)
            window_size = 20;
            s_err = local_error(gt_s, est_s, mask, window_size, floor(window_size/2));
            r_err = local_error(gt_r, est_r, mask, window_size, floor(window_size/2));
            err = [s_err r_err];
        end
        
        function [err] = almse(gt_s, gt_r, est_s, est_r, mask)
            window_size = 20;
            s_err = alocal_error(gt_s, est_s, mask, window_size, floor(window_size/2));
            r_err = alocal_error(gt_r, est_r, mask, window_size, floor(window_size/2));
            err = [s_err r_err];
        end
        
        function [err] = bothmse(gt_s, gt_r, est_s, est_r, mask)
            ssq_s   = ssq_error(gt_s, est_s, mask);
            total_s = sum(mask(:) .* gt_s(:) .^ 2);
            ssq_r   = ssq_error(gt_r, est_r, mask);
            total_r = sum(mask(:) .* gt_r(:) .^ 2);
            
            assert(total_s > eps && total_r > eps);
            % err = 0.5*ssq_s/total_s + 0.5*ssq_r/total_r;
            s_err = ssq_s/total_s;
            r_err = ssq_r/total_r;
            err = [s_err r_err];
        end
        
        function [err] = ncorr(gt_s, gt_r, est_s, est_r, mask)
            s_err = 1 - corr2(gt_s(mask), est_s(mask));
            r_err = 1 - corr2(gt_r(mask), est_r(mask));
            err = [s_err r_err];
            % cor = 0.5*corr2(gt_s(mask), est_s(mask)) + 0.5*corr2(gt_r(mask), est_r(mask));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following functions are error metrics defined in:
%
%   Ground truth dataset and baseline evaluations for intrinsic image algorithms.
%   Roger Grosse et. al. iccv 2009 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [score] = score_image(true_shading, true_refl, estimate_shading, estimate_refl, mask, window_size)
    if nargin < 6
        window_size = 20;
    end
    score = 0.5 * local_error(true_shading, estimate_shading, mask, window_size, floor(window_size/2)) + ...
            0.5 * local_error(true_refl, estimate_refl, mask, window_size, floor(window_size/2));
end

function [error_val] = local_error(correct, estimate, mask, window_size, window_shift)
% Returns the sum of the local sum-squared-errors, where the estimate may
% be rescaled within each local region to minimize the error. The windows are
% window_size x window_size, and they are spaced by window_shift.
    M = size(correct,1);
    N = size(correct,2);

    ssq = 0;
    total = 0;
    for i = 1:window_shift:(M - window_size + 1)
        for j = 1:window_shift:(N - window_size + 1)
            correct_curr = correct(i:i+window_size-1, j:j+window_size-1);
            estimate_curr = estimate(i:i+window_size-1, j:j+window_size-1);
            mask_curr = mask(i:i+window_size-1, j:j+window_size-1);
            ssq = ssq + ssq_error(correct_curr, estimate_curr, mask_curr);
            total = total + sum(mask_curr(:) .* correct_curr(:) .^ 2);
        end
    end

    assert(~isnan(ssq/total));
    error_val = ssq / total;
end

function [error_val] = ssq_error(correct, estimate, mask)
% Compute the sum-squared-error for an image, where the estimate is
% multiplied by a scalar which minimizes the error. Sums over all pixels
% where mask is True. If the inputs are color, each color channel can be
% rescaled independently.
    assert(ndims(correct) == 2);
    correct = correct(:); estimate = estimate(:); mask = mask(:);

    if sum(estimate.^2 .* mask) > 1e-5
        alpha = sum(correct .* estimate .* mask) / sum(estimate.^2 .* mask);
    else
        alpha = 0.;
    end
    error_val = sum(mask .* (correct - alpha*estimate).^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the following functions are error metrics defined in:
%
%   Correlation-based Intrinsic Image Extraction from a Single Image.
%   Xiaoyue Jiang et. al. eccv 2010 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ascore] = ascore_image(true_shading, true_refl, estimate_shading, estimate_refl, mask, window_size)
    if nargin < 6
        window_size = 20;
    end
    ascore = 0.5 * alocal_error(true_shading, estimate_shading, mask, window_size, floor(window_size/2)) + ...
             0.5 * alocal_error(true_refl, estimate_refl, mask, window_size, floor(window_size/2));
end

function [error_val] = alocal_error(correct, estimate, mask, window_size, window_shift)
    M = size(correct,1);
    N = size(correct,2);

    assq = 0;
    total = 0;
    for i = 1:window_shift:(M - window_size + 1)
        for j = 1:window_shift:(N - window_size + 1)
            correct_curr = correct(i:i+window_size-1, j:j+window_size-1);
            estimate_curr = estimate(i:i+window_size-1, j:j+window_size-1);
            mask_curr = mask(i:i+window_size-1, j:j+window_size-1);
            assq = assq + assq_error(correct_curr, estimate_curr, mask_curr);
            total = total + sum(mask_curr(:) .* correct_curr(:) .^ 2);
        end
    end

    assert(~isnan(assq/total));
    error_val = assq / total;
end

function [error_val] = assq_error(correct, estimate, mask)
    correct = correct - mean(correct(:));
    estimate = estimate - mean(estimate(:));
    error_val = ssq_error(correct, estimate, mask);
end