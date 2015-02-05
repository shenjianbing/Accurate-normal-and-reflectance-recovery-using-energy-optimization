% 简单的阈值计算
% function [t] = getThreshold(img, mask)
% 
% m_img = img(mask);
% bins = 100;
% [freq val] = hist(m_img, bins);
% 
% d_p = freq(2:bins-1) - freq(1:bins-2);
% d_n = freq(2:bins-1) - freq(3:bins);
% flag = (d_p < 0) & (d_n < 0);
% 
% valley = find(flag);
% if isempty(valley)
%     t = 0.1*max(img(:));
% else
%     t = val(valley(1));
% end

% 使用高斯混合模型计算阈值
function [t] = getThreshold(img, mask)

m_img = img(mask);
x = m_img(:)';

nModel = 2;
k_centers = [0; 1];
[Priors, Mu, Sigma] = EM_init_kmeans(x, nModel, k_centers);
[Priors, Mu, Sigma] = EM(x, Priors, Mu, Sigma);

% 两个高斯分量的处理
if Mu(1) < Mu(2)
    t1 = min(Mu(1)+3*Sigma(1), Mu(2)-3*Sigma(2));
    t2 = max(Mu(1)+3*Sigma(1), Mu(2)-3*Sigma(2));
else
    t1 = min(Mu(2)+3*Sigma(2), Mu(1)-3*Sigma(1));
    t2 = max(Mu(2)+3*Sigma(2), Mu(1)-3*Sigma(1));
end

t = (t1+t2)/2;

end

