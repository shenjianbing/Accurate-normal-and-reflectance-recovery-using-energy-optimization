% Takes in a vector x = [mu, nu, lambda] and returns the entropy.
%
% ============
% Neil Alldrin
%
function e = cost_entropy(x,uB,binWidth,showHist)

G = [1 0 0; 0 1 0; x(1) x(2) x(3)];

if ~exist('binWidth') binWidth = range(sqrt(sum(uB.^2,2)))/128; end;
if ~exist('showHist') showHist = false; end;

rhoG = sqrt(sum((uB*G).^2,2));
rhoG = rhoG./mean(rhoG);

%[p bins] = hist(rhoG,256);
[p bins] = hist(rhoG,min(rhoG):binWidth:max(rhoG));
p = p/length(rhoG); % probability density
ind = find(p~=0); 
plogp = p(ind).*log2(p(ind));
e = -sum(plogp);% + log2(bins(2)-bins(1));

if (showHist) bar(bins,p); end;
