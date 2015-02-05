function [r_new] = updateReflCF(r_init, normal, ims, mask, lit_mask, tr, ts)
% closed-form solution to update reflectance
% note this optimization is in log domain

if nargin < 6
    tr = 0.008;
    ts = 0.003;
end

% use 3x3 neighbor
wndSize = 9;
findNeighbor = @(p,h) [ ...
    p-h-1 p-h p-h+1 ...
    p-1 p+1 ...
    p+h-1 p+h p+h+1 ...
]';

% % use 5x5 local window
% wndSize = 25;
% findNeighbor = @(p,h) [...
%     p-2*h-2 p-2*h-1 p-2*h p-2*h+1 p-2*h+2 ...
%     p-h-2 p-h-1 p-h p-h+1 p-h+2 ...
%     p-2 p-1 p+1 p+2 ...
%     p+h-2 p+h-1 p+h p+h+1 p+h+2 ...
%     p+2*h-2 p+2*h-1 p+2*h p+2*h+1 p+2*h+2 ...
% ]';

mask = logical(mask);
[h w] = size(mask);
F = size(ims,4);
P = sum(mask(:));
logI = log( squeeze(mean(ims,3))+eps );
logI = reshape(logI, [h*w F]);
lit_mask = reshape(lit_mask, [h*w F]);
normal = reshape(normal, [h*w 3]);

chroma = sum(ims, 4);
chroma = chroma ./ repmat( max(eps,mean(chroma,3)), [1 1 3]);
c = reshape(chroma, [h*w 3]);
c = c/max(c(:));

ind = find(mask);
inv_ind = zeros(h,w);
inv_ind(mask) = 1:P;

nz_max = P*wndSize;
A_val = zeros(nz_max,1);
A_row = zeros(size(A_val));
A_col = zeros(size(A_val));
b = zeros(P,1);

nz = 0;
for i = 1:P
    p = ind(i);
    neigInd = findNeighbor(p,h);
    neigInd = neigInd(mask(neigInd));
    neigNum = length(neigInd);
    
    wr = 100*ones(neigNum,1);
    ws = ones(neigNum,F);
    
    for n = 1:neigNum
        q = neigInd(n);
        if sum( (c(p,:)-c(q,:)).^2 ) > tr*tr            % wr by chroma
            wr(n) = 0;
        end
        
        ws(n, xor(lit_mask(p,:), lit_mask(q,:)) ) = 0;  % ws by shadow mask
        if 1 - sum(normal(p,:).*normal(q,:)) > ts       % ws by normal
            ws(n,:) = 0.01;
            wr(n) = wr(n)+0.01;
        end
    end
    
    A_val(nz+1:nz+neigNum) = -(F*wr + sum(ws,2));
    A_row(nz+1:nz+neigNum) = i*ones(neigNum,1);
    A_col(nz+1:nz+neigNum) = inv_ind(neigInd);
    nz = nz+neigNum+1;
    A_val(nz) = sum(F*wr + sum(ws,2));
    A_row(nz) = i;
    A_col(nz) = i;
    
    tmp = ws .* ( repmat(logI(p,:),[neigNum 1]) - logI(neigInd,:) );
    b(i) = sum(tmp(:));
end

A_val = A_val(1:nz);
A_row = A_row(1:nz);
A_col = A_col(1:nz);

A = sparse(A_row, A_col, A_val, P, P, nz);

warning off;
r = A\b;
warning on;

r = exp(r);
r = r/max(r(:));

r_new = zeros(size(mask));
r_new(mask) = r;

end