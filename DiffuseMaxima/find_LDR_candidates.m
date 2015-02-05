function [LDR]= find_LDR_candidates(I,mask)
% Produces a first estimate of the spatial locations for the LDR maxima 
% Inputs: 
%   I       :   NxM Image Data; N is the number of pixels per image, M is the number of images.
%   MASK    :   A binary mask to segment the object;
% Outputs:
%   LDR     :   A MxN binary matrix ( where is M is the number of pixels of each
%   image and N is the number of images) having ones at the local maxima spatial locations and zeros elsewhere
%============
% Author: Thoma Papadhimitri
% http://www.cvg.unibe.ch/staff/thoma-papadhimitri
% November 2011
% All rights reserved
% Research use is allowed provided that the following work is cited:
% "A Closed-Form Solution to Uncalibrated Photometric Stereo via Diffuse
% Maxima" by Paolo Favaro and Thoma Papadhimitri, CVPR 2012.
%============


nrows=size(mask,1);ncols=size(mask,2);num_images=size(I,3);
% We find the initial LDR points and also perform a Gaussian smoothing:
LDR=zeros(nrows,ncols,num_images);
[x,y] = meshgrid(-11:11,-11:11);
kernel = exp(-(x.^2+y.^2)/2/2^2);
weights = conv2(ones(nrows,ncols),kernel,'same');
for k=1:num_images
    LDR(:,:,k)=(imregionalmax(conv2(I(:,:,k).*double(mask),kernel,'same')./weights)); 
end

% As described in the paper we consider a region around each LDR maximum
% with radius 1:
se = strel('disk',1);
for i=1:num_images
    LDR(:,:,i) = imdilate(LDR(:,:,i),se); 
end

% We do not consider LDR maxima at the same spatial location in at least
% two images:
indm=sum(LDR,3)>1;
LDR=reshape(LDR,nrows*ncols,num_images);
for k=1:num_images
    LDR(indm,k)=0;                           
end

%We discard low intensity LDR maxima:
thresholdil=zeros(num_images);
for k=1:num_images
    thresholdil(k)=.5*(max(max(I(:,:,k)))-min(min(I(:,:,k))));     
end
filter=zeros(nrows,ncols,num_images);
for k=1:num_images
    filter(:,:,k)=I(:,:,k)>=thresholdil(k);
end
LDR=reshape(LDR,nrows,ncols,num_images).*filter;   % The final selected LDR maxima