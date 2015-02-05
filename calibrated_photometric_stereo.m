% Computes the facet vector (surface normal times albedo) at each pixel given a set of images
% illuminated by a varying point light source and the light source directions / intensities.
%
%  function [B,S] = calibrated_photometric_stereo(I,S,mask)
%
% I     : NxM Image data; N is the number of pixels per image, M is the number of images.
% S     : 3xM Light source directions/intensities. The ith column corresponds to the light source
%         for the ith image. The direction of the vector specifies the light source direction and
%         the magnitude represents the light source intensity.
% mask  : Nx1 Mask image, specifying which pixels to operate on - for points not
%         in the mask, behavior is undefined
% B     : Nx3 Facet matrix (N=W*H). Each row corresponds to the surface normal times the albedo of the ith pixel.
%
% OR
%
%  function [B,S] = calibrated_photometric_stereo(ims,S,mask)
%
% where ims is a cell array of grayscale images and output B is a Nx3 matrix with N=H*W
%
% ============
% Neil Alldrin
%
function [B,S] = calibrated_photometric_stereo(I,S,mask)

% Note: mask isn't used for this version of PS, since the solution is purely
%       local and very little overhead is incurred by computing a solution at
%       every pixel.

if iscell(I)
	Ivec = cell2mat(cellfun(@(x)(x(:)),I(:)','UniformOutput',false));
	[B,S] = calibrated_photometric_stereo(Ivec, S);
	%B = reshape(B,[size(I{1}),3]);
else
	% I = B*S ==> B = I*pinv(S)
	B = I*pinv(S);
end
