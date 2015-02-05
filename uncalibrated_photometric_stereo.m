% function [B,S] = uncalibrated_photometric_stereo(Icell,mask,gbrinit)
%
% An implementation of the uncalibrated photometric stereo technique described in 
% "Shape and Albedo from Multple Images using Integrability" by Yuille and Snow.
%
% I     : NxM Image data; N is the number of pixels per image, M is the number of images.
%          OR
%         Mx1 Cell array of HxW grayscale images
% Isize : 1x2 Size of a single image.
% Imask : WxH Mask image to exclude pixels from being considered.
% B     : Nx3 Facet matrix (N=W*H). Each row corresponds to the surface normal times the albedo of the ith pixel.
% S     : 3xM Light source directions/intensities. The ith column corresponds to the light source
%         for the ith image. The direction of the vector specifies the light source direction and
%         the magnitude represents the light source intensity.
%
% ============
% Neil Alldrin
%
function [B,S] = uncalibrated_photometric_stereo(I,mask,gbrinit)

Isize = size(mask);
if ~exist('gbrinit','var') gbrinit = [1 1 1]; end;

if iscell(I)
	I = cell2mat(cellfun(@(x)(x(:)),I(:)','UniformOutput',false)); % NxM
end

I = im2double(I); % I is required to be of type double for later operations to work

Iinds = find(mask(:)>0); % indices of non-masked pixels
M = size(I,2);          % # of images / illuminations
N = size(I,1);  % # of non-masked pixels

% Transpose I so that it's consistent with the yuille and snow paper
I_ = I'; % non-masked
I = I_(:,Iinds); % masked

% We should have the following relation at the end
%  I_' = B*S;

% get the 4 most significant singular values
%  note: we are using the notation in the paper
%  note 2: need 4 singular values b/c of ambient term
[f D e] = svds(I,4);
e_ = I_' * pinv(D*f');

% f : Mx4
% e : Nx4
% D : 4x4 diagonal matrix

% We now need to restrict f and e by enforcing integrability

%  compute de/dx and de/dy
smooth_kernel = fspecial('gaussian',[25 25],2.5);
%smooth_kernel = fspecial('disk',3.0);
%smooth_kernel = 1;
for i = 1:4
  e_smooth = imfilter(reshape(e_(:,i).*mask(:),Isize),smooth_kernel,'same');
  [dedx_ dedy_] = gradient(e_smooth);
  dedy_ = -dedy_;
  dedx(:,i) = dedx_(Iinds);
  dedy(:,i) = dedy_(Iinds);
end

%  solve the following system of equations for x_ij and y_ij
%   sum_{j=2:4} sum_{i=1:j-1} ( a_ij*x_ij - b_ij*y_ij) = 0
%  where,
%   a_ij = e(i)*dedx(j)-e(j)*dedx(i)
%   c_ij = P_3i*P_2j - P_2i*P_3j
%   b_ij = e(i)*dedy(j)-e(j)*dedy(i)
%   d_ij = P_3i*P_1j - P_1i*P_3j

a12 = e(:,1).*dedx(:,2) - e(:,2).*dedx(:,1);
a13 = e(:,1).*dedx(:,3) - e(:,3).*dedx(:,1);
a23 = e(:,2).*dedx(:,3) - e(:,3).*dedx(:,2);
a14 = e(:,1).*dedx(:,4) - e(:,4).*dedx(:,1);
a24 = e(:,2).*dedx(:,4) - e(:,4).*dedx(:,2);
a34 = e(:,3).*dedx(:,4) - e(:,4).*dedx(:,3);

b12 = e(:,1).*dedy(:,2) - e(:,2).*dedy(:,1);
b13 = e(:,1).*dedy(:,3) - e(:,3).*dedy(:,1);
b23 = e(:,2).*dedy(:,3) - e(:,3).*dedy(:,2);
b14 = e(:,1).*dedy(:,4) - e(:,4).*dedy(:,1);
b24 = e(:,2).*dedy(:,4) - e(:,4).*dedy(:,2);
b34 = e(:,3).*dedy(:,4) - e(:,4).*dedy(:,3);

% Solve the system A*x = 0 (terms with index of 4 correspond to ambient term)
%   A = [a12 a13 a23 -b12 -b13 -b23 a14 a24 a34 -b14 -b24 -b34]
%   x = [c12 c13 c23  d12  d13  d23 c14 c24 c34  d14  d24  d34];
A = [a12 a13 a23 -b12 -b13 -b23 a14 a24 a34 -b14 -b24 -b34];
%[AU AS AV] = svds(A,12);
%x = AV(:,12);
[AU,AS,AV] = svd(A(:,1:6),0);
x = AV(:,size(AV,2));

%  construct the "co-factor matrix" (eqn 12) from the paper
%  P3inv = [-c23  d23 ???;
%            c13 -d13 ???;
%           -c12  d12 ???];
P3inv = [-x(3)  x(6)  gbrinit(1); ...
	  x(2) -x(5)  gbrinit(2); ...
	 -x(1)  x(4)  gbrinit(3)];

%  invert P3inv to get P3 (up to GBR)
P3 = inv(P3inv);

%  solve for Q3 using the relation P3'*Q3 = D3
%Q3 = (inv(P3')*D(1:3,1:3));
Q3 = (P3') \ D(1:3,1:3);

% %  compute the remaining portions of P4 / Q4
% w = sum(f,1)';
% %  set up linear system to solve p4 using the remaining cross-product terms
% %   A*p4 = b
% %   p4 = [P14 P24 P34]'
% %   b  = [c14 c24 c34 d14 d24 d34]'
% A = [      0 P3(3,1) -P3(2,1);
%            0 P3(3,2) -P3(2,2);
%            0 P3(3,3) -P3(2,3);
%      P3(3,1)       0 -P3(1,1);
%      P3(3,2)       0 -P3(1,2);
%      P3(3,3)       0 -P3(1,3)];
% b = x(7:12);
% p4 = pinv(A)*b;

% Compute the B and S matrices from P3 and Q3
B = e_(:,1:3)*P3';
S = Q3*f(:,1:3)';

% Adjust the sign so that nz is positive
if mean(B(mask(:),3))<0
  B = B*[1 0 0; 0 1 0; 0 0 -1];
  S = [1 0 0; 0 1 0; 0 0 -1]*S;
end

% Remove pixels in shadow / highlight by estimating a visibility matrix V

% I = B*S; [Nx4] = [Nx3]*[3x4]

% [I1; I2; I3; I4] = [B*S1; B*S2; B*S3; B*S4]
% [I1; I2; I3; I4] = [(V1.*B)*S1; (V2.*B)*S2; (V3.*B)*S3; (V4.*B)*S4]
% [I1; I2; I3; I4] = [(V1.*B); (V2.*B); (V3.*B); (V4.*B)] * S_
% [4*Nx1] = [4*Nx3]*[]

% What constitutes an outlier???
%  [1x1] = [1x3]*[3x1]; A single surface point with one light source
