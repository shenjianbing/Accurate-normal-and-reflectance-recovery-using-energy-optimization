% Converts an Nx3 facet matrix into unit normal and albedo images
% 
%  function [nx,ny,nz,a] = B2normals(B,sizeI)
%
% ============
% Neil Alldrin
%
function [nx,ny,nz,a] = B2normals(B,sizeI)

dim=size(B,1);


B = reshape(B,dim,3);

a = reshape(sqrt(sum(B.^2,2)),sizeI);
nx = reshape(B(:,1),sizeI)./(a+(a==0));
ny = reshape(B(:,2),sizeI)./(a+(a==0));
nz = reshape(B(:,3),sizeI)./(a+(a==0));
