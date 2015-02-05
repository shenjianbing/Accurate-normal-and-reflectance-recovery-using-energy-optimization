function [mu,nu,lambda] = loc_search(Np,Lk)
% Produces the estimated GBR parameters by finding the median of all the
% intersections in the GBR parameter space as described in the paper
%
% Inputs 
%   Np              :   The pseudo-normals matrix Nx3 (N is the number of pixels for each image).
%   Lk              :   The pseudo-lights matrix 3xM (M is the number of light sources);
% Outputs    
%   mu, nu, la      :   The 3 GBR paramters of the GBR matrix [1 0 0; 0 1 0; mu nu la]

%============
% Author: Paolo Favaro
% November 2011
% All rights reserved
% Research use is allowed provided that the following work is cited:
% "A Closed-Form Solution to Uncalibrated Photometric Stereo via Diffuse
% Maxima" by Paolo Favaro and Thoma Papadhimitri, CVPR 2012.
%============


tol = 0;
% compute maxima for each segment
segmentCost = zeros(1,size(Np,2));
for ii=1:size(Np,2)
    a11 = -Np(1,ii)./Np(3,ii);
    a12 = -Np(2,ii)./Np(3,ii);
    a21 = (-Lk(2,ii)^2*Np(1,ii)+Lk(1,ii)*Lk(2,ii)*Np(2,ii)+Lk(1,ii)*Lk(3,ii)*Np(3,ii))./( ...
        Np(3,ii).*(Lk(1,ii).^2+Lk(2,ii).^2) );
    a22 = ( Lk(1,ii)*Lk(2,ii)*Np(1,ii)-Lk(1,ii)^2*Np(2,ii)+Lk(2,ii)*Lk(3,ii)*Np(3,ii))./( ...
        Np(3,ii).*(Lk(1,ii).^2+Lk(2,ii).^2) );
    
    mu = (a11+a21)/2;
    nu = (a12+a22)/2;
    la = (...
        ((Np(1,ii)+mu*Np(3,ii)).^2+(Np(2,ii)+nu*Np(3,ii)).^2).*(mu*Lk(1,ii)+nu*Lk(2,ii)-Lk(3,ii)).^2 ...
        ./( Np(3,ii).^2.*(Lk(1,ii).^2+Lk(2,ii).^2) )...
        ).^(1/4);
    
    Ghat = [1 0 0; 0 1 0; mu nu la];
    Gihat = [1 0 0; 0 1 0; -mu/la -nu/la 1/la];
    
    segmentCost(ii) = (Np(:,ii)'*Lk(:,ii))/sqrt(sum((Ghat'*Np(:,ii)).^2))./sqrt(sum((Gihat*Lk(:,ii)).^2));
end



% we consider every LDR point
listIntersections = zeros(3,size(Np,2));
count = 0;
for ii=1:size(Np,2)
    a11 = -Np(1,ii)./Np(3,ii);
    a12 = -Np(2,ii)./Np(3,ii);
    a21 = (-Lk(2,ii)^2*Np(1,ii)+Lk(1,ii)*Lk(2,ii)*Np(2,ii)+Lk(1,ii)*Lk(3,ii)*Np(3,ii))./( ...
        Np(3,ii).*(Lk(1,ii).^2+Lk(2,ii).^2) );
    a22 = ( Lk(1,ii)*Lk(2,ii)*Np(1,ii)-Lk(1,ii)^2*Np(2,ii)+Lk(2,ii)*Lk(3,ii)*Np(3,ii))./( ...
        Np(3,ii).*(Lk(1,ii).^2+Lk(2,ii).^2) );
    
    jj=ii+1:size(Np,2);
    % and find intersections with every other LDR point
    b11 = -Np(1,jj)./Np(3,jj);
    b12 = -Np(2,jj)./Np(3,jj);
    b21 = (-Lk(2,jj).^2.*Np(1,jj)+Lk(1,jj).*Lk(2,jj).*Np(2,jj)+Lk(1,jj).*Lk(3,jj).*Np(3,jj))./( ...
        Np(3,jj).*(Lk(1,jj).^2+Lk(2,jj).^2) );
    b22 = ( Lk(1,jj).*Lk(2,jj).*Np(1,jj)-Lk(1,jj).^2.*Np(2,jj)+Lk(2,jj).*Lk(3,jj).*Np(3,jj))./( ...
        Np(3,jj).*(Lk(1,jj).^2+Lk(2,jj).^2) );
    
    x = (b12-a12-b11.*(b22-b12)./(b21-b11)+a11.*(a22-a12)./(a21-a11))./((a22-a12)./(a21-a11)-(b22-b12)./(b21-b11));
    y = (x-a11)./(a21-a11).*(a22-a12)+a12;
    z = (  ((Np(1,ii)+x*Np(3,ii)).^2+(Np(2,ii)+y*Np(3,ii)).^2).*(x*Lk(1,ii)+y*Lk(2,ii)-Lk(3,ii)).^2 ...
        ./( Np(3,ii).^2.*(Lk(1,ii).^2+Lk(2,ii).^2) )...
        ).^(1/4);
    
    % handle division by zero
    discard = (abs(a21-a11)<1e-6) | (abs(b21-b11)<1e-6);
    discard2 = abs(((a22-a12)./(discard+(a21-a11))-(b22-b12)./(discard+(b21-b11)))<1e-6);
    % intersection is out of the segment limits
    discard3 = (x<(1+tol)*max(min(a11,a21),min(b11,b21)))|(x>(1-tol)*min(max(a21,a11),max(b21,b11)))|(y<(1+tol)*max(min(a12,a22),min(b12,b22)))|(y>(1-tol)*min(max(a12,a22),max(b12,b22)));
    keep = (1-discard)&(1-discard2)&(1-discard3);
    
    if sum(keep)>0
        listIntersections(:,count+[1:sum(keep)]) = [x(keep);y(keep);z(keep)];
        count = count+sum(keep);
    end
    
end
% We find the median of all intersections
mu = median(listIntersections(1,:));
nu = median(listIntersections(2,:));
lambda = median(listIntersections(3,:));


return
