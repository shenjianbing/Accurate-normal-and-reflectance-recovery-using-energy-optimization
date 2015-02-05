function [normal light albedo] = photometricStereo(I, mask, method, sign)
% method can be 'me', 'dm'
    
    if nargin < 4
        sign = 1;
    end
    % uncalibrated photometric stereo to get initial B and S
    gbrInit = [1 1 1];
    [uB uS] = uncalibrated_photometric_stereo(I, logical(mask), gbrInit);
    
    if strcmp(method, 'me')
        % resolve GBR ambiguity by minimize entropy (ME)
        display('UPS by minimizing entropy');
        [mu_me nu_me lambda_me e] = solve_gbr(uB(mask(:),:), uS);
        G_me = [1 0 0; 0 1 0; mu_me nu_me lambda_me*sign]; %
        B_me = uB * G_me;
        [nx_me ny_me nz_me a_me] = B2normals(B_me,size(mask));
        n_me = cat(2, nx_me(:), ny_me(:), nz_me(:));
        l_me = G_me\uS;
        %l_me = G_me\uS;

        normal = normr(n_me);
        light = normc(l_me);
        a_me(mask) = mat2gray(a_me(mask));
        albedo = a_me;
    else if strcmp(method, 'dm')
            % resolve GBR ambiguity via diffuse maxima (DM)
            display('UPS by diffuse maxima');
            nImages = size(I, 2);

            I = reshape(I, [size(mask) nImages]);
            ldr = find_LDR_candidates(I, mask);

            Si = [];
            for k = 1:nImages
                Si = [Si repmat(uS(:,k),[1,sum(sum(ldr(:,:,k)))])];
            end
            N_curr = repmat(uB,[nImages 1]);
            Ni = N_curr(ldr==1,:)';
            [mu_dm nu_dm lambda_dm] = loc_search(Ni,Si);
            G_dm = [1 0 0; 0 1 0; mu_dm nu_dm lambda_dm*sign]; %
            B_dm = uB * G_dm;
            [nx_dm ny_dm nz_dm a_dm] = B2normals(B_dm,size(mask));
            n_dm = cat(2, nx_dm(:), ny_dm(:), nz_dm(:));
            l_dm = G_dm\uS;

            normal = normr(n_dm);
            light = normc(l_dm);
            a_dm(mask) = mat2gray(a_dm(mask));
            albedo = a_dm;
        end
    end
        
end