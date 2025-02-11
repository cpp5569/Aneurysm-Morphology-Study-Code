function [cal_mat]=vortex_computation_3D(vortex_method,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
% function [cal_mat,grid_y,grid_x]=vortex_computation(U,V,vortex_method,gradient_method,mag,vector_skip,vector_spacing)
% Compute possible vortex perimeters
% Analyzes output from POD reconstructed fields or vector fields in general
% Input must be 2D vector field
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code:
%   1. Computes isocontours of possible vortex locations
%   2. Gives results in physical units in user-defined grid
% 
% Required inputs:
%   1. 3D vector field (U,V,W)
%   1. Method to compute vortex field (vortex_method, Qcrit, Dcrit, or Lamda)
%   2. 3D velocity gradient fields:
%   (dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
%
% Adapted by Melissa Brindise (11.05.2018) from:
% Chris Weiland (orig by Mike Brady)
% 10.02.2007
% v1.0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Define grid operating on
[y_L,x_L,z_L]=size(dudx);

if strcmpi(vortex_method,'Lamda')
    vortex_method = 'Lambda';
end

% 2. Initialize vortex matrix
cal_mat=zeros(size(dudx));

% 3. Compute based on user-defined algorithmn
for k=3:1:z_L-2
    for j=3:1:x_L-2
        for i=3:1:y_L-2
            velgradient=[dudx(i,j,k) dudy(i,j,k) dudz(i,j,k);...
                dvdx(i,j,k) dvdy(i,j,k) dvdz(i,j,k); ...
                dwdx(i,j,k) dwdy(i,j,k) dwdz(i,j,k)]; % velocity gradient tensor
            % Compute the symmetric and anti-symmetric parts of the
            % velocity gradient tensor
            S=0.5*(velgradient+transpose(velgradient)); % this is symmetric
            Omega=0.5*(velgradient-transpose(velgradient)); % this is anti-symmetric
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 4. Choose computation method for vortex identification & do
            % computation     
            switch lower(vortex_method)
                case 'lambdaci'
                    %velgradient=[dudx(i,j) dudy(i,j); dvdx(i,j) dvdy(i,j)];
                    try
                        E = eig(velgradient);
                    catch
                        keyboard
                    end
                    cal_mat(i,j,k) = max(imag(E));

                case 'qcrit'
                    % Compares magnitude of straining & rotational motions
                    Omegan=norm(Omega,'fro'); % why normalizing these?
                    Sn=norm(S,'fro');
                    Q=0.5*(Omegan.^2-Sn.^2);
                    if Q>0
                        cal_mat(i,j,k)=Q;
                    end

                case 'dcrit'
                    % Cantwell's method
                    Omegan=norm(Omega,'fro');
                    Sn=norm(S,'fro');
                    Q=0.5*(Omegan.^2-Sn.^2);
                    Delta=(Q/3)^3+(det(velgradient)/2)^2;
                    if Delta>0
                        cal_mat(i,j,k)=Delta;
                    end

                case 'lambda'
                    % Lambda_2 method
                    try
                        T=eig(S^2+Omega^2); % principal axes?
                        Lambda = sort(T,'descend');
                        if Lambda(2) < 0
                            cal_mat(i,j,k)=Lambda(2);
                        else
                            cal_mat(i,j,k)=0;
                        end
                    catch
                        cal_mat(i,j,k) = 0;
                    end
                    
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end


end % End of function