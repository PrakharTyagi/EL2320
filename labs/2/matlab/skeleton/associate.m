% function [outlier,Psi] = associate(S_bar,z,W,Lambda_psi,Q)
%           S_bar(t)            4XM
%           z(t)                2Xn
%           W                   2XN
%           Lambda_psi          1X1
%           Q                   2X2
% Outputs: 
%           outlier             1Xn
%           Psi(t)              1XnXM
function [outlier,Psi] = associate(S_bar,z,W,Lambda_psi,Q)
    if isempty(z)
        outlier = [0];
        np = size(S_bar,2);
        Psi = repmat(1/np, 1,np);
    else
        % sizes for all the relevant bits
        nlandmarks = size(W,2);
        nobservations = size(z,2);
        nparticles = size(S_bar,2);
        % precalculate the normalisation constant
        normalisation = 1/(2*pi*sqrt(det(Q)));
        Z=zeros(2,size(S_bar,2),nlandmarks);
        Psi = zeros(1,max(1,nobservations),nparticles);
        % precompute values for k so that they do not have to be recomputed
        for k=1:nlandmarks
            Z(:,:,k)=observation_model(S_bar,W,k);
        end
        for i=1:nobservations
            for k=1:nlandmarks
                nu = repmat(z(:,i),1,nparticles) - Z(:,:,k);
                nu(2,:)=mod(nu(2,:)+pi,2*pi)-pi;
                % We only care about the diagonal part of the computed matrix,
                % which contains the likelihoods of the individual particles.
                % If we were calculating on a per-particle basis this would
                % give us a scalar - diag allows us to extract that value.
                psi(k,:) = diag(normalisation*exp(-0.5*nu'*inv(Q)*nu));
            end
            % We are only interested in the landmark which maximises the
            % likelihood of the current observation for each particle
            Psi(1,i,:) = max(psi);
        end
        % Check whether the mean value of each row of psi is an outlier or not
        outlier = mean(reshape(Psi, nobservations, nparticles),2) <= Lambda_psi;
    end
end