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

    % sizes for all the relevant bits
    nlandmarks = size(W,2);
    nobservations = size(z,2);
    nparticles = size(S_bar,2);
    % precalculate the normalisation constant
    normalisation = 1/(2*pi*sqrt(det(Q)));
    Z=zeros(2,size(S_bar,2),nlandmarks);
    % precompute values for k so that they do not have to be recomputed
    for k=1:nlandmarks
       Z(:,:,k)=observation_model(S_bar,W,k);
    end
    for i=1:nobservations
        for k=1:nlandmarks
            nu = repmat(z(:,i),1,nparticles) - Z(:,:,k);
            nu(2,:)=mod(nu(2,:)+pi,2*pi)-pi;
            % We only care about the diagonal part of the computed matrix,
            % which contains the likelihoods of             
            psi(k,:) = diag(normalisation*exp(-0.5*nu'*inv(Q)*nu));
        end
        Psi(1,i,:) = max(psi);
        d(i) = mean(Psi(1,i,:));
    end
    
    outlier = d <= Lambda_psi;
end