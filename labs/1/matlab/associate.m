% function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,M,Lambda_m,Q)
% This function should perform the maximum likelihood association and outlier detection.
% Note that the bearing error lies in the interval [-pi,pi)
%           mu_bar(t)           3X1
%           sigma_bar(t)        3X3
%           R                   3X3
%           Q                   2X2
%           z_i(t)              2X1
%           M                   2XN
%           Lambda_m            1X1
% Outputs: 
%           c(t)                1X1
%           outlier             1X1
%           nu^i(t)             2XN
%           S^i(t)              2X2XN
%           H^i(t)              2X3XN
function [c,outlier, nu, S, H] = associate(mu_bar,sigma_bar,z_i,M,Lambda_m,Q)

% go over all landmarks in m
S = [];
H = [];
psi = [];
nu = [];
Z = [];
for k = 1:size(M,2)
    % measurement for this landmark
    Z(:,:,k) = observation_model(mu_bar,M,k);
    H(:,:,k) = jacobian_observation_model(mu_bar,M,k,Z(:,:,k),1);
    S(:,:,k) = H(:,:,k)*sigma_bar*H(:,:,k)' + Q;
    nu(:,:,k) = z_i - Z(:,:,k);

    psi(k) = det(2*pi*S(:,:,k))^(-1/2)*exp((-1/2)*nu(:,:,k)'*inv(S(:,:,k))*nu(:,:,k));
end
[cmag c] = max(psi); % c is the arg, cmag is the value of that maximal index
if cmag > Lambda_m
    outlier = 1;
else
    outlier = 0;
end
%nu = nu(:,:,c);
%S = S(:,:,c);
%H = H(:,:,c);

end