function [J, Jinv, dNr,N] = Jacobian_H2(r,s,r0)
% INPUTS:
% r,s,t: coordinates of integration point in element coordinates
% r0   : nodal position vector (global coordinates) (!!! there is no relation between r0 and r !!!)
% 
% OUTPUTS:
% jacobian:  Jacobian matrix with d[x; y;z] = jacobian d[r;s;t]
% Jinv:      inverse of jacobian
% dNr :      derivative of shape functions wrt r,s,t
% N   :      shape function values

% Write derivatives of shape functions wrt r,s,t
% Table 3.5.1 in The finite element method by T.J.R. Hughes:
gp = [-1 -1; 1 -1; 1 1; -1 1];

dNdr = gp(:,1)         .*  (1+gp(:,2)*s) /4;
dNds = (1+gp(:,1)*r)   .*  gp(:,2)       /4;

dNr = [dNdr dNds];

N = (1+gp(:,1)*r).*(1+gp(:,2)*s)/4;

J = reshape(r0,2,4)*dNr;
Jinv = inv(J);
end