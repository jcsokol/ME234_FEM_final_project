function [F, BL, BNL, J, P,dNx,N] = DeformationStrain_H2(r0, r, s, v_i)

% INPUTS:
%
% r0:       element nodal position vector
% (r,s,t):  integration point position in element coordinates
% v_i:      displacement vector for this element
% 
% OUTPUTS:
% F     :   Deformation gradient
% BL    :   BL matrix
% BNL   :   BNL matrix such that dstrain/dv = (BL+BNL), strain =(BL+0.5BNL)*v
% J     :   Jacobian matrix
% P     :   Struct with P1,P2,P3,P23,P13,P12
% dNx   :   Matrix that contains dN(i)/dx(j) in i,j entry
% N     :   Vector with N(i) in entry i

        
    % Compute jacobian, Jinv, Ndr, Nds for this integration point
    [J, Jinv, dNr, N] = Jacobian_H2(r,s,r0);
    
    % Calculate derivatives of N wrt global coordinates
    dNx = dNr * Jinv;
    
    % Compute deformation gradient
    r_c = r0+v_i;
    F = reshape(r_c,2,4)*dNx;
    
    Ndx = zeros(2,8); Ndy = Ndx;
    for k=1:4
        Ndx(:,(k-1)*2+1:2*k) = dNx(k,1)*eye(2);
        Ndy(:,(k-1)*2+1:2*k) = dNx(k,2)*eye(2);
    end
    
    % Compute P matrices
    P1 = Ndx'*Ndx;      P2 = Ndy'*Ndy;     
    P12 = Ndx'*Ndy+Ndy'*Ndx;
    
    % Store P matrices in a struct
    P = struct('P1',P1,'P2',P2,'P12',P12);
    
    % Compute BL and BNL
    BL =  [r0.'*P1; r0.'*P2; r0.'*P12];
    BNL = [v_i.'*P1; v_i.'*P2; v_i.'*P12];
    
    %NOTE: A more compact way to do this (for C++):
    
    x0 = r0(1:2:end);
    y0 = r0(2:2:end);
    
    vx_i = v_i(1:2:end);
    vy_i = v_i(2:2:end);
    
    Ndx = dNx(:,1)'; Ndy=dNx(:,2)';
    P1 = Ndx'*Ndx;      P2 = Ndy'*Ndy;
    P12 = Ndx'*Ndy+Ndy'*Ndx;
    %{
    BL = zeros(6,24);
    BNL = zeros(6,24);
    
    BL(:,1:3:end) =  [x0.'*P1; x0.'*P2; x0.'*P3; x0.'*P23 ; x0.'*P13 ; x0.'*P12];
    BL(:,2:3:end) =  [y0.'*P1; y0.'*P2; y0.'*P3; y0.'*P23 ; y0.'*P13 ; y0.'*P12];
    BL(:,3:3:end) =  [z0.'*P1; z0.'*P2; z0.'*P3; z0.'*P23 ; z0.'*P13 ; z0.'*P12];
    
    BNL(:,1:3:end) =  [vx_i.'*P1; vx_i.'*P2; vx_i.'*P3; vx_i.'*P23 ; vx_i.'*P13 ; vx_i.'*P12];
    BNL(:,2:3:end) =  [vy_i.'*P1; vy_i.'*P2; vy_i.'*P3; vy_i.'*P23 ; vy_i.'*P13 ; vy_i.'*P12];
    BNL(:,3:3:end) =  [vz_i.'*P1; vz_i.'*P2; vz_i.'*P3; vz_i.'*P23 ; vz_i.'*P13 ; vz_i.'*P12]; 
    
    xC = x0+vx_i;
    yC = y0+vy_i;
    zC = z0+vz_i;
    
    B=zeros(6,24);
    B(:,1:3:end) =  [xC.'*P1; xC.'*P2; xC.'*P3; xC.'*P23 ; xC.'*P13 ; xC.'*P12];
    B(:,2:3:end) =  [yC.'*P1; yC.'*P2; yC.'*P3; yC.'*P23 ; yC.'*P13 ; yC.'*P12];
    B(:,3:3:end) =  [zC.'*P1; zC.'*P2; zC.'*P3; zC.'*P23 ; zC.'*P13 ; zC.'*P12];
    
    errorr = max(max(abs(B-BL-BNL)))
    
    
    %max(max(abs(BL-BL0)))
    %}
end