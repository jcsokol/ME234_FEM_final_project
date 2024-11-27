% Anisotropic prescribed growth: Fg = I + (theta-1)*e1*e1
function psolid = PSOLID11(material,var)
    % Initialize variables
    psolid.material = material;
    psolid.th_rate = var.th_rate;
    
    % Functions
    psolid.Piola2Stiffness = @piola2Stiffness;
    function [Svoigt,Cvoigt, Fg] = piola2Stiffness(F,Fg0,dt)
        error('Piola2Stiffness not implemented for PSOLID21');
    end


    psolid.Piola1Stiffness = @piola1Stiffness;
    function [P, A, Fg] = piola1Stiffness(F,Fg0,dt,ndim)
        theta = Fg0(1,1)+ dt*psolid.th_rate;
        
        Fg = eye(ndim); Fginv = Fg;
        Fg(1,1) = theta; Fginv(1,1) = 1/theta;
        Fe = F * Fginv;
        
        [Pe,Ae] = material.material.Piola1Stiffness(Fe,ndim);
        
        e1 = zeros(1,ndim); e1(1,1)=1;
        thRat = (1-theta)/theta;
        
        P = Pe; A = zeros(ndim,ndim,ndim,ndim);
        for i=1:ndim; for j=1:ndim
            P(i,j) = P(i,j) + thRat * Pe(i,1)*e1(j);
        end;end
    
        for i=1:ndim; for j=1:ndim for k=1:ndim; for l=1:ndim
            A(i,j,k,l)  = Ae(i,j,k,l) + thRat * (Ae(i,j,k,1)*e1(l) + Ae(i,1,k,l)*e1(j)) ...
                                      + thRat^2 * Ae(i,1,k,1)*e1(j)*e1(l);
        end;end;end;end
    end

    psolid.Diffusivity = @diffusivity;
    function [q,D,D3]  = diffusivity(F,gradRho,ndim)
        D = material.material.D;% = dq/dgradRho 
        q = D*gradRho;
        D3 = zeros(ndim,ndim,ndim);              % = dq/dF
    end

end