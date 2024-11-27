% Anisotropic prescribed area growth: Fg = theta*I + (1-theta)*e1*e1
function psolid = PSOLID12(material,var)
    % Initialize variables
    psolid.material = material;
    psolid.th_rate = var.th_rate;
    
    % Functions
    psolid.Piola2Stiffness = @piola2Stiffness;
    function [Svoigt,Cvoigt, Fg] = piola2Stiffness(F,Fg0,dt)
        error('Piola2Stiffness not implemented for PSOLID12');
    end


    psolid.Piola1Stiffness = @piola1Stiffness;
    function [P, A, Fg] = piola1Stiffness(F,Fg0,dt,ndim)
        theta = Fg0(end,end)+ dt*psolid.th_rate;
        
        Fg = theta*eye(ndim); Fginv = 1/theta*eye(ndim);
        Fg(1,1) = 1; Fginv(1,1) = 1;
        Fe = F * Fginv;
        
        [Pe,Ae] = material.material.Piola1Stiffness(Fe,ndim);
        
        e1 = zeros(1,ndim); e1(1,1)=1;
        th_1 = theta-1;
        
        P = Pe/theta; A = zeros(ndim,ndim,ndim,ndim);
        for i=1:ndim; for j=1:ndim
            P(i,j) = P(i,j) + th_1/theta * Pe(i,1)*e1(j);
        end;end
    
        for i=1:ndim; for j=1:ndim for k=1:ndim; for l=1:ndim
            A(i,j,k,l) = 1/theta^2*( Ae(i,j,k,l) + th_1 * ...
                                    (Ae(i,j,k,1)*e1(l) + Ae(i,1,k,l)*e1(j)) ...
                                      + th_1^2 * Ae(i,1,k,1)*e1(j)*e1(l));
        end;end;end;end
    end

    psolid.Diffusivity = @diffusivity;
    function [Q,D]  = diffusivity(F,GradRho,ndim)
        D = material.material.D;
        Q = -D*GradRho;
    end

end