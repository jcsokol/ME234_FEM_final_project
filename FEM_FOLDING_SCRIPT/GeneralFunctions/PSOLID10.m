% Isotropic prescribed growth: Fg = theta*I
function psolid = PSOLID10(material,var)
    % Initialize variables
    psolid.material = material;
    psolid.th_rate = var.th_rate;
    
    % Functions
    psolid.Piola2Stiffness = @piola2Stiffness;
    function [Svoigt,Cvoigt, Fg] = piola2Stiffness(F,Fg0,dt)
        theta = Fg0(1,1)+ dt*psolid.th_rate;
        Fe = F/theta;
        Fg = theta*eye(3);
        
        [Sevoigt,Cevoigt] = material.material.Piola2Stiffness(Fe);
        
        Svoigt = Sevoigt/(theta^2);
        Cvoigt = Cevoigt/(theta^4);
    end


    psolid.Piola1Stiffness = @piola1Stiffness;
    function [P, A, Fg] = piola1Stiffness(F,Fg0,dt,ndim)
        theta = Fg0(1,1)+ dt*psolid.th_rate;
        Fe = F/theta;
        Fg = theta*eye(ndim);
        
        [Pe,Ae] = material.material.Piola1Stiffness(Fe,ndim);
        
        P = Pe/theta;
        A = Ae/theta^2;
    end
    
    psolid.Diffusivity = @diffusivity;
    function [Q,D]  = diffusivity(F,GradRho,ndim)
        D = material.material.D;
        Q = -D*GradRho;
    end

end