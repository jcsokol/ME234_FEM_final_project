% Isotropic strain driven growth
function psolid = PSOLID20(material,var)
    % Initialize variables
    psolid.material = material;
    psolid.Gs = var.Gs;
    psolid.Jcrit = var.Jcrit;
    
    % Functions
    psolid.Piola2Stiffness = @piola2Stiffness;
    function [Svoigt,Cvoigt, Fg] = piola2Stiffness(F,Fg0,dt)
        error('Piola2Stiffness not implemented for PSOLID20');
    end

    psolid.Piola1Stiffness = @piola1Stiffness;
    function [P, A, Fg] = piola1Stiffness(F,Fg0,dt,ndim)
        th0 = Fg0(1,1);
        
        Finv = inv(F);
        
        if(dt>0)
            [th,K] = updateFg(F,th0,dt,ndim);
        else
            th = th0;
        end
        
        Fg=th*eye(ndim);
        Fe= F/th;
        
        % Compute stress, and stress derivative (A)
        [Pe,Ae] = material.material.Piola1Stiffness(Fe,ndim);
        P = Pe/th;
        
        if(dt>0)
            % Compute dP_ij/dth = -P_ij/th - Ae_ijrl*Fe_rl / th^2
            dPdTH = -P/th;
            for i=1:ndim; for j=1:ndim; for r=1:ndim; for l=1:ndim; 
                dPdTH(i,j) = dPdTH(i,j) - Ae(i,j,r,l)*Fe(r,l)/th^2;
            end; end; end;end;
    
            % Compute dTHdF
            dTHdF = psolid.Gs*dt/K*det(F)*Finv'/th^ndim;
            %dTHdF = psolid.Gs*dt/K*det(F)*Finv'/2/th^3;
           
            
            % Compute A = dPdF_ijkl + dPdTH_ij*dTHdF_kl
            A = zeros(ndim,ndim,ndim,ndim);
            for i=1:ndim;for j=1:ndim;for k=1:ndim; for l=1:ndim	
                A(i,j,k,l) = Ae(i,j,k,l)/th^2 + dPdTH(i,j)*dTHdF(k,l);
            end, end, end, end
            
            
            
        end
        
    end

    psolid.Diffusivity = @diffusivity;
    function [Q,D]  = diffusivity(F,GradRho,ndim)
        D = material.material.D;
        Q = -D*GradRho;
    end

    function [th,K] = updateFg(F,th0,dt,ndim)
        th = th0;
        J = det(F);
        
        it = 0;
        epsilon = 1;
        while(epsilon>10^(-12) && it < 50)
            it = it+1;
           
            R = th-th0-psolid.Gs*max(J/th^ndim-psolid.Jcrit,0)*dt;
            K = 1+ndim*psolid.Gs*J/th^(ndim+1)*dt;
            epsilon = abs(R/K);
            
            th = th - R/K;
        end
    end
end