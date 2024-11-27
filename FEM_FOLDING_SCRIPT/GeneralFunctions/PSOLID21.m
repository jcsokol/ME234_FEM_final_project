% Anisotropic stretch driven growth: Fg = I + (theta-1)*e1*e1, theta_dot =
% Gs*(lambdaE-lambdaCrit
function psolid = PSOLID21(material,var)
    % Initialize variables
    psolid.material = material;
    psolid.Gs = var.Gs;
    psolid.stretchCrit = var.lambdaCrit;
    
    % Functions
    psolid.Piola2Stiffness = @piola2Stiffness;
    function [Svoigt,Cvoigt, Fg] = piola2Stiffness(F,Fg0,dt)
        error('Piola2Stiffness not implemented for PSOLID21');
    end

    psolid.Piola1Stiffness = @piola1Stiffness;
    function [P, A, Fg] = piola1Stiffness(F,Fg0,dt,ndim)
        th0 = Fg0(1,1);
        
        stretch = norm(F(:,1));
        if(dt>0)
            [th,K] = updateFg(F,th0,dt);
        else
            th = th0;
        end
        
        Fg = eye(ndim); Fginv = Fg;
        Fg(1,1) = th; Fginv(1,1) = 1/th;
        Fe = F * Fginv;
        
        %{
        F
        stretch
        Fg
        (stretch/th-psolid.stretchCrit)*psolid.Gs
        %}
        e1 = zeros(ndim,1); e1(1,1)=1;
        thRat = (1-th)/th;
        
        % Compute stress, and stress derivative (A)
        [Pe,Ae] = material.material.Piola1Stiffness(Fe,ndim);
        
        % Compute P_ij = Pe_ij + thRat*Pe_i1*e1_j
        P = Pe + thRat*Pe(:,1)*e1';
        
        if(dt>0)
            % Compute A_ijkl = dP_ij/dF_kl |Fg (ONLY FIRST CONTRIBUTION)
            A = zeros(ndim,ndim,ndim,ndim);
            for i=1:ndim; for j=1:ndim for k=1:ndim; for l=1:ndim
                A(i,j,k,l)  = Ae(i,j,k,l) + thRat * (Ae(i,j,k,1)*e1(l) + Ae(i,1,k,l)*e1(j)) ...
                                          + thRat^2 * Ae(i,1,k,1)*e1(j)*e1(l);
            end;end;end;end
            
            % Compute dP_ij/dth = -P_ij/th - Ae_ijrl*Fe_rl / th^2
            dPdTH = zeros(ndim,ndim);
            for i=1:ndim; for j=1:ndim; 
                dPdTH(i,j) = -P(i,1)*e1(j)/th;    
                for r=1:ndim; 
                    dPdTH(i,j) = dPdTH(i,j) - Fe(r,1)*( Ae(i,j,r,1)/th + ...
                                    Ae(i,1,r,1)*e1(j)*(1-th)/th^2);
                end; 
            end; end;
    
            % Compute dTHdF_pq = ...*F_p1*f_q
            dTHdF = psolid.Gs*dt/K/th/stretch*F(:,1)*e1';
            
            % Compute A = A + dPdTH_ij*dTHdF_kl
            for i=1:ndim;for j=1:ndim;for k=1:ndim; for l=1:ndim	
                A(i,j,k,l) = A(i,j,k,l) + dPdTH(i,j)*dTHdF(k,l);
            end, end, end, end
        end
        
    end

    psolid.Diffusivity = @diffusivity;
    function [q,D,D3]  = diffusivity(F,gradRho,ndim)
        q = material.material.D*gradRho;
        D = material.material.D;                % = dq/dgradRho
        D3 = zeros(ndim,ndim,ndim);              % = dq/dF
    end

    function [th,K] = updateFg(F,th0,dt)
        th = th0;
        stretch = norm(F(:,1));
        
        it = 0;
        epsilon = 1;
        while(epsilon>10^(-12) && it < 50)
            it = it+1;
           
            R = th-th0-psolid.Gs*(stretch/th-psolid.stretchCrit)*dt;
            K = 1+psolid.Gs*stretch/th^2*dt;
            epsilon = abs(R/K);
            
            th = th - R/K;
        end
    end
end


