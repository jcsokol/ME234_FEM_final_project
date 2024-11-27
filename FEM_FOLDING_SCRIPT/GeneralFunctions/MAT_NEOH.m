function mat_neoh = MAT_NEOH(matProp)
    % Initialize variables
    mat_neoh.E          = matProp.E;
    mat_neoh.nu         = matProp.nu;
    mat_neoh.alpha      = matProp.alpha;
    mat_neoh.tref       = matProp.tref;
    mat_neoh.rho        = matProp.rho;
    
    mat_neoh.mu     = mat_neoh.E / 2.0 / (1.0+mat_neoh.nu);
    mat_neoh.lambda = mat_neoh.E * mat_neoh.nu / (1.0+mat_neoh.nu) / ( 1.0-2.0*mat_neoh.nu ); 
    
    try
       mat_neoh.D = matProp.D;
    catch
       mat_neoh.D = eye(3);
    end
    
    
    
    % Functions
    mat_neoh.Piola2Stiffness = @piola2Stiffness;
    function [Svoigt,Cvoigt] = piola2Stiffness(Fe,ndim)
        
        mu = mat_neoh.mu;
        lambda = mat_neoh.lambda;
        
        J = det(Fe);
        Ce = Fe'*Fe;
        Cinv = inv(Ce);
        
        if(ndim==2)
            S = (lambda*log(J)-mu)*Cinv+mu*eye(2);
            Svoigt = [S(1,1);S(2,2);S(1,2)]; 

            vId1 = [1 2 1];
            vId2 = [1 2 2];

            Cvoigt = zeros(3);
            for row = 1:3
                for col = 1:3
                    i = vId1(row); j=vId2(row); k=vId1(col); l=vId2(col);

                    Cvoigt(row,col) =       lambda               * Cinv(i,j)*Cinv(k,l) ...
                                         +  (mu-lambda*log(J))   * (Cinv(i,k)*Cinv(j,l)+Cinv(i,l)*Cinv(k,j));  
                end
            end
        elseif(ndim==3)
            S = (lambda*log(J)-mu)*Cinv+mu*eye(3);
            Svoigt = [S(1,1);S(2,2);S(3,3);S(2,3);S(1,3);S(1,2)];

            vId1 = [1 2 3 2 1 1];
            vId2 = [1 2 3 3 3 2];

            Cvoigt = zeros(6);
            for row = 1:6
                for col = 1:6
                    i = vId1(row); j=vId2(row); k=vId1(col); l=vId2(col);

                    Cvoigt(row,col) =       lambda               * Cinv(i,j)*Cinv(k,l) ...
                                         +  (mu-lambda*log(J))   * (Cinv(i,k)*Cinv(j,l)+Cinv(i,l)*Cinv(k,j));  
                end
            end
        end
    end



    mat_neoh.Piola1Stiffness = @piola1Stiffness;
    function [Pe,Ae] = piola1Stiffness(Fe,ndim)
        mu = mat_neoh.mu;
        lambda = mat_neoh.lambda;
        
        delta = eye(ndim);
        Je = det(Fe);
        Feinv = inv(Fe);
        
        % 1st piola kirchhoff stress
        Pe = mu*Fe+(lambda*log(Je)-mu)*Feinv';
        
        % 4th order elastic stiffness Ae = dPe/dFe
        Ae = zeros(ndim,ndim,ndim,ndim);
        for i=1:ndim; for j=1:ndim; for k=1:ndim; for l=1:ndim
            Ae(i,j,k,l) =  lambda                  * Feinv(j,i) * Feinv(l,k) ...
                         - (lambda * log(Je) - mu) * Feinv(l,i) * Feinv(j,k) ...
                         +                  mu  * delta(i,k) * delta(j,l);     
        end; end; end; end
    end

    mat_neoh.Piola1StiffnessGrowth = @piola1StiffnessGrowth;
    function [P,dPdF,dPdFg] = piola1StiffnessGrowth(Fe,Fg,ndim)
        mu = mat_neoh.mu;
        lambda = mat_neoh.lambda;
        
        delta = eye(ndim);
        Je = det(Fe);
        Fginv = inv(Fg);
        Feinv = inv(Fe);
        
        % 1st piola kirchhoff stress
        Pe = mu*Fe+(lambda*log(Je)-mu)*Feinv';
        P = Pe*Fginv';
        
        % 4th order elastic stiffness Le = dPe/dFe
        Ae = zeros(ndim,ndim,ndim,ndim);
        for i=1:ndim; for j=1:ndim; for k=1:ndim; for l=1:ndim
            Ae(i,j,k,l) =  lambda                  * Feinv(j,i) * Feinv(l,k) ...
                         - (lambda * log(Je) - mu) * Feinv(l,i) * Feinv(j,k) ...
                         +                  mu  * delta(i,k) * delta(j,l);     
        end; end; end; end

        % 4th order dP_ij/dFg_kl
        dPdFg = zeros(ndim,ndim,ndim,ndim);
        for i=1:ndim; for j=1:ndim; for k=1:ndim; for l=1:ndim
            % 1st contribution
            dPdFg(i,j,k,l) = -P(i,l)*Fginv(j,k);
            
            % 2nd contribution, summation over p, q, r, s
            for p=1:ndim; for r=1:ndim; for s=1:ndim
                 dPdFg(i,j,k,l) = dPdFg(i,j,k,l) - Fginv(j,p)*Fginv(l,s) ...
                                                 * Ae(i,p,r,s)*Fe(r,k);
            end; end; end;
        end; end; end; end

        % 4th order dP_ij/dF_kl
        dPdF = zeros(ndim,ndim,ndim,ndim);
        for i=1:ndim; for j=1:ndim; for k=1:ndim; for l=1:ndim
            for p=1:ndim; for r=1:ndim
                dPdF(i,j,k,l) = dPdF(i,j,k,l) + Fginv(j,p)*Fginv(l,r)*Ae(i,p,k,r);
                end; end;
        end; end;end;end
    end

    mat_neoh.Piola2StiffnessGrowth2d = @piola2StiffnessGrowth2d;
    function [Svoigt,dSdCvoigt,dSdFg] = piola2StiffnessGrowth2d(Fe,Fg)
        mu = mat_neoh.mu;
        lambda = mat_neoh.lambda;
        
        J = det(Fe);
        Ce = Fe'*Fe;
        Cinv = inv(Ce);
        Fginv = inv(Fg);
        
        % 2nd piola kirchhoff stress
        Se = (lambda*log(J)-mu)*Cinv+mu*eye(2);
        S = Fginv*Se*Fginv';
        Svoigt = [S(1,1);S(2,2);S(1,2)];
        
        % 4th order elastic stiffness Le = dSe/dCe
        Le = zeros(2,2,2,2);
        for i=1:2; for j=1:2; for k=1:2; for l=1:2
            Le(i,j,k,l) = lambda               * Cinv(i,j)*Cinv(k,l) ...
                            +  (mu-lambda*log(J)) * (Cinv(i,k)*Cinv(j,l)+Cinv(i,l)*Cinv(k,j));     
        end; end; end; end
        %%{
        % !!! CHECK !!! dSdTH_C
        th = Fg(1,1);
        dSdTH_C = zeros(2);
        for r=1:2; for s=1:2
            dSdTH_C(i,j) = dSdTH_C(i,j)+Le(i,j,r,s)*Ce(r,s);
        end;end
        dSdTH_C = -dSdTH_C/th^3-2/th*S;
        %}
        % 4th order dS_ij/dFg_kl
        dSdFg = zeros(2,2,2,2);
        for i=1:2; for j=1:2; for k=1:2; for l=1:2
            % 1st contribution
            dSdFg(i,j,k,l) = dSdFg(i,j,k,l) -(Fginv(i,k)*S(l,j)  + Fginv(j,k)*S(i,l));
            
            % 2nd contribution, summation over p, q, r, s
            for p=1:2; for q=1:2; for r=1:2; for s=1:2
                 dSdFg(i,j,k,l) = dSdFg(i,j,k,l) - 0.5 * Fginv(i,p)*Fginv(j,q)*...
                                  Le(p,q,r,s)*(Fginv(l,r)*Ce(k,s)+Fginv(l,s)*Ce(r,k));
            end; end; end; end
        end; end; end; end

        %{
        % !!! CHECK !!! 4th order dS_ij/dFg_kl
        dSdFg_C = zeros(2,2,2,2);
        I = eye(2);
        for i=1:2; for j=1:2; for k=1:2; for l=1:2
            % 1st contribution
            dSdFg_C(i,j,k,l) = -( I(i,k)*S(l,j) + I(j,k)*S(i,l) )/th;    
            
            % 2nd contribution, summation over p, q, r, s
            for r=1:2; for s=1:2
                 dSdFg_C(i,j,k,l) = dSdFg_C(i,j,k,l) - 0.5 /th^3 *...
                                  Le(i,j,r,s)*(I(l,r)*Ce(k,s)+I(l,s)*Ce(r,k));
            end; end;
        end; end; end; end
        dSdFg_C
    
    
        % !!! CHECK !!! dSdTH_C2
        th = Fg(1,1);
        dSdTH_C2 = zeros(2);
        for i=1:2;for j=1:2;for k=1:2;for l=1:2;
            dSdTH_C2(i,j) = dSdTH_C2(i,j) + dSdFg(i,j,k,l)*I(k,l);
        end;end;end;end
        dSdTH_C2
        %}
        dSdFg = dSdTH_C;

        % Write dSdC in voigt notation
        vId1 = [1 2 1]; vId2 = [1 2 2];
        dSdCvoigt = zeros(3);
        for row = 1:3; for col = 1:3
                i = vId1(row); j=vId2(row); k=vId1(col); l=vId2(col);
                
                for p=1:2; for q=1:2; for r=1:2; for s=1:2
                dSdCvoigt(row,col) = dSdCvoigt(row,col) + 0.5*Fginv(i,p)*...
                                     Fginv(j,q)*Le(p,q,r,s)*Fginv(k,r)*Fginv(l,s);
                end; end; end; end
        end; end
    end


    mat_neoh.Piola2StiffnessGrowth = @piola2StiffnessGrowth;
    function [Svoigt,dSdCvoigt,dSdFg] = piola2StiffnessGrowth(Fe,Fg)
        mu = mat_neoh.mu;
        lambda = mat_neoh.lambda;
        
        J = det(Fe);
        Ce = Fe'*Fe;
        Cinv = inv(Ce);
        Fginv = inv(Fg);
        
        % 2nd piola kirchhoff stress
        Se = (lambda*log(J)-mu)*Cinv+mu*eye(3);
        S = Fginv*Se*Fginv';
        Svoigt = [S(1,1);S(2,2);S(3,3);S(2,3);S(1,3);S(1,2)];
        
        % 4th order elastic stiffness Le = dSe/dCe
        Le = zeros(3,3,3,3);
        for i=1:3; for j=1:3; for k=1:3; for l=1:3
            Le(i,j,k,l) = lambda               * Cinv(i,j)*Cinv(k,l) ...
                            +  (mu-lambda*log(J)) * (Cinv(i,k)*Cinv(j,l)+Cinv(i,l)*Cinv(k,j));     
        end; end; end; end

        % 4th order dS_ij/dFg_kl
        dSdFg = zeros(3,3,3,3);
        for i=1:3; for j=1:3; for k=1:3; for l=1:3
            % 1st contribution
            dSdFg(i,j,k,l) = -(Fginv(i,k)*S(l,j)  + Fginv(j,k)*S(i,l));    
            
            % 2nd contribution, summation over p, q, r, s
            for p=1:3; for q=1:3; for r=1:3; for s=1:3
                 dSdFg(i,j,k,l) = dSdFg(i,j,k,l) - 0.5 * Fginv(i,p)*Fginv(j,q)*...
                                  Le(p,q,r,s)*(Fginv(l,r)*Ce(k,s)+Fginv(l,s)*Ce(r,k));
            end; end; end; end
        end; end; end; end

        % Write dSdC in voigt notation
        vId1 = [1 2 3 2 1 1]; vId2 = [1 2 3 3 3 2];
        dSdCvoigt = zeros(6);
        for row = 1:6; for col = 1:6
                i = vId1(row); j=vId2(row); k=vId1(col); l=vId2(col);
                
                for p=1:3; for q=1:3; for r=1:3; for s=1:3
                dSdCvoigt(row,col) = dSdCvoigt(row,col) + 0.5*Fginv(i,p)*...
                                     Fginv(j,q)*Le(p,q,r,s)*Fginv(k,r)*Fginv(l,s);
                end; end; end; end
        end; end
    end


end
