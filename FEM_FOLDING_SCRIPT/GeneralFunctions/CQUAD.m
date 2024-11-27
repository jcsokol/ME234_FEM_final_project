function cquad = CQUAD(nodes,element)
    % Initialize variables
    cquad.nodes = nodes(1:4);
    cquad.element = element;
    
    % Rotation matrix T
    ex = element.direction;
    th = atan2(ex(2),ex(1));
    ey = [-sin(th) cos(th)]';
    
    TNode = [ex';ey'];
    cquad.T = blkdiag(TNode,TNode,TNode,TNode);
    
    % Functions
    cquad.buildElementMatrix = @buildElementMatrix;
    function [mat, var] = buildElementMatrix(u,dt,Fg0)
        ve = cquad.element.dof(u);
        veE = cquad.T*ve;                   % ve in element coordinates
        r0 = cquad.element.loc();
        r0E = cquad.T*r0(:);                % r0 in element coordinates
        
        gp = [-1 -1; 1 -1; 1 1; -1 1];
        wp = [1 1 1 1];
        ndim=2;
        mat = zeros(length(ve));
        var.Fg={};
        for ip = 1:4
            
            % Extract (r,s) of this integration point (coordinates in
            % element coordinates
            r = gp(ip,1)/sqrt(3);
            s = gp(ip,2)/sqrt(3);
            
            % Compute deformation gradient and strain derivatives
            [F, ~,~,J,~,dNx] = DeformationStrain_H2(r0E(:), r, s, veE(:));
            
            % Compute first p
            [~,A, var.Fg{ip}]  = element.property.property.Piola1Stiffness(F,Fg0{ip},dt,ndim);
            
            % Compute using A
            dNx = dNx';
            for i=1:4; for j=1:i
                eni=(i-1)*ndim;
                enj=(j-1)*ndim;
                
                for p=1:ndim; for q=1:ndim; for r=1:ndim; for s=1:ndim
                        mat(eni+p,enj+q) = mat(eni+p,enj+q) + dNx(r,i)*A(p,r,q,s)*dNx(s,j) *det(J)*wp(ip);
                end;end;end;end
            end;end
        end
        % Make symmetric
        mat = tril(mat)+tril(mat,-1)';
        % Transform back into global x,y,z coordinats
        mat = cquad.T'*mat*cquad.T;
        
    end

    cquad.buildInternalForceVector = @buildInternalForceVector;
    function vec = buildInternalForceVector(u,Fg)
        ve = cquad.element.dof(u);
        veE = cquad.T*ve;                   % ve in element coordinates
        r0 = cquad.element.loc();
        r0E = cquad.T*r0(:);                % r0 in element coordinates
        
        gp = [-1 -1; 1 -1; 1 1; -1 1];
        wp = [1 1 1 1];
        ndim = 2;
        vec = zeros(length(ve),1);
        for ip = 1:4
            
            % Extract (r,s,t) of this integration point (coordinates in
            % element coordinates
            r = gp(ip,1)/sqrt(3);
            s = gp(ip,2)/sqrt(3);
            
            % Compute deformation gradient and strain derivatives
            [F, ~, ~,J,~,dNx] = DeformationStrain_H2(r0E(:), r, s, veE(:));
            
            % Compute 1st piola kirchhof stress
            P  = element.property.property.Piola1Stiffness(F,Fg{ip},0,ndim);
            
            % Compute vec using P
            dNx = dNx';
            for i=1:4
                en=(i-1)*ndim;
                for p=1:ndim; for q=1:ndim;
                    vec(en+p) = vec(en+p) + P(p,q)*dNx(q,i)'*det(J)*wp(ip);
                end;end
            end
        end
        
        % Transform back into global x,y,z coordinats
        vec = cquad.T'*vec;
        
    end
end
