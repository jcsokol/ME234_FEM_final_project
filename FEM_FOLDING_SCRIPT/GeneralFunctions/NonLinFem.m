%% Computations
%matlabpool ('open',4);      % Call to open the distributed processing
ndof = nnode*param.numDofPerNode;
u = zeros(ndof,1);
u0 = u;
outputVar.time = 0;
outputVar.frameId = 0;
outputMov = plotOutput(elements,u,folders,outputVar);
tic

step = 1;
time = dt0;
dt   = dt0;
timeFrac = time/tEnd;
dtFrac = dt/tEnd;
while(dt>0 && step < max_step)

    iter = 0;
    residual = 1;
    reset = false;
    while (iter<1 || (residual>param.tolerance && iter < max_iter))
        iter = iter+1;

        %%% Assemble Global Matrix
        K = sparse(ndof,ndof);
        du = sparse(ndof,1);
        
        %%{
        for i=1:nel
            dofID = elements(i).dofID();
            
            [kel, var_el] = elements(i).element.buildElementMatrix(u,dt,elements(i).Fg0);
            elements(i).Fg = var_el.Fg;
            
            K(dofID,dofID) = K(dofID,dofID) + kel;
        end
        %}
        %{
        KEL= cell(nel,1);
        parfor i=1:nel
            dofID = elements(i).dofID();
            
            [kel, var_el] = elements(i).element.buildElementMatrix(u,dt,elements(i).Fg0);
            elements(i).Fg = var_el.Fg;
            
            KEL{i} = kel;
        end
        for i=1:nel
            dofID = elements(i).dofID();
            K(dofID,dofID) = K(dofID,dofID) + KEL{i};
        end
        %}
        %%% Assemble global residual in first iteration
        if (iter==1)
            R = sparse(ndof,1);
            % External loads
            for i=1:length(loads)
                dofID = loads(i).node.dofID;
                R(dofID) = R(dofID)+timeFrac*loads(i).loadVec(:);
            end

            % Minus internal loads
            for i=1:nel
                dofID = elements(i).dofID();
                R(dofID) = R(dofID) - elements(i).element.buildInternalForceVector(u,elements(i).Fg);
            end
        end
        %%% Apply SPC
        dofp = [];
        dup   = [];
        for i =1:length(spc)
            dof = spc(i).node.dofID(spc(i).dofIDs);
            numdof = length(dof);
            value = dtFrac*spc(i).value;

            dofp = [dofp; dof(:)];
            if(iter==1)
                dup   = [dup; value*ones(numdof,1)]; 
            else
                dup   = [dup; zeros(numdof,1)];
            end
        end

        dofr = [1:ndof]'; dofr(dofp)=[];

        Krr = K(dofr,dofr);
        Krp = K(dofr,dofp);

        Rr = R(dofr);


        %%% Solve system
        dur = Krr\(Rr-Krp*dup);
        du(dofp) = dup;
        du(dofr) = dur;
        u = u+du;


        %%% Assemble global residual
        R = sparse(ndof,1);
        % External loads
        for i=1:length(loads)
            dofID = loads(i).node.dofID;
            R(dofID) = R(dofID)+timeFrac*loads(i).loadVec(:);
        end

        % Minus internal loads
        %%{
        for i=1:nel
            dofID = elements(i).dofID();
            R(dofID) = R(dofID) - elements(i).element.buildInternalForceVector(u,elements(i).Fg);
        end
        %}
        %{
        REL= cell(nel,1);
        parfor i=1:nel
            REL{i} = - elements(i).element.buildInternalForceVector(u,elements(i).Fg);
        end
        for i=1:nel
            dofID = elements(i).dofID();
            R(dofID) = R(dofID) + REL{i};
        end
        %}
        
        
        
        
        
        residual = norm(R(dofr))/ndof;
        sim_deflection = max(u(2:2:end));
        fprintf(['Step: ', num2str(step),'  \t Iter: ', num2str(iter),'  \t dt: ', num2str(dt),'  \t t: ', num2str(time),'  \t Residual: ',num2str(residual),'  \t SimDef: ',num2str(sim_deflection),'\n']);
    end
    
    
    % Compute new time step
    if(iter==max_iter && iter~=1 && residual > param.tolerance)
        fprintf(2,'Too many iterations, time step decreased.\n');
        reset = true;
    end
    
    
    if(reset)
        
        time = time -dt;
        dt = dt /2 ;
        time = time +dt;
        u = u0;
        if(dt<dtMin) 
            fprintf(2,'ERROR: time increment too small');
            break;
        end
        
    else
        outputVar.time = time;
        outputVar.frameId = length(outputMov);
        outputMov(end+1) = plotOutput(elements,u,folders,outputVar);
        
        if(iter<=max_iter_inc)
            dt = dt*1.25;
        end
        if(dt>dtMax) dt = dtMax; end
        if(time+dt>tEnd); dt = tEnd-time;end
        
        step = step+1;
        time = min(time + dt,tEnd);
        
        u0 = u;

        % Update Fg0 in all elements
        for i=1:nel; elements(i).Fg0 = elements(i).Fg;end
    end
    timeFrac = time/tEnd;
    dtFrac = dt/tEnd;
    
    toc
    
end

%matlabpool ('close');      % Call to close the distributed processing

% Create AVI file for stress plot
writeObj = VideoWriter([folders.Output,'/folding.avi']);
set(writeObj,'FrameRate',step/outputVar.movDuration);
open(writeObj);
writeVideo(writeObj,outputMov);
close(writeObj);
