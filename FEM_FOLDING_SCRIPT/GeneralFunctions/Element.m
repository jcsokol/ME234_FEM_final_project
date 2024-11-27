function element = Element(elementType,nodes,elemDir,property)
    

    %Initialize variables
    element.matrix     = [];
    element.geomMatrix = [];
    element.type = elementType;
    
    
    element.property = property;
    element.localID = -1;
    
    % Function declarations
    element.dofID = @dofID;
    element.dofPhiID = @dofPhiID;
    element.dof = @dof;
    element.loc = @loc;
    element.direction = elemDir;
    
    % Variables that change during the simulation
    element.Fg0 = {};                       %Fg at beginning of step (element coordinates)
    element.Fg  = {};                       %Fg during iteration (element coordinates)
    for j=1:8
        element.Fg0{j} = eye(3);
        element.Fg{j} = eye(3);
    end
    
    element.X  = {};                       % heat flux during iteration (element coordinates)
    for j=1:8
        element.X{j} = zeros(3,1);
    end
    
    element.nc0  = {};                      % number of connectivities (for multifield models)
    element.nc  = {};                       % number of connectivities (for multifield models)
    for j=1:8
        element.nc0{j} = 0;
        element.nc{j} = 0;
    end

    
    % Create the specific element type
    switch element.type
        case 'CHEXA'
            element.element = CHEXA(nodes,element);
        case 'CQUAD'
            element.element = CQUAD(nodes,element);
        case 'CQUAD_MF'
            element.element = CQUAD_MF(nodes,element);
        case 'CQUAD_DF'
            element.element = CQUAD_DF(nodes,element);
        case 'CQUAD7'
            element.element = CQUAD7(nodes,element);
        case 'CQUAD8'
            element.element = CQUAD8(nodes,element);
        case 'CTRIA_DF'
            element.element = CTRIA_DF(nodes,element);
        otherwise
            element.element = [];
            fprintf(2,['ERROR: element type ', element.type, ' is not supported.\n']);
    end
    
    
    
    
    
    % Functions implementations
    function dofID = dofID()
        dofID = [];
        for i=1:length(element.element.nodes)
            dofID = [dofID element.element.nodes(i).dofID];
        end  
    end
    function dofID = dofPhiID()
        dofID = [];
        for i=1:length(element.element.nodes)
            dofID = [dofID element.element.nodes(i).dofPhiID];
        end  
    end


    function dof = dof(u)
        dof = u(dofID());
    end

    function loc = loc()
        loc = [];
        for i=1:length(element.element.nodes)
            loc = [loc element.element.nodes(i).loc];
        end  
    end

end
