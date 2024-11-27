function node = Node(loc)

    node.loc = loc;
    node.x = loc(1);
    node.y = loc(2);
    if(length(loc)==3) node.z=loc(3); else node.z=0; end
    
    node.dofID = [];        % ID of nodal dof in global matrices/vectors
    node.localID = -1;
    
    
    % Function declarations
end