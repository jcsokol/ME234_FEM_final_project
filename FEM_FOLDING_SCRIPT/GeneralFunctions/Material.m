function material = Material(materialType,matProp)
    

    %Initialize variables
    material.localID = -1;
    material.type = materialType;
    
    % Create the specific element type
    switch material.type
        case 'MAT_NEOH'
            material.material = MAT_NEOH(matProp);
        case 'MAT_LE'
            material.material = MAT_LE(matProp);
        otherwise
            material.material = [];
            error(['ERROR: material type ', material.type, ' is not supported']);
    end

end
