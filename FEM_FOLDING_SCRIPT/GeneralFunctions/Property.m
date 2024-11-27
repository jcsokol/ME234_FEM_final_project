function property = Property(propertyType,material,var)
    

    %Initialize variables
    property.localID = -1;
    property.type = propertyType;
    
    % Create the specific element type
    switch property.type
        case 'PSOLID10'
            property.property = PSOLID10(material,var);
        case 'PSOLID11'
            property.property = PSOLID11(material,var);
        case 'PSOLID12'
            property.property = PSOLID12(material,var);
        case 'PSOLID20'
            property.property = PSOLID20(material,var);
        case 'PSOLID21'
            property.property = PSOLID21(material,var);
        case 'PSOLID30'
            property.property = PSOLID30(material,var);
        case 'PSOLID71'
            property.property = PSOLID71(material,var);
        case 'PSOLID81'
            property.property = PSOLID81(material,var);
        otherwise
            property.property = [];
            error(['ERROR: property type ', property.type, ' is not supported ']);
    end

end