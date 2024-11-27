function elemDir = getElemDir(x,y,elemDir_option)

    % Create the specific element type
    switch elemDir_option
        case 'random'
            th = 2*pi*rand;
            elemDir = [cos(th); sin(th)];
        case 'horizontal'
            elemDir = [1;0];
        case 'vertical'
            elemDir = [0;1];
        case 'diagonal'
            elemDir = [1;1]/sqrt(2);
        case 'radial'
            elemDir = [x;y];
            elemDir = elemDir/norm(elemDir);
        otherwise
            error(['ERROR: elemDir option ', elemDir_option, ' is not supported ']);
    end


end