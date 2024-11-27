function [nodaba eldaba cortex_elements subcortex_elements] = readAbaqusFile(nelc)

    % Extract nodal locations from input file
    filename =(['InputFiles/nc_',num2str(nelc),'.inp']);

    % Open file
    fid = fopen(filename);
    INP = fread(fid,'*char');
    INP = INP';

    % in the .inp file a node definition line always starts with '*Node'. So
    % find this string in the file
        iNODE       = strfind(INP,'*Node');
        iELEM       = strfind(INP,'*Element'); iELEM1 = iELEM+strfind(INP(iELEM+1:end),'*');
        iCORTEX     = strfind(INP,'*Elset, elset=RECTANGLE_FACE-CORTEX, generate');
        iCORTEX1 = iCORTEX+strfind(INP(iCORTEX+1:end),'*');
        iSUBCORTEX  = strfind(INP,'*Elset, elset=RECTANGLE_FACE-SUBCORTEX, generate');
        iSUBCORTEX1 = iSUBCORTEX+strfind(INP(iSUBCORTEX+1:end),'*');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% NODES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Take the part of the file which contains nodal information
        B = INP(iNODE+6:iELEM-1);
        B_data = textscan(B, '%d %f %f ','delimiter',',');
        nnode = size(B_data{1,1},1);
    % Create nodal database: [nID,X,Y]
        nodaba = zeros(nnode,3);
        nodaba(:,1) = B_data{1,1};
        nodaba(:,2) = B_data{1,2};
        nodaba(:,3) = B_data{1,3};

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ELEMENTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Take the part of the file which contains element information
    B = INP(iELEM:iELEM1-1);
    % Remove the first line from B, which contains *Element, type=CPS4R
        B_top = textscan(B, '%s', 2, 'delimiter', ',');
        B_top = B_top{1,1}(end);
        iELEM_1 = strfind(B,char(B_top))+length(char(B_top));
        B = B(iELEM_1:end);

    % Get the element id's, with the element nodes
        B_data = textscan(B, '%d %d %d %d %d','delimiter',',');
        nelem = size(B_data{1,1},1);

    % Create element database: [eID n1 n2 n3 n4]
        eldaba = zeros(nelem,5);
        eldaba(:,1) = B_data{1,1};
        eldaba(:,2) = B_data{1,2};
        eldaba(:,3) = B_data{1,3};
        eldaba(:,4) = B_data{1,4};
        eldaba(:,5) = B_data{1,5};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CORTEX ELEMENTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Take the part of the file which contains element information
    B = INP(iCORTEX:iCORTEX1-1);
    % Remove the first line from B, which contains *Elset, elset=ELLIPSE_FACE-CORTEX, generate
        B_top = textscan(B, '%s', 3, 'delimiter', ',');
        B_top = B_top{1,1}(end);
        iELEM_1 = strfind(B,char(B_top))+length(char(B_top));
        B = B(iELEM_1:end);
        
    % Get the cortex element id's, with the element nodes
        B_data = textscan(B, '%d %d %d','delimiter',',');
        cortex_elements = B_data{1}:B_data{3}:B_data{2};
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SUBCORTEX ELEMENTS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Take the part of the file which contains element information
    B = INP(iSUBCORTEX:iSUBCORTEX1-1);
    % Remove the first line from B, which contains *Elset, elset=ELLIPSE_FACE-SUBCORTEX, generate
        B_top = textscan(B, '%s', 3, 'delimiter', ',');
        B_top = B_top{1,1}(end);
        iELEM_1 = strfind(B,char(B_top))+length(char(B_top));
        B = B(iELEM_1:end);
        
    % Get the cortex element id's, with the element nodes
        B_data = textscan(B, '%d %d %d','delimiter',',');
        subcortex_elements = B_data{1}:B_data{3}:B_data{2};
end