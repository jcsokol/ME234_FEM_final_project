classdef printFunctions
    methods (Static = true)
        % Nodes
        function Nodes(nodes)
            for i=1:length(nodes)
                fprintf(['Node ',num2str(nodes(i).localID),':  \t x = ', num2str(nodes(i).x,'%1.3f'),'\t y = ', num2str(nodes(i).y,'%1.3f'),'\t z = ', num2str(nodes(i).z,'%1.3f'),'\n']); 
            end
        end
        
        % Elements
        function Elements(elements)
            for i=1:length(elements)
                fprintf(['Element ',num2str(elements(i).localID),':  \t Type: ',elements(i).type,'\t Nodes: ']); 
                
                for n = 1:length(elements(i).element.nodes)
                    fprintf([num2str(elements(i).element.nodes(n).localID), ' ']);
                end

                fprintf(['\t Property: ',num2str(elements(i).property.localID),'\n']);


            end
        end
        
        % SPC
        function SPC(spc)
            for i=1:length(spc)
                fprintf(['SPC ',num2str(i),':  \t Node: ', num2str(spc(i).node.localID),'\t presValue: ', num2str(spc(i).value,'%1.3f'), '\t presDOF: ']); 
                
                for n = 1:length(spc(i).dofIDs)
                    fprintf([num2str(spc(i).dofIDs(n)), ' ']);
                end

                fprintf('\n');


            end
        end
        
        % LOAD
        function LOAD(load)
            for i=1:length(load)
                fprintf(['LOAD ',num2str(i),':  \t Node: ', num2str(load(i).node.localID),'\t loadVec: ']); 
                
                for n = 1:length(load(i).loadVec)
                    fprintf([num2str(load(i).loadVec(n),'%1.3f'), ' ']);
                end

                fprintf('\n');


            end
        end
        
        
        
        
    end
end