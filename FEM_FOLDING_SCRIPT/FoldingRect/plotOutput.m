function frame = plotOutput(elements,u,folders,var)

    scrsz = get(0,'ScreenSize');
    Fig1 = figure(1);
    clf
    set(Fig1,'color','white');
    wFig= scrsz(3)/1.5;
    hFig= 0.5*wFig;
    set(Fig1,'Position',[scrsz(3)/2-wFig/2 scrsz(4)/2-hFig/2 wFig hFig])
    hold on
    axis off
    axis equal
    ylim([0 1.5*var.sizeY])
    
    % All quad elements
    for i=1:length(elements)
        el = elements(i);
        n1 = el.element.nodes(1);   n2 = el.element.nodes(2);
        n3 = el.element.nodes(3);   n4 = el.element.nodes(4);

        plotFace(el,n1,n2,n3,n4);
    end
    text(-var.sizeX,-var.sizeY/10,0,['Time = ',num2str(var.time),' [days]'],'fontsize',14)
    
    
    % Store frame and save as eps
    print('-dpng', '-r300', ...
      [folders.Output,'/frame_', num2str(var.frameId,'%5.5d')])
  
    frame = getframe(gcf);
    
        
     function plotFace(el,n1,n2,n3,n4)
        
        % Compute color based on property
        if(~isempty(findstr(el.property.type,'PSOLID1')))
            faceColor = [255 102 102]/255;
        else
            faceColor = [255 204 153]/255;
        end
        
        % Deformed
        fill( [n1.x+u(n1.dofID(1)),n2.x+u(n2.dofID(1)),n3.x+u(n3.dofID(1)),n4.x+u(n4.dofID(1))],...
              [n1.y+u(n1.dofID(2)),n2.y+u(n2.dofID(2)),n3.y+u(n3.dofID(2)),n4.y+u(n4.dofID(2))],...
              faceColor,'FaceAlpha',1.0,'EdgeAlpha',0.1);
          
        % Plot element direction for subcortical elements (driven growth)
        if(~isempty(findstr(el.property.type,'PSOLID2')))
            % Element direction in current deformation
            elemDir = currentElemDir(el);
            
            % Center of element en size of vector
            locN = [n1.loc+u(n1.dofID)';
                    n2.loc+u(n2.dofID)';
                    n3.loc+u(n3.dofID)';
                    n4.loc+u(n4.dofID)'];
            locAvg = mean(locN);
            dloc = [max(abs(locAvg(1)-locN(:,1))); max(abs(locAvg(2)-locN(:,2)))];
            lVec = min(abs(dloc./elemDir))/3;
            
            % Plot
            plot( [locAvg(1)-lVec*elemDir(1) locAvg(1)+lVec*elemDir(1)], ...
              [locAvg(2)-lVec*elemDir(2) locAvg(2)+lVec*elemDir(2)],'k');
        end
        
        
    end

    function elemDir = currentElemDir(el)
       % Compute deformation gradient in center of element
       ve = el.dof(u);
       veE = el.element.T*ve;                   % ve in element coordinates
       r0 = el.loc();
       r0E = el.element.T*r0(:);                % r0 in element coordinates
       
       F = DeformationStrain_H2(r0E(:), 0,0, veE(:));
       elemDirE = F(:,1);
       elemDir = el.element.T(1:2,1:2)'*elemDirE;
        
    end


end