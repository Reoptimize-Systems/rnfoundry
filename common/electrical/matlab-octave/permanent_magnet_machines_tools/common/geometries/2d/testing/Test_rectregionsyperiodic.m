% Test_rectregionsyperiodic

x = 1;
y1 = 0.8;
y2 = 0.2;
xoffset = 0;
tol = 1e-5;

%% cases

test = 3;

switch test

    case 1
        % no displacement
        ydisp = 0;

    case 2
        % B top on top
        ydisp = y2/2;

    case 3
        % A bot on bottom
        ydisp = -y2/2;

    case 4
        % halway B
        ydisp = (y1/2 + y2/2);

    case 5
        % halway A
        ydisp = -(y1/2 + y2/2);

    case 6
        % B bot on bot 
        ydisp = (y1 + y2/2);

    case 7
        % A top on top
        ydisp = -(y1 + y2/2);
        
    case 8
        % one period translation
        % A bot B top
        ydisp = 2*(y1 + y2);       
        
    case 9
        % -ve one period translation
        % A bot B top
        ydisp = -2*(y1 + y2);   
        
    case 10
        % half period translation
        % A top B bot
        ydisp = (y1 + y2);     
        
    case 11
        % -ve half period translation
        % A top B bot
        ydisp = -(y1 + y2);           
        
        
end

%%

[nodes, nodeids, links, rectcentres, spacecentres] = rectregionsyperiodic(x, y1, y2, xoffset, ydisp, 'Tol', tol);

%%

close all 

linkstart = 1;
linkend = size(links,1);
plotnodelinks(nodes, links(linkstart:linkend, :));
hold on
for i = 1:size(rectcentres, 1)
    if rectcentres(i,3)
        text(rectcentres(i,1), rectcentres(i,2), 'B')
    else
        text(rectcentres(i,1), rectcentres(i,2), 'A')
    end
end
plot(spacecentres(:,1), spacecentres(:,2), '+m');
hold off



%% multi-pole

x = 1;
y1 = 0.8;
y2 = 0.2;
xoffset = 0;
tol = 1e-5;

%% cases

test = 13;

switch test

    case 1
        % no displacement
        n = 3
        ydisp = n*2*(y1 + y2);

    case 2
        % B top on top
        n = 2
        ydisp = n*2*(y1 + y2) + y2/2;

    case 3
        % A bot on bottom
        n = 1
        ydisp = n*2*(y1 + y2) -y2/2;

    case 4
        % halway B
        n = 1
        ydisp = n*2*(y1 + y2) + (y1/2 + y2/2);

    case 5
        % halway A
        n = 2
        ydisp = n*2*(y1 + y2) -(y1/2 + y2/2);

    case 6
        % B bot on bot
        n = 1;
        ydisp = n*2*(y1 + y2) + (y1 + y2/2);

    case 7
        % A top on top
        n = 2
        ydisp = n*2*(y1 + y2) -(y1 + y2/2);
        
    case 8
        % one period translation
        % A bot B top
        ydisp = 2*(y1 + y2);       
        
    case 9
        % -ve one period translation
        % A bot B top
        ydisp = -2*(y1 + y2);   
        
    case 10
        % half period translation
        % A top B bot
        n = 0
        ydisp = n*2*(y1 + y2) + (y1 + y2);     
        
    case 11
        % -ve half period translation
        % A top B bot
        ydisp = -(y1 + y2);
        
        
    case 12
        % error
        ydisp = 1.28;
        
    case 13
        ydisp = 960.0000e-003;
        
        
end

%%

[nodes, nodeids, links, rectcentres, spacecentres] = rectregionsyperiodic(x, y1, y2, xoffset, ydisp, 'Tol', tol, 'NodeCount', 0, 'NY1Pairs', 4);

%%

close all 

linkstart = 1;
linkend = size(links,1);
plotnodelinks(nodes, links(linkstart:linkend, :));
hold on
for i = 1:size(rectcentres, 1)
    if rectcentres(i,3)
        text(rectcentres(i,1), rectcentres(i,2), 'B')
    else
        text(rectcentres(i,1), rectcentres(i,2), 'A')
    end
end
plot(spacecentres(:,1), spacecentres(:,2), '+m');
hold off



%% video


close all 

Ny1p = 3;
totaldisp =  Ny1 * (y1 + y2);
ydisp = -totaldisp:totaldisp/100:totaldisp;
% figure
for n = 1:numel(ydisp)
    
    [nodes, nodeids, links, rectcentres, spacecentres] = rectregionsyperiodic(x, y1, y2, xoffset, ydisp(n), 'Tol', tol, 'NodeCount', 0, 'NY1Pairs', Ny1p);
    
    linkstart = 1;
    linkend = size(links,1);
%     figure
    plotnodelinks(nodes, links(linkstart:linkend, :));
    hold on
    for i = 1:size(rectcentres, 1)
        if rectcentres(i,3) == 1
            text(rectcentres(i,1), rectcentres(i,2), 'B')
        else
            text(rectcentres(i,1), rectcentres(i,2), 'A')
        end
    end
    plot(spacecentres(:,1), spacecentres(:,2), '+m');
    hold off
    
    xlim ([-1, 1]);
    ylim ([-1, Ny1+1]);
    
    pause (0.2)
    
    cla
    
    

end