% Test_rectregionsyperiodic

x = 1;
y1 = 0.8;
y2 = 0.2;
xoffset = 0;
tol = 1e-5;

%% cases

test = 1;

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

[nodes, nodeids, links, rectcentres, spacecentres] = rectregionsyperiodic(x, y1, y2, xoffset, ydisp, tol);

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



