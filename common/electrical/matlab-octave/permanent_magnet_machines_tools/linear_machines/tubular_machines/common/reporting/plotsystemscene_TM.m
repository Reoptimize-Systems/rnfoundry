function plotsystemscene_TM(design, simoptions, xB, yB, yT, maxyT, yOffset, xOffset, mwscale)

    if nargin < 9
        mwscale = 1;
    end
    
    hTrans = design.Wp * design.Poles(1);
    
    if isfield(design, 'minArmLength')
        hArm = design.minArmLength;
    else
        hArm = design.Poles(2) * design.Wp;
    end
    
    wArm = mwscale * design.Ra - design.Ri;
    
    wTrans = mwscale * 2*design.Rm;
    
    % draw left armature part
    x = -wArm - design.g - wTrans/2 + xOffset;
    
    y = 0 - (hArm/2) + yOffset;
    
    w = wArm;
    
    h = hArm;
    
    position = [x,y,w,h];
    
    rectangle('Position',position,...
          'Curvature', [0,0],...
         'LineWidth', 1)
     
    hold on;
    
    % draw right armature part
    x = wTrans/2 + design.g + xOffset;
    
    y = 0 - (hArm/2) + yOffset;
    
    w = wArm;
    
    h = hArm;
    
    position = [x,y,w,h];
    
    rectangle('Position',position,...
          'Curvature', [0,0],...
         'LineWidth', 1)
     
     
    % draw the translator 
    x = -wTrans/2 + xOffset;
    
    y = yT - (hTrans/2) + yOffset;
    
    w = wTrans;
    
    h = hTrans;
    
    position = [x,y,w,h];
    
    rectangle('Position',position,...
          'Curvature', [0,0],...
         'LineWidth', 1)
     
    % draw the buoy 
    
    hBuoy = simoptions.BuoySim.BuoyParameters.draft * 2;
    wBuoy = simoptions.BuoySim.BuoyParameters.a * 2;
    
    hBStart = simoptions.BuoySim.tether_length + maxyT + yOffset;
        
    x = -wBuoy/2 + xB + xOffset;
    
    y = yB + hBStart;
    
    w = wBuoy;
    
    h = hBuoy;
    
    position = [x,y,w,h];
    
    rectangle('Position',position,...
          'Curvature', [0.1,0.4],...
          'LineWidth', 1)
     

    % Finally draw the tether 
    tetherPoints = [xOffset, (hTrans/2) + yT + yOffset;
                    xOffset, (hTrans/2) + maxyT + yOffset
                    xOffset + xB, hBStart + yB];
                
    plot(tetherPoints(:,1), tetherPoints(:,2));
    
    hold off;
    
end
