function [Qx, Qy] = bezierpoints(Px,Py,n)
% Px contains x-coordinates of control points [Px0,Px1,Px2,Px3]
% Py contains y-coordinates of control points [Py0,Py1,Py2,Py3]
% n is number of intervals
%


% % % --------------------------------
% % % Author: Dr. Murtaza Khan
% % % Email : drkhanmurtaza@gmail.com
% % % Modified By Richard Crozier
% % % --------------------------------

    dt  =  1/n;
    t = (2:n-1) * dt;
    
    order = numel(Px)-1;

    switch order
        
        case 1
            
            Px0 = Px(1);
            Py0 = Py(1);
            Px1 = Px(2);
            Py1 = Py(2);
            
            Qx = [Px(1), Px0 + t.*(Px1 - Px0), Px(end)];
            Qy = [Py(1), Py0 + t.*(Py1 - Py0), Py(end)];

        case 2
            
            Px0 = Px(1);
            Py0 = Py(1);
            Px1 = Px(2);
            Py1 = Py(2);
            Px2 = Px(3);
            Py2 = Py(3);

            Qx = [Px(1), ((1 - t).^2) .*Px0 + 2.*(1 - t).*t.*Px1 + (t.^2).*Px2, Px(end)];
            Qy = [Py(1), ((1 - t).^2) .*Py0 + 2.*(1 - t).*t.*Py1 + (t.^2).*Py2, Py(end)];
            
        case 3
            % Equation of cubic Bezier Curve, utilizes Horner's rule for efficient computation.
            % Q(t) = (-P0 + 3*(P1-P2) + P3)*t^3 + 3*(P0-2*P1+P2)*t^2 + 3*(P1-P0)*t + Px0
            Px0 = Px(1);
            Py0 = Py(1);
            Px1 = Px(2);
            Py1 = Py(2);
            Px2 = Px(3);
            Py2 = Py(3);
            Px3 = Px(4);
            Py3 = Py(4);

            cx3 = -Px0 + 3*(Px1-Px2) + Px3;
            cy3 = -Py0 + 3*(Py1-Py2) + Py3;
            cx2 = 3*(Px0-2*Px1+Px2); 
            cy2 = 3*(Py0-2*Py1+Py2);
            cx1 = 3*(Px1-Px0);
            cy1 = 3*(Py1-Py0);
            cx0 = Px0;
            cy0 = Py0;

            Qx = [Px(1), ((cx3 .* t + cx2) .* t + cx1) .* t + cx0, Px(end)];
            Qy = [Py(1), ((cy3 .* t + cy2) .* t + cy1) .* t + cy0, Py(end)];

        otherwise
            
            error('Bezier curves higher than order three not supported');

    end

end