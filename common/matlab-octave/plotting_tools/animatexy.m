function animatexy(T, X, Y, speedfac, fps, tailtime)

    duration = max(T) - min(T);
    
    if nargin < 4
        speedfac = 1;
    end
    
    if nargin < 5
        fps = 25;
    end
    
    if nargin < 6
        tailtime = 0.05 * duration;
    end
    
    xmin = 1.1*min(X(:));
    
    xmax = 1.1*max(max(X(:)), eps);
    
    ymin = 1.1*min(Y(:));
    
    ymax = 1.1*max(max(Y(:), eps));

    figure('visible', 'on', 'units', 'normalized', 'outerposition', [0 0.03 1 0.97], 'Color', 'w');
    
    axis([xmin, xmax, ymin, ymax]);
    
    tstep = 1 / fps;
    
    linearT = min(T):tstep:max(T);

    X = interp1(T, X, linearT);
    
    Y = interp1(T, Y, linearT);
    
    if size(X,1) == 1
        X = X';
        Y = Y';
    end
    
    for i = 1:numel(linearT)
        
        tailcount = max(1, i - ceil(tailtime / tstep));
        
        plot(X(tailcount:i, :), Y(tailcount:i, :));
        
        axis([xmin, xmax, ymin, ymax]);
        
        pause(1 / speedfac / fps);
        
        cla;
    
    end
    

end