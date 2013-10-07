function animatesytem_TM(design, simoptions, T, Y, results, frames, fps, avifilename, mwscale)
% Animates the snapper device given some vectors of displacement
%
% Input:
%
%   Poles - The number of armature Poles to draw in the simulation
%
%   Wp - The pole pitch of the machine
%
%   xA - vector of armature positions at the corresponding times in 'T'
%
%   xT - vector of translator positions at the corresponding times in 'T'
%
%   frames - the total number of frames to create in the animation
%
%   fps - (optional) scalar number of frames per second in the animation if it
%         is to be an avi file
%
%   avifilename - (optional) string containing the name of a file in which
%                 to store a video of the animation. If provided, the
%                 animation will not be displayed, but frames taken from a
%                 hidden figure will instead be stored in an avi file.
%
% Output:
%
%   No output
%

    if nargin > 7
        %scrsz = get(0,'ScreenSize');
        h=figure('visible','off','units','normalized','outerposition',[0 0 1 1],'Color','w'); %'Position',[1, 1, 3*scrsz(3)/4, 3*scrsz(4)/4]); %turns visibility of figure off
    else
        h=figure('visible','on','units','normalized','outerposition',[0 0.03 1 0.97],'Color','w');
    end  
    
    if nargin < 9
        mwscale = 1;
    end
    
    axes;
    
    maxyT = max(abs(results.xT));
    
    maxxB = max(abs(Y(:,2)));
    
    maxyB = max(Y(:,1));
    
    xOffset = 0;
    yOffset = 0;
    
    linearT = min(T):(max(T)-min(T))/frames:max(T);
    
    xT = interp1(T, results.xT, linearT);
    
    xBh = interp1(T, Y(:,1), linearT);
    
    xBs = interp1(T, Y(:,3), linearT);
    
    % plotscene_snapper(h, design, simoptions, xB, yB, yT, yA, maxyA, maxyT, yOffset, xOffset)
    plotsystemscene_TM(design, simoptions, xBs(1), xBh(1), xT(1), maxyT, yOffset, xOffset, mwscale);
    
    % now determine suitible axis limits to keep everything to scale
    Axmax = max(simoptions.BuoyParameters.a + max(simoptions.BuoyParameters.a/2, 1.1*max(abs(Y(:,1)))), ...
        simoptions.BuoyParameters.draft*2 + simoptions.tether_length + abs(maxyT) + max(abs(Y(:,3))));
    
    Aymax = Axmax;
    
    Axmin = -Axmax;
    
    hTrans = design.Poles(2) * design.Wp;

    Aymin = -maxyT - hTrans/2;
    
    % set the axis limits
    axis([Axmin Axmax Aymin Aymax]+xOffset);
    
    if nargin > 7
        aviobj=avifile(avifilename, 'fps', fps, 'compression', 'Cinepak'); %creates AVI file
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'MP43'); %creates AVI file
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'MSVC'); %creates AVI file
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'DIVX'); %creates AVI file
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'RLE'); %creates AVI file
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'MRLE'); %creates AVI file
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'cvid');
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'MPG4');
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'UYVY');
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'XVID');
        %aviobj=avifile(avifilename, 'fps', fps, 'compression', 'none');
    end
    
    for i = 2:length(linearT)
        
        % Draw the snapper scene
        % plotscene_snapper(h, design, simoptions, xB, yB, yT, yA, maxyA,
        % maxyT, yOffset, xOffset)
        plotsystemscene_TM(design, simoptions, xBs(i), xBh(i), xT(i), maxyT, yOffset, xOffset, mwscale);
        %pause;
        
        if nargin > 7
            % Set the axis height and width to suitible levels
            axis([Axmin Axmax Aymin Aymax] + xOffset);
            
            try
                aviobj=addframe(aviobj,h); %adds frames to the AVI file
            catch ERR
                disp('Error encountered, closing avi file')
                aviobj=close(aviobj) %closes the AVI file
                rethrow(ERR);
            end
            
        else
%             % Pause for some time
%             pause(linearT(i)-linearT(i-1));
            pause(1/fps);
        end

        cla
    end
    
    if nargin > 7
        disp('Closing avi file')
        
        aviobj=close(aviobj) %closes the AVI file
        
        close(h); %closes the invisible figure
    end

end