function IVars = beamvars(IMethod, n)

    switch IMethod
        
        case '1.3'
            
            % Load the database of available beams
            % beams columns: t, d, b
            beams = rectangularsections;
            
            % Change to metres from mm
            beams = beams ./ 1000;

            t = beams(:,1);
            d = beams(:,2);
            di = d - (2*t);
            b = beams(:,3);
            bi = b - (2*t);
            
            IVars = [d, b, di, bi];
            
        case '1.6'
            
            % Load the database of available I-Beams
            % beams columns: h, b, tw, t, r, d
            beams = Ibeamsections;

            % Change to metres from mm
            beams = beams ./ 1000;

            b = beams(:,2);
            t = beams(:,4);
            tw = beams(:,3);
            d = beams(:,1)-beams(:,4);
            
            IVars = [b, t, tw, d];
            
        otherwise
                
            error('Beam type not supported');
            
    end
    
    if nargin > 1
        IVars = IVars(n,:);
    end
    
end