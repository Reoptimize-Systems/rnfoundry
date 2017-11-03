function position = planarrotorpos(ypole, position, fracpolepos, rotoranglepos)
% calculates the relative rotor position in units of fractions of a pole
% width
%
% Syntax
%
% position = planarrotorpos(ypole, position, fracpolepos, rotoranglepos)
%
% Input
%
%  ypole - pole width in the same units as 'position'
%
%  position - scalar relative position of field component to armature. Can
%    be empty, in which case fracpolepos is checked next. If not empty,
%    position is simply passed through to the output.
%
%  fracpolepos - scalar relative position of field component to armature as
%    a ratio of pole widths. Can be empty, in which case rotoranglepos is
%    checked next. Position is calculated as fracpolepos * ypole.
%
%  rotoranglepos - two element vector. The first number should be the
%    number of poles in the machine, the second the absolute position of a
%    rotor in radians.
%
% Output
%
%  position - the relative rotor position in units of fractions of a pole
%    width
%

% Copyright 2012 - 2014 Richard Crozier 


    if ~isempty(fracpolepos) && ~isempty(rotoranglepos) 
        
        error('ROTARY:axfluxrotor:badpos', 'Both fracpolepos and rotoranglepos specified.')
        
    elseif ~isempty(fracpolepos)
       
        position = fracpolepos * ypole;
        
    elseif ~isempty(rotoranglepos)
        
        if isvector(rotoranglepos) && numel(rotoranglepos) == 2
            
            if ~check.isint2eps(rotoranglepos(1))
                
                error('ROTARY:axfluxrotor:integerpoles', ...
                      ['The first value in rotoranglepos should be', ...
                       ' an integer denoting the total number of Poles in the machine.']);
                  
            else
                position = ypole * round(rotoranglepos(1)) * (rotoranglepos(2) / tau);
            end
        
        else
            error('ROTARY:axfluxrotor:badpos', ...
                 ['If specified, rotoranglepos must be a vector of two \n', ...
                  'values containing the number of Poles in the machine \n', ...
                  'and the rotor angle in radians respectively.']);
        end
        
    elseif isempty(position)
        
        error('Could not calculate rotor position');
        
    end
    
    
end