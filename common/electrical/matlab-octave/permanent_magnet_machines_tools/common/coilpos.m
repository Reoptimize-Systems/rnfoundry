function pos = coilpos(Phases)
% chooses the spacing between Phases for a group of n Phases in an
% electrical machine
%
% Syntax
%
% pos = coilpos(Phases)
%
% Description
%
% pos = coilpos(Phases) produces a vector of normalised coil spacing for a
% set of 'Phases' Phases in an electrical machine. These positions are
% normalised to a single magnetic pole width.
%
% Example
%
% pos = coilpos(3)
% 
% pos =
% 
%      0.0000e+000
%    666.6667e-003
%      1.3333e+000
%
% 

    switch Phases
        
        case 1
            
            pos = 0;
            
        case 2
            
            pos = [0, 0.5]';
            
        case 3
            
            pos = [0, 2/3, 4/3]';
            
        otherwise
            
            error('This number of Phases has not yet been implemented');
            
    end

end