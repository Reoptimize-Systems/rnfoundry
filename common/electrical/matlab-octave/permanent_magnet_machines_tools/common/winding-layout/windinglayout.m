function [coillayout, phaselayout] = windinglayout (NPhases, Qs, Poles, SL)
% calculates the layout of a winding using the star-of-slots method
%
% Syntax
%
% [coillayout, phaselayout] = windinglayout (NPhases, Qs, Poles, SL)
%
% Input
%
%  NPhases - scalar number of phases in the machine
%
%  Qs - scalar number of slots in the winding
%
%  Poles - scalar number of poles in the machine (note NOT the number of
%    pole-pairs, but 2x the number of pole-pairs)
%
%  SL - true/false flag indicating if this is single layered or
%    double-layered winding
%
% Output
%
%  coillayout - matrix containing the coil layout for each phase in terms
%    of the locations of each coil side. Each phase is represented as two
%    columns where each row is a coil, and the first column is the slot
%    number the first coil side lies in and the second the slot the other
%    coil side lies in.
%
%  phaselayout - matrix (or column vector) describing how the phases are
%    wound in the machine. Each column represents a the layers of all
%    slots, and each row the phase number and winding direction in that
%    layer of the slot.
%
% Example
%
% find the winding layout for a 4 pole 3 phase integral slot machine with
% overlapping windings. In this case there will be three slots per pole and
% the same number of slots as coils. coillayout will contain six columns,
% two for each phase, and phaselayout will contain two columns, one for
% each coil layer in the slots.
%
% >> [coillayout, phaselayout] = windinglayout (3, 4*3, 4, 0)
% 
% coillayout =
% 
%      1.0000e+000     4.0000e+000     3.0000e+000     6.0000e+000     5.0000e+000     8.0000e+000
%      7.0000e+000    10.0000e+000     9.0000e+000    12.0000e+000    11.0000e+000     2.0000e+000
%      4.0000e+000     7.0000e+000     6.0000e+000     9.0000e+000     2.0000e+000     5.0000e+000
%     10.0000e+000     1.0000e+000    12.0000e+000     3.0000e+000     8.0000e+000    11.0000e+000
% 
% 
% phaselayout =
% 
%      1.0000e+000    -1.0000e+000
%      3.0000e+000    -3.0000e+000
%      2.0000e+000    -2.0000e+000
%      1.0000e+000    -1.0000e+000
%      3.0000e+000    -3.0000e+000
%      2.0000e+000    -2.0000e+000
%      1.0000e+000    -1.0000e+000
%      3.0000e+000    -3.0000e+000
%      2.0000e+000    -2.0000e+000
%      1.0000e+000    -1.0000e+000
%      3.0000e+000    -3.0000e+000
%      2.0000e+000    -2.0000e+000
%
%
% NB: windinglayout is based on a mex interface to the C++ star-of-slots
% based winding layout calculator 'Koil' developed at the Electric Drives
% Laboratory (EDLab) of University of Padova in collaboration with prof.
% Nicola Bianchi.
%
%

    % call the underlying mex winding layout calculator
    coillayout = mexmPhaseWL (NPhases, Qs, Poles/2, SL);
    
    % from the coil layout, create the slot fill scheme, i.e. how the
    % phases are distributed around the slots. This is done by using the
    % contents of coillayout as indices into the phaselayout matrix
    if SL
        
        % in a single layered case, there is one column with the same
        % number of rows as slots (there is the same number of coils as
        % slots)
        phaselayout = nan * ones (Qs, 1);
        
        for phasen = 1:NPhases
            % for each phase the slots in the the pair of columns
            % representing the layout for the phase in coillayout are the
            % positively wound slots and negatively wound slots
            % respectively.
            
            phaselayout( coillayout(:,2*(phasen-1)+1) ) = phasen;
            
            phaselayout( coillayout(:,2*(phasen-1)+2) ) = -phasen;
        
        end
        
    else
        % in a double layered case, thereare two columns with the same
        % number of rows as slots (there are twice as many coils as
        % slots)
        phaselayout = nan * ones (Qs, 2);
        
        for phasen = 1:NPhases
            % for each phase the slots in the the pair of columns
            % representing the layout for the phase in coillayout are the
            % positively wound slot layers and negatively wound slots layer
            % respectively.
            
            phaselayout( coillayout(:,2*(phasen-1)+1),1 ) = phasen;
            
            phaselayout( coillayout(:,2*(phasen-1)+2),2 ) = -phasen;
        
        end
        
        
    end
    
end