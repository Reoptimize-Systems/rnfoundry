function PMSM_I = TranslatorMomentOfArea_DS_PMSM(IVars)
% TranslatorMomentOfArea_DS_PMSM
%
% A function for calculating the second moment of area of a single pole of
% the double sided linear permanent magnet machine.
%
% Input:
%
%   IVars - (1 x 4) matrix of values for calculating the second moment of
%           area of the teeth and support/central translator section
%
%           IVars(1,1): dbi, half thickness of central section (back 
%                       iron thickness of armature)
%           IVars(1,2): Taup, the pole width
%           IVars(1,3): hs, the tooth height
%           IVars(1,4): bt, the tooth width
%
% Output:
%
%   PMSM_I - second moment of area of single pole

    
    % IVars(1,1) is the back iron thickness, i.e. half the thickness of the
    % central part of the translator. Therefore we double this to get the
    % total thickness/height of the section
    dbi = IVars(:,1);
    dc = dbi * 2;
    bc = IVars(:,2);
    dt = IVars(:,3);
    bt = IVars(:,4);
    
    % Calculate second moment of area of the translator central section and
    % a single tooth
    Ic = MomentOfInertiaY1([dc, bc], '1.2');
    It = MomentOfInertiaY1([dt, bt], '1.2');
    
    % Teeth neutral axis
    Y2 = dbi + (dt .* 0.5);
    
    % Total second moment of area = I1 + 6(I2 + A2*Y2^2)
    PMSM_I = Ic + 6 .* (It + ( (dt .* bt) .* (Y2^2) )); 
    
end