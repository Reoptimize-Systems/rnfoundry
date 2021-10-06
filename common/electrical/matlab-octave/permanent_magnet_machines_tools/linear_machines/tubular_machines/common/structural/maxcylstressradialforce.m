function MaxSigma2 = maxcylstressradialforce(Fi, Fo, cVars)
% maxcylstressradialforce: determines the maximum stresses in the body of a
% cylindrical coil with a linearly varying radial body force from Fi to Fo
%
% The function determines the maximum strain on the coil using formulas
% from roark's formulas for stress and strain. If the maximum strain
% exceeds the strength of the material in the coil the coil will collapse
% or buckle.
%
% To perform this calculation we will superimpose two radial body forces on
% the coil to approximate the forces. the first is a force which is uniform
% throughout the depth of the coil. This force will be the force at the
% outer turn of the coil. The second force is a linearly varying force from
% the inner to the outer part of the coil. This force will vary from the
% forces at the inner turn minus the force at the outer turn (the uniform
% force used) to zero. 
%
% The superposition of these two sets of strains will give us the total
% maximum strain in the coil. If this exceeds the coil strength the coil is
% not structurally sound.
%
% Syntax:
%
%   maxcylstressradialforce(Fi, Fo, )
%
% Input:
%
%   Fi - Force at the inner turns of the coil
%
%   Fo - Force at the outer turns of the coil
%
%   cVars - (n x 4) matrix:-
%           Col 1. a, the outer radius
%           Col 2. b, the inner radius
%           Col 3. v, Poisson's ratio for the material  

    Fvar = Fi - Fo;
    
    MaxSigma2 = Table32r1eMaxSigma2([repmat(Fo,size(cVars,1)) cVars]) + ...
                Table32r1fMaxSigma2([repmat(Fvar,size(cVars,1)) cVars]);

end