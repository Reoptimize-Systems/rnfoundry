classdef globalref < handle
% provides a reference which always refers to the global coordinate frame
%
% See als: mbdyn.pre.reference
%
   
   properties (GetAccess = public, SetAccess = private )
       
       name = 'global'; % name for the reference (used in plots)
       pos = [0; 0; 0]; % position of the reference in the global frame: [0; 0; 0]
       orientm = mbdyn.pre.orientmat ('orientation', eye (3)); % orientation of the reference in the global frame (3 x3 identity matrix)
       v = [0; 0; 0]; % velocity of the reference in the global frame: [0; 0; 0]
       omega = [0; 0; 0]; % angular velocity of the reference in the global frame :[0; 0; 0]
       
   end

end