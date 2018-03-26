classdef globalref < handle
% provides a reference which always refers to the global coordinate frame
%
% See als: mbdyn.pre.reference
%
   
   properties (GetAccess = public, SetAccess = private )
       
       name = 'global';
       pos = [0; 0; 0];
       orientm = mbdyn.pre.orientmat ('orientation', eye (3));
       v = [0; 0; 0];
       omega = [0; 0; 0];
       
   end

end