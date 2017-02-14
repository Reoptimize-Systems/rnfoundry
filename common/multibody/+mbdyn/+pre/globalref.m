classdef globalref < handle
   
   properties (GetAccess = public, SetAccess = private )
       
       pos = [0; 0; 0];
       orientm = mbdyn.pre.orientmat ('orientation', eye (3));
       v = [0; 0; 0];
       omega = [0; 0; 0];
       
   end

end