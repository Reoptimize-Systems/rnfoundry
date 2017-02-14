classdef reference < handle
    
   properties (GetAccess = public, SetAccess = private)
       
       positionParent;
       orientParent;
       velocityParent;
       omegaParent;
       
       dpos;
       dorientm;
       dv;
       domega;
       
   end
   
   properties (Dependent)
       
       pos;
       orientm;
       v;
       omega;
       
   end
   
   methods
       
       function this = reference (dpos, dorientm, dv, domega, varargin)
           
           options.Parent = mbdyn.pre.globalref ();
           options.PositionParent = [];
           options.OrientParent = [];
           options.VelParent = [];
           options.OmegaParent = [];
           
           options = parse_pv_pairs (options, varargin);
           
           if ~ ( isa (options.Parent, 'mbdyn.pre.reference') ...
                   || isa (options.Parent, 'mbdyn.pre.globalref') )
               
               error ('RENEWNET:mbdyn:badrefparent', ...
                   'The parent must be another reference object or the global object');
               
           end
           
           if isempty (options.PositionParent)
               options.PositionParent = options.Parent;
           else
               if ~ ( isa (options.PositionParent, 'mbdyn.pre.reference') ...
                       || isa (options.PositionParent, 'mbdyn.pre.globalref') )
                   
                   error ('RENEWNET:mbdyn:badrefparent', ...
                       'The PositionParent must be another reference object or the global object');
                   
               end
           end
           
           if isempty (options.OrientParent)
               options.OrientParent = options.Parent;
           else
               if ~ ( isa (options.OrientParent, 'mbdyn.pre.reference') ...
                       || isa (options.OrientParent, 'mbdyn.pre.globalref') )
                   
                   error ('RENEWNET:mbdyn:badrefparent', ...
                       'The OrientParent must be another reference object or the global object');
                   
               end
           end
           
           if isempty (options.VelParent)
               options.VelParent = options.Parent;
           else
               if ~ ( isa (options.VelParent, 'mbdyn.pre.reference') ...
                       || isa (options.VelParent, 'mbdyn.pre.globalref') )
                   
                   error ('RENEWNET:mbdyn:badrefparent', ...
                       'The OrientParent must be another reference object or the global object');
                   
               end    
           end
           
           if isempty (options.OmegaParent)
               options.OmegaParent = options.Parent;
           else
               if ~ ( isa (options.OmegaParent, 'mbdyn.pre.reference') ...
                       || isa (options.OmegaParent, 'mbdyn.pre.globalref') )
                   
                   error ('RENEWNET:mbdyn:badrefparent', ...
                       'The OrientParent must be another reference object or the global object');
                   
               end
           end
           
           this.positionParent = options.PositionParent;
           this.orientParent = options.OrientParent;
           this.velocityParent = options.VelParent;
           this.omegaParent = options.OmegaParent;
           
           if isempty (dpos)
               dpos = [ 0; 0; 0 ];
           end
           
           if isempty (dorientm)
               this.dorientm = mbdyn.pre.orientmat ('orientation', eye ());
           elseif ~isa (dorientm, 'mbdyn.pre.orientmat')
               error ('RENEWNET:mbdyn:badreforientation', ...
                  'dorientm should be a mbdyn.pre.orientmat  object or empty' );
           else
               this.dorientm = dorientm;
           end
           
           if isempty (dv)
               dv = [ 0; 0; 0 ];
           end
           
           if isempty (domega)
               domega = [ 0; 0; 0 ];
           end
           
           check.multicheck ( @(x) (isnumeric(x) && (size (x,1) == 3) && (size (x,2) == 1) ), ...
                              'dpos, dorientm, dv, domega must all be numeric column vectors 3 elements (or empty)', ...
                              'RENEWNET:mbdyn:badrefvalues', ...
                                dpos, dv, domega );
                            
           this.dpos = dpos;
           this.dv = dv;
           this.domega = domega; 
           
       end
       
       function value = get.pos (this)
           value =  this.get_pos ();
       end
       
       function value = get.orientm (this)
           value = this.get_orientm ();
       end
       
       function value = get.v (this)
           value = this.get_v ();
       end
       
       function value = get.omega (this)
           value = this.get_omega ();
       end
       
   end
   
   methods (Access = private)
       % need this complication as you cannot access dependent properties
       % of a class property in Maltab, i.e. the x, orientm, v and omega
       % properties of the parent reference from this class
       
       function value = get_pos (this)
           % return cartesian position in global frame
           
           % rotate translation vector in parent's orientation
           newdpos = this.dpos.' * this.orientParent.orientm.orientationMatrix;
           
           % translate points
           value = this.positionParent.pos + newdpos';
       end
       
       function value = get_orientm (this)
           % get the orientation of this reference relative to parent
           %
           value =  this.dorientm * this.orientParent.orientm;
       end
       
       function value = get_v (this)
           value = this.velocityParent.v + this.dv;
       end
       
       function value = get_omega (this)
           value = this.omegaParent.omega + this.domega;
       end
       
   end
    
end