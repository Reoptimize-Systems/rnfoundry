classdef reference < handle
    
   properties (GetAccess = public, SetAccess = private)
       
       name;
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
           options.Name = '';
           
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
           
           if isempty (dpos) || strcmp (dpos, 'null')
               dpos = [ 0; 0; 0 ];
           end
           
           if isempty (dorientm) || strcmp (dorientm, 'null')
               this.dorientm = mbdyn.pre.orientmat ('orientation', eye (3));
           elseif ~isa (dorientm, 'mbdyn.pre.orientmat')
               error ('RENEWNET:mbdyn:badreforientation', ...
                  'dorientm should be a mbdyn.pre.orientmat  object or empty' );
           else
               this.dorientm = dorientm;
           end
           
           if isempty (dv) || strcmp (dv, 'null')
               dv = [ 0; 0; 0 ];
           end
           
           if isempty (domega) || strcmp (domega, 'null')
               domega = [ 0; 0; 0 ];
           end
           
           check.multicheck ( @(x) (isnumeric(x) && (size (x,1) == 3) && (size (x,2) == 1) ), ...
                              'dpos, dv, domega must all be numeric column vectors 3 elements (or empty)', ...
                              'RENEWNET:mbdyn:badrefvalues', ...
                                dpos, dv, domega );
                            
           if ischar (options.Name)
               this.name = options.Name;
           else
               error ('Name must be a char array.')
           end
                            
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
       
       function [hax, hquiv] = draw (this, varargin)
           
           options.PlotAxes = [];
           options.Title = true;
           options.DrawGlobal = true;
           options.Scale = 1;
           
           options = parse_pv_pairs (options, varargin);

           [hquiv, hax] = draw ( this.orientm, ...
                                 'PlotAxes', options.PlotAxes, ...
                                 'Title', false, ...
                                 'Offset', this.pos, ...
                                 'DrawGlobal', false, ...
                                 'Scale', options.Scale );
                             
           x = options.Scale * [1;0;0];
           y = options.Scale * [0;1;0];
           z = options.Scale * [0;0;1];
           
           if options.DrawGlobal
               % global frame
               hquiv(4) = vect.plotvec3 (x, [], 'Properties', {'Color', 'r'}, 'PlotAxes', hax);
               hquiv(5) = vect.plotvec3 (y, [], 'Properties', {'Color', 'g'}, 'PlotAxes', hax);
               hquiv(6) = vect.plotvec3 (z, [], 'Properties', {'Color', 'b'}, 'PlotAxes', hax);
           end
           
           if options.Title
               title ('Reference Plot')
           end
           
           axis (hax, 'equal');
           
       end
       
   end
   
   methods (Access = private)
       % need this complication as you cannot access dependent properties
       % of a class property in Maltab, i.e. the x, orientm, v and omega
       % properties of the parent reference from this class
       
       function value = get_pos (this)
           % return absolute cartesian position in global frame
           %
           
           value = this.positionParent.pos ...
               + this.orientParent.orientm.orientationMatrix * this.dpos;

       end
       
       function value = get_orientm (this)
           % get the absolute orientation of this reference in the global
           % frame
           %
           
           value =  this.orientParent.orientm * this.dorientm;
           
       end
       
       function value = get_v (this)
           value = this.velocityParent.v + this.dv;
       end
       
       function value = get_omega (this)
           value = this.omegaParent.omega + this.domega;
       end
       
   end
    
end