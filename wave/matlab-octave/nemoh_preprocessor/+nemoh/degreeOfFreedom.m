classdef degreeOfFreedom < nemoh.base
    
   properties

       directionVector;
       comment;
       
   end
   
   
   methods
       
       function  self = degreeOfFreedom (direction, varargin)
           
           options.Comment = '';
           
           options = parse_pv_pairs (options, varargin);
           
           self.directionVector = direction;
           
           self.comment = options.Comment;
           
       end
       
       function set.directionVector (self, direction)
           
           assert ( isnumeric (direction) ...
                      && numel (direction) == 3 ...
                      && isreal (direction) ...
                    , 'direction must be a real numeric vector of length 3' );
                
            self.directionVector = direction;
            
       end
       
       function set.comment (self, new_comment)
           
           assert ( ischar (new_comment) || (isstring (new_comment) && isscalar (new_comment)) ...
                    , 'comment must be a character vector or scalar string' );
                
            self.comment = new_comment;
            
       end
       
   end
    
end