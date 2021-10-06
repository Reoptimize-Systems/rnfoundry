function varargout = listsearch(funname, varlist, funargs)
% listsearch: finds the member of a sorted list who just satisfies some
% specified criteria.
%
% Syntax:
%
%   r = listsearch(funname, varlist)
%   r = listsearch(..., funargs)
%   [r,V] = listsearch(...)
%
% Description:
%
% This function takes a list of variables sorted according to some order
% (e.g. price of the list member) and finds the lowest member of the list
% who satisfies some other criteria specified by the function in funname.
% The member is found using a simple bisection algorithm.
%
% N.B. listsearch assumes the output of the function specified in funname
% is a step function with a single 'step-up' point. I.e. if every member of
% a list was evaluated for success the output would look something like
% this:
%
%           [ 0; 0; 0; 0; 0; 1; 1; 1; 1; 1; 1 ]
%
% and not this:
%
%           [ 0; 0; 1; 0; 0; 1; 1; 0; 1; 1; 1 ]
%
% so once the lowest successful member is found, all previous members are
% assumed to be unsuccessful, and all subsequent members assumed
% successful.
%
%
% Example of listsearch use:
%
% You have a number of (equally sized) oranges you wish to pack in a box
% for transportation. You need to buy a box in which all the oranges fit,
% but don't want to waste money buying a box too big for the job. When you
% go to buy a box you find there are 20 standard sizes available, each at
% different prices. The prices are directly proportional to the volume of
% the boxes.
%
% To use listsearch to determine which is the best to buy, you must first
% write a function which accepts the dimensions of a box as input and
% outputs true if the oranges all fit in or false if they will not. You
% would then put the dimensions of each box into a matrix in which each row
% described the dimensions of the box being evaluated. You would then sort
% this list of box dimensions according to the price (e.g. using the MATLAB
% function 'sortrows'). 
%
% Having done this simply pass in the name of your evaluation function, the
% sorted list, and any other required variables and listsearch will locate
% the lowest cost box which satisfies the criteria.
%
% Input:
%
%   funname - name or function handle of the function which will evaluate
%             the success of each member of varlist. This function must
%             return true (1) or false (0) depending on the success of the
%             list member
%
%   varlist - (n x p) matrix each row of which is a vector of input
%             variables for the function specified by funname. The rows
%             will be sorted in ascending order according to the criteria
%             you wish to base your choice on, e.g. the price of the box in
%             the example above. Alternatively, this could be a column
%             vector of structures containing the appropriate information.
%
%   funargs - optional cell array of additional arguments to pass to the
%             evaluation function for every evaluation
%
% Output:
%
%   r - the row number of the successful list member earliest in list, this
%       will be -1 if no list members satisfy the criteria
%
%   V - (optional) the contents of the winning list member, i.e.
%       varlist(r,:). If no winner is found, returns an empty matrix.
%

    if isempty(varlist)  
        error('varlist is empty matrix')
    end
    
    if nargin < 3
        funargs = {};
    end
    
    % evaluate the first member of the list, we can stop at this point if
    % it is successful, i.e. we will skip the while loop and return 1
    fargs = [{varlist(1,:)}, funargs];
    success1 = feval(funname, fargs{:});
    
    if success1
        varargout{1} = 1;
        varargout{2} = varlist(1,:);
        return
    end
    
    % evaluate the last member of the list as if it is not successful, we
    % will waste time trying to search for one which works. If this is the
    % case we will return -1
    max_bound2 = size(varlist, 1);
    
    fargs = [{varlist(max_bound2,:)}, funargs];
    
    if ~feval(funname, fargs{:})
        varargout{1} = -1;
        varargout{2} = [];
        return
    end
    
    % we will ensure we do not get locked in a endless loop due to
    % incorrectly sorted data or some other problem
    maxloops = max_bound2 + 2;
    loopcount = 0;
    
    % If we enter the loop we know that the first member of the list
    % was unsuccessful, but the last member is successful
    bound2 = ceil(max_bound2 / 2);
    bound1 = 1;
    
    fargs = [{varlist(bound2,:)}, funargs];
    success2 = feval(funname, fargs{:});
    
    while loopcount <= maxloops
        % we test the half way point between the two bounds in each loop,
        % until the bounds are next to each other
        if success2 && ~success1
            % if the two bounds are next to each other, we have found the
            % right member, otherwise, we move the upper bound down again
             if bound2 - bound1 == 1 || bound1 == bound2
                 varargout{1} = bound2;
                 varargout{2} = varlist(bound2,:);
                 return
             else
                 % Still searching, so lower the new upper bound of
                 % successful members
                 max_bound2 = bound2;
                 % Now choose a test point half way between the lower bound
                 % and maximum upper bound
                 bound2 = bound1 + ceil((max_bound2 - bound1) / 2);
                 % Test this point for success
                 fargs = [{varlist(bound2,:)}, funargs];
                 success2 = feval(funname, fargs{:});
             end
             
        elseif ~success2 && ~success1
            % move the lower bound up to the new lowest unsuccessful member
            bound1 = bound2;
            % set the upper bound to the member half way between the old
            % upper bound and the new lower bound
            bound2 = bound1 + ceil((max_bound2 - bound1)/2);
            % evaluate the new upper bound for success
            fargs = [{varlist(bound2,:)}, funargs];
            success2 = feval(funname, fargs{:});
        end
        
        loopcount = loopcount + 1;

    end
    
    if loopcount > maxloops
        error('Crozier:InfiniteLoop','More iterations than list members have been performed,\ncheck input data is correctly sorted and evaluation function\nis operating correctly')
    end

end