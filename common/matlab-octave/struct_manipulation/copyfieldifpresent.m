function S2 = copyfieldifpresent (fnames, S1, S2)
% copies a specified filed from S1 to S2 if it is exists in S1
%
% Syntax
%
% S2 = copyfieldifpresent (fnames, S1)
% S2 = copyfieldifpresent (fnames, S1, S2)
%
% Descrption
%
% copyfieldifpresent copies a field in an input structure to a second
% structure, if the field exists in the first structure. If the field is
% not present there is no change to either structure. If only one structure
% is provided a new structure is creatd to hold the field copies.
% 
% Inputs
%
%  fnames - field name(s) to be copied from S1 to S2, can be a string, or
%    cell array of strings
%
%  S1 - structure from which fields in fnames will be copied if present
%
%  S2 - (optional) structure into which the field will be copied, if not
%    supplied, a new structure will be created into which the copy will be
%    added.
%
% Output
%
%  S2 - structure conating fields copied from S1
%
% See also: setfieldifabsent, isfield
%

    if nargin < 3
        S2 = struct ();
    end
    
    if ischar (fnames)
        fnames = {fnames};
    elseif ~iscellstr (fnames)
        error ('fnames ust be a string, or cell array of strings containing field names');
    end
    
    for ind = 1:numel(fnames)
        if isfield (S1, fnames{ind})
            S2.(fnames{ind}) = S1.(fnames{ind});
        end
    end

end