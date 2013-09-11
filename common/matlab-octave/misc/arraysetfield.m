function s = arraysetfield(s, field, v)
% Sets all the members of an array of structure to have the same value in a
% field
%
%
        
    for i = 1:numel(s)

        s(i).(field) = v;

    end

end
    