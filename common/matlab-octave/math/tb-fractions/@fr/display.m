function display(FRAC)
% fr/display: Display a fraction object, calls disp
%
%  See also: disp

% Author: Ben Petschel 25/7/09
%
% Version history:
%   25/7/09 - first release (using vpi/display as a template)


name = inputname(1);
if ~isempty(name)
  display([name,' ='])
else
  display('ans =')
end
disp(FRAC)
