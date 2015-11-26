function tf = doproceed (msg, default)
% Interactively ask a user whether to proceed or stop on the command line
%
% Syntax
%
% tf = ui.doproceed (msg, default)
%
% Input
%
%   msg - string containing the message to be given to the user when
%     prompting for a reply. This will be appended with the either the
%     string ' [Y/n]', or ' [y/N]' depending on the desire default value.
%
%   default - optional boolean value determining the default value of the
%     user's response. This value is applied if the user hits enter without
%     supplying any input. If default is true, the string ' [Y/n]' is
%     appended to the message string, and the default answer is true. If
%     default is false, the string ' [y/N]' is appended to the message
%     string and the default response is false.
%
%     If not supplied, false is used as the default value.
%
% Output
%
%   tf - boolean value representing the user's chosen response, true for
%     yes and false for no.
%
%
   
    if nargin < 2
        default = false;
    end
    
    if default
        endstr = ' [Y/n]: ';
    else
        endstr = ' [y/N]: ';
    end
    
    if ~ischar (msg)
        error ('msg must be a string');
    end
    
    str = input ([msg, endstr],'s');
    
    while ischar (str)
        
        switch str

            case 'Y'
                str = 1;
                tf = true;

            case 'y'
                str = 1;
                tf = true;

            case 'N'
                str = 1;
                tf = false;

            case 'n'
                str = 1;
                tf = false;

            case ''
                str = 1;
                tf = default;

            otherwise
                % repeat the prompt to the user, they must answer!!!
                str = input ([msg, endstr],'s');
        end
    
    end
    
end