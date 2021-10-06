function [pos] = tricenter(Pa,Pb,Pc,type,doplot)
% tricenter calculates and optionally shows the orthocenter, circumcenter,
% barycenter and incenter of a triangle, given their vertex's coordinates
% Pa, Pb and Pc
%
% Syntax
%
% [pos] = tricenter(Pa,Pb,Pc,type,doplot)
%
% Input
%
%   Pa - 3D coordinates of first vertex
%   Pb - 3D coordinates of second vertex
%   Pc - 3D coordinate of thrid vertex
%   
%   type - string determining type of center desired, options are:
%     'incenter', 'barycenter, 'circumcenter', 'orthocenter'. If not
%     supplied, the 'incenter' is returned.
% 
%   doplot - flag determining if the result and triangle should be plotted.
%     If evaluating to true, the plot is produced. Defaults to false if not
%     supplied.
%
% Example: tricenter([0 0.5 0], [1 0 0], [1 1 3], 'orthocenter')
%
% Made by: Ing. Gustavo Morales, University of Carabobo, Venezuela.
% 09/14/09
%
% Modified By: Richard Crozier, University of Edinburgh, UK
% 16 March 2012

    if nargin < 5
        doplot = false;
    end
    
    if nargin < 4
        type = 'incenter';
    end
    
    Pa = Pa(:); Pb = Pb(:); Pc = Pc(:); % Converting to column vectors (if needed)
    AB = Pb - Pa; AC = Pc - Pa; BC = Pc - Pb; % Side vectors
    
    switch type
        case 'incenter'%
            uab = AB./norm(AB); uac = AC./norm(AC); ubc = BC./norm(BC); uba = -uab;
            L1 = uab + uac; L2 = uba + ubc; % directors
            P21 = Pb - Pa;
            P1 = Pa;
        case 'barycenter'
            L1 = (Pb + Pc)/2 -Pa; L2 = (Pa + Pc)/2 - Pb; % directors
            P21 = Pb - Pa;
            P1 = Pa;
        case 'circumcenter'
            N = cross(AC,AB);
            L1 = cross(AB,N); L2 = cross(BC,N); % directors
            P21 = (Pc - Pa)/2;
            P1 = (Pa + Pb)/2;
        case 'orthocenter'
            N = cross(AC,AB);
            L1 = cross(N,BC); L2 = cross(AC,N); % directors
            P21 = Pb - Pa;
            P1 = Pa;
        otherwise
            error('Unknown Center Type');
    end
    ML = [L1 -L2]; % Coefficient Matrix
    lambda = ML\P21;  % Solving the linear system
    pos = P1 + lambda(1)*L1; % Line Equation evaluated at lambda(1)

    if doplot
        X = [Pa(1); Pb(1); Pc(1)]; Y = [Pa(2); Pb(2); Pc(2)];
        if numel(Pa)== 3 % Tridimensional Case
            Z = [Pa(3); Pb(3); Pc(3)];
            patch(X,Y,Z,'b','FaceAlpha',0.5); hold on;

            plot3(pos(1),pos(2),pos(3),'o','Color','r'); view(3);

        else             % Bidimensional case
            patch(X,Y,'b','FaceAlpha',0.5); hold on;
            plot(pos(1),pos(2),'o','Color','r');
        end
        grid on; daspect([1 1 1]);
    end

end
 
% Original FIle Exchange Licence:
%
% Copyright (c) 2009, Gustavo Morales
%
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Universidad de Carabobo nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.
