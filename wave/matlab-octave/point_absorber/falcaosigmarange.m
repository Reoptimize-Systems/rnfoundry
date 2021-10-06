function range = falcaosigmarange()
% returns the range of angular frequencies recommended by Falcoe for sea
% wave simulation using a spectrum of frequencies
%
% Syntax
%
% range = falcaosigmarange()
%
% Description
%
% falcaosigmarange returns a two element vector containing the minimum and
% maximum angular frequency recommended by Falcao for sea wave spectrum
% simulation. The paper in question is Phase control through load control
% of oscillating-body wave energy converters with hydraulic PTO system ,
% Ocean Engineering 2008.
%

    range = [0.1*0.6^0.5, 0.1*(0.6)^0.5 + 2.24];
    
end