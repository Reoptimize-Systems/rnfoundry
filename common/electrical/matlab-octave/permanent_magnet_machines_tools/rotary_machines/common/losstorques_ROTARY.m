function [TqLtot, TqLiron, TqLeddy] = losstorques_ROTARY(design, simoptions, thetaR, omegaR)
% calculates the forces due to eddy current losses in a rotary machine at a
% given angular velocity
%
% Syntax
%
% [TqLtot, TqLiron, TqLeddy] = losstorques_ROTARY(design, simoptions, thetaR, omegaR)
%
% 

    [Pironloss, Pexteddy] = losspower_AM(design, simoptions, thetaR, omegaR);
    
    TqLiron = power2force(omegaR, Pironloss);
    
    TqLeddy = power2force(omegaR, Pexteddy);
    
    TqLtot = TqLiron + TqLeddy;
    
end