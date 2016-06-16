function [FesVec, Fes, FfricVec, Ffric, unitv] = buoysystemforces_AM(design, simoptions, xT, vT, xBh, xBs)
% determines forces acting on a buoy system other than wave interaction and
% PTO forces, e.g. end stops
%
% Syntax
%
% [FesVec, Fes, unitv] = buoysystemforces_AM(simoptions, xT, vT, xBh, xBs)
%
% 

    % find a unit vector pointing from the hause hole to the buoy
    unitv = [simoptions.BuoySim.tether_length+1e-6+xBh, xBs] / norm([simoptions.BuoySim.tether_length+1e-6+xBh, xBs]);
    
    % calculate the end-stop forces
    [FesVec, Fes] = endstopforces(simoptions, xT, vT, unitv);
    
    % calculate a frictional force
    [FfricVec, Ffric] = buoysystemfriction(vT, design.WECFriction, unitv);

end


function [FesVec, Fes] = endstopforces(simoptions, xT, vT, unitv)
% calcualtes the force transmitted to a heaving buoy via a tether due to
% PTO end-stops
%
% 

    % add end stop forces
    if (xT > simoptions.maxAllowedxT && vT > 0) || (xT < -simoptions.maxAllowedxT && vT < 0)

        Fes = - sign(xT) * simoptions.EndStopks * (abs(xT) - simoptions.maxAllowedxT);
        
    elseif (xT > simoptions.maxAllowedxT && vT <= 0) || (xT < -simoptions.maxAllowedxT && vT >= 0)
        
        % get the end stop force
        Fesmag = simoptions.EndStopks * (abs(xT) - simoptions.maxAllowedxT);
        
        % modify the force so it falls faster if we are moving in the
        % opposite direction, ensuring we are not adding energy by
        % exceeding the end stop force
        Fvelmod = min(Fesmag, abs(vT) * (Fesmag / 1));
        
        Fesmag = Fesmag - Fvelmod;
        
        % add it to the total force
        Fes = - sign(xT) * Fesmag;
       
    else
        Fes = 0;
    end
    
    FesVec = Fes * unitv;
    
end
    
function [FfricVec, Ffric] = buoysystemfriction(v, frictionmag, unitv)
% calculates frictional forces on the buoy system
%
% Syntax
%
% [FfricVec, Ffric] = buoysystemfriction(v, frictionmag, unitv)
%
% Description
%
% The force calculated by buoysystemfriction is scaled if the velocity is
% less than a threshold velocity
%

    vthresh = 0.01;
    
    vabs = abs(v);

    % the function v^2 / vthresh^2 is a quadratic which rises from zero to
    % 1 between v = 0 and v = vthresh. 
    frictionmag(vabs < vthresh) = frictionmag(vabs < vthresh) .* vabs(vabs < vthresh).^2 ./ (vthresh^2);
    
    Ffric = -sign(v) * frictionmag;
    
    FfricVec = Ffric * unitv;

end
