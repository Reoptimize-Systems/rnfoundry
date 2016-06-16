function [results, design] = splitsystemresfun_ACTM(T, Y, design, simoptions, results)

        design.minArmLength = results.minArmLength;
        
        design.minArmPoles = ceil(design.minArmLength ./ design.Wp);

        design.minArmLength = design.minArmPoles * design.Wp;
    
        design.Irms = sqrt(results.Isquaredsum / results.Isquaredn);
        
        design.Jrms = design.Irms / design.conductorArea;
        
        design.EMFrms = sqrt(results.EMFsquaredsum / results.EMFsquaredn);

        design.Ipeak = results.Ipeak;
        design.maxJ = design.Ipeak / design.conductorArea;

        design.EMFpeak = results.EMFpeak;

        design.TotalGridEnergy = results.TotalGridEnergy;
        design.AverageEnergy = design.TotalGridEnergy ./ (simoptions.ODESim.TimeSpan(2) - simoptions.ODESim.TimeSpan(1));

        design.GridMeanPower = results.GridPowersum / results.GridPowern;

        design.peakPower = results.peakPower;
        
end
