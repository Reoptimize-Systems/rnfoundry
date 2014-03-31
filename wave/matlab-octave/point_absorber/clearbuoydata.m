function simoptions = clearbuoydata(simoptions)

    simoptions = rmiffield(simoptions, 'buoylibdir');
    
    simoptions = rmiffield(simoptions, 'BuoyParameters');
    
    simoptions = rmiffield(simoptions, 'ExcitationFile');
    
    simoptions = rmiffield(simoptions, 'HeaveFile');
    
    simoptions = rmiffield(simoptions, 'SurgeFile');
    
    simoptions = rmiffield(simoptions, 'HydroCoeffsFile');

end