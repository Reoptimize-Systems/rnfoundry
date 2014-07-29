% check_flux_linkage_sampling

% setup design
design.g = 5/1000; 

design.tc = 0.02744;

design.hm = 0.3182;
design.Rmm = 1.005;
design.Rmo = design.Rmm + design.hm/2;
design.Rmi = design.Rmm - design.hm/2;  
Npoles = 28;


design.taupm = 0.226;
design.taumm = 0.15;
design.tauco = 4 * design.taupm / 3;
design.Wc = 0.112;
design.tauci = design.tauco - 2 * design.Wc;
design.Dc = design.taumm / 1000;
design.CoilFillFactor = 0.8;
design.Hc = design.tc;
design.CoilTurns =  75;
design.Dc = 6.04 / 1000;
design.tm = 0.15 * design.taumm;
design.tbi = 0.045;

design.NCoilsPerPhase = Npoles;
design.RlVRp = 0.1;
design.LgVLc = 0;
design.Phases = 1;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.HcMag = mgoe2hc(br2mgoe(1.32));

design.MagSimMaterials.Magnet = 'NdFeB 32 MGOe';

design.MagSimMaterials.FieldIron = '1117 Steel';
design.MagSimMaterials.CoilWinding = '36 AWG';


pos = linspace(0.5, 1.5, 15);
circprops = zeros(numel(pos), 3);
% for now, with femm
openfemm;
main_minimize;

for i = 1:numel(pos)
    
    design.Hc = design.tc;
    design.Wc = (design.tauco - design.tauci) / 2;

    [design, simoptions] = simfun_ROTARY(design, simoptions);
    
    design.gap = 2*design.g + design.tc;
    
    % do inductance sim 
%     design.CoilResistance, design.CoilInductance
    
    % Draw the main sim to extract the vector potential
    [design.FemmProblem, design.outermagsep] = corelessfemmprob_torus(design, ...
                                            'NStages', design.NStages, ...
                                            'DrawCoils', true, ...
                                            'Position', pos(i)*design.taupm);
    
    % make an appropriate name for the .fem file from the base name created
    % by simfun_ROTARY
%     femfilename = [simoptions.filenamebase, '.fem'];

    if isfield(design, 'HcMag')
        design.FemmProblem.Materials(2).H_c = design.HcMag;
%         design.FemmProblem.Materials(2).Name = regexprep(design.FemmProblem.Materials(2).Name, '\d+', sprintf('%3.0d', hc2mgoe(design.HcMag)));
    end
    
    femfilename = 'temp_simfun_TORUS_CORELESS.fem';
    
    % write the fem file to disk
    writefemmfile(femfilename, design.FemmProblem);
    
    opendocument(femfilename);
    mi_analyse(1);
    mi_loadsolution;
    
    circprops(i,:) = mo_getcircuitproperties('A');
    
    mo_close;
    mi_close;
    
end