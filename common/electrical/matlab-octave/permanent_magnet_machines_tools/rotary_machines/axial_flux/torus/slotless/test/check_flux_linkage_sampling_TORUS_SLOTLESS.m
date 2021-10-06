% check_flux_linkage_sampling

% setup design
design.g = 5/1000; 

design.tc = 0.1;

design.Rmo = 10;
design.Rmi = 9.5;  
Npoles = 400;

design.taupm = (pi * (design.Rmo + design.Rmi)) / Npoles;
design.taumm = 0.85 * design.taupm;
design.tauco = 0.333 * design.taupm;
design.Hc = design.tc;
design.Wc = design.tauco;
design.Dc = design.taumm / 1000;
design.CoilFillFactor = 0.8;
design.Hc = design.tc;
design.Wc = design.tauco;
[design.CoilTurns, design.Dc] = CoilTurns(design.Hc * design.Wc, design.CoilFillFactor, design.Dc);

design.tm = 0.15 * design.taumm;
design.tbi = 0.05;
design.ty = 2 * design.tbi;

design.NPhaseCoils = Npoles;
design.RlVRp = 0.1;
design.LgVLc = 0;
design.Phases = 3;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.MagFEASimMaterials.Magnet = 'NdFeB 32 MGOe';
design.MagFEASimMaterials.FieldBackIron = '1117 Steel';
design.MagFEASimMaterials.ArmatureCoil = '36 AWG';


pos = linspace(0, 1, 15);
circprops = zeros(numel(pos), 3);
% for now, with femm
openfemm;
main_minimize;

for i = 1:numel(pos)
    
    design.Hc = design.tc;
    design.Wc = design.tauco / 2;

    [design, simoptions] = simfun_ROTARY(design, simoptions);
    
    design.gap = design.g + design.tc;
    
    % do inductance sim 
%     design.CoilResistance, design.CoilInductance
    
    % Draw the main sim to extract the vector potential
    [design.FemmProblem, design.outermagsep] = slotlessfemmprob_torus(design, ...
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
    
    femfilename = 'temp_simfun_TORUS_SLOTLESS.fem';
    
    % write the fem file to disk
    writefemmfile(femfilename, design.FemmProblem);
    
    opendocument(femfilename);
    mi_analyse(1);
    mi_loadsolution;
    
    circprops(i,:) = mo_getcircuitproperties('A');
    
    mo_close;
    mi_close;
    
end

plot(circprops(:,3))