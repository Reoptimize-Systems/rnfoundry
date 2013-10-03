%% Test_designAndEvaluate_ACTIAM

%% Exactly 1/3 pole coils
clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();


%% Less than 1/3 pole coils

clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/4;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();

%% Exactly 1/3 pole coils mode = 0;
clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

design.mode = 0;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();

%% Less than 1/3 pole coils mode = 0

clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/4;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

design.mode = 0;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();

%% Exactly 1/3 pole coils mode = 1;
clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

design.mode = 1;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();

%% Less than 1/3 pole coils mode = 1

clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/4;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

design.mode = 1;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();

%% Exactly 1/3 pole coils mode = 2;
clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

design.mode = 2;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();

%% Less than 1/3 pole coils mode = 2

clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/4;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

design.mode = 2;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();

%% Exactly 1/3 pole coils mode = 3;
clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

design.mode = 3;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();

%% Less than 1/3 pole coils mode = 3

clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/4;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.fillfactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

design.mode = 3;

FemmProblem = femmprob_ACTIAM(design);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();



%% Exactly 1/3 pole coils
clear design
design.phases = 3;         % Number of phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.CoilTurns = 100;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RgVRc = 10; % Ratio of machine resistance to grid resistance

design.poles = [1 1];

design = ratios2dimensions_ACTIAM(design);

design.FEMMTol = 1e-5;

FemmProblem = femmprob_ACTIAM(design, 'CoilCurrents', [1,10,15]);

openprobleminfemm_mfemm(FemmProblem);

main_maximize();

