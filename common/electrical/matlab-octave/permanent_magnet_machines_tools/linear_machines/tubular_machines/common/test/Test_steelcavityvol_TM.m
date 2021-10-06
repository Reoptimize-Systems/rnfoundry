
design.WmVWp = 0.5;
design.WpVRm = 0.5;
design.RoVRm = 1.2;
design.RaVRo = 1.05;
design.RsoVRm = 0.2;
design.Rm = 0.1;
design.WcVWp = 1/3;
design.RiVRm = 1.01;
design.RsiVRm = 0.01;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;
design.mode = 1;

%%
design = ratios2dimensions_ACTIAM(design);
%%
if ~isfemmopen, openfemm; end
RunStructFEMMSimNew_ACTIAM(design);

%%
vol = steelcavityvol_TM(design)