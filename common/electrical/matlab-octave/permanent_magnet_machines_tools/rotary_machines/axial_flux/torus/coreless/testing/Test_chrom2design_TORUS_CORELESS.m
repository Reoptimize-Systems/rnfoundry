% Test_chrom2design_TORUS_CORELESS

% number of Phases
options.Phases = 3;
% number of coils per pole and phase
options.qc = fr(3,4*3);
options.RlVRp = 10;
options.CoilFillFactor = 0.9;
options.BranchFac = 0;
options.ModuleFac = 0;
options.NStages = 1;
options.taucmoVtaucsm = 0.99;

% convert machine ratios to actual dimensions
taucmiVtaucmo = 0.2;
tcVtm = 2;
gVtm = 0.05;
tbiiVtbio = 0.9;
tbioVtm = 0.5;
tmVtaumm = 0.1;
taummVtaupm = 0.85;
RsVRbi = 0.95;
RbiVRmi = 0.75;
RmiVRmo = 0.975;
Rmo = 2;
NBasicWindings = 7;
DcAreaFac = 0.01;
BranchFac = 1;
ModuleFac = 0.8;
NStages = 2;

%     design.taucmiVtaucmo = Chrom(1);
%     design.tcVtm = Chrom(2);
%     design.gVtm = Chrom(3);
%     design.tbiiVtbio = Chrom(4);
%     design.tbioVtm = Chrom(5);
%     design.tmVtaumm = Chrom(6);
%     design.taummVtaupm = Chrom(7);
%     design.RsVRbi = Chrom(8);
%     design.RbiVRmi = Chrom(9);
%     design.RmiVRmo = Chrom(10);
%     design.Rmo = Chrom(11);
%     design.NBasicWindings = Chrom(12);
%     design.DcAreaFac = Chrom(13);
%     design.BranchFac = Chrom(14);
%     design.ModuleFac = Chrom(15);

Chrom = [ taucmiVtaucmo, ...
          tcVtm, ...
          gVtm, ...
          tbiiVtbio, ...
          tbioVtm, ...
          tmVtaumm, ...
          taummVtaupm,...
          RsVRbi, ...
          RbiVRmi, ...
          RmiVRmo, ...
          Rmo, ...
          NBasicWindings, ...
          DcAreaFac, ...
          BranchFac, ...
          ModuleFac ];

[design, simoptions] = chrom2design_TORUS_CORELESS(struct(), Chrom);



