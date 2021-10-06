% Test_checkdesign_am

clear design

design.a = 1;

design.b = 2;

design.c = 3;

design.d = 4;

design.aVb = design.a / design.b;

valfields = {'a', 'b', 'c'};

ratiofields = {'aVb', 'cVb', 'dVc'};

basefield = 'b';

design = checkdesign_am(design, valfields, ratiofields, basefield)



