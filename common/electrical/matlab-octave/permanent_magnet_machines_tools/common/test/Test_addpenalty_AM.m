
clear ptestdesign ptestsimoptions

ptestdesign.testquant = 100;
ptestsimoptions.max_testquant = 90;
type = 'upper'
quantityfname = 'testquant';
ptestdesign.BaseScore = 1;

[score, ptestdesign, ptestsimoptions] = addpenalty_AM(ptestdesign, ptestsimoptions, type, quantityfname);

dispstruct(ptestdesign, 20);
dispstruct(ptestsimoptions, 20);


%% lower +ve +ve test
clear ptestdesign ptestsimoptions

ptestdesign.testquant = 10;
ptestsimoptions.min_testquant = 30;
type = 'lower'
quantityfname = 'testquant';
ptestdesign.BaseScore = 1;

[score, ptestdesign, ptestsimoptions] = addpenalty_AM(ptestdesign, ptestsimoptions, type, quantityfname, score);

dispstruct(ptestdesign, 20);
dispstruct(ptestsimoptions, 20);

%% lower -ve -ve test
clear ptestdesign ptestsimoptions

ptestdesign.testquant = -30;
ptestsimoptions.min_testquant = -10;
type = 'lower'
quantityfname = 'testquant';
ptestdesign.BaseScore = 1;

[score, ptestdesign, ptestsimoptions] = addpenalty_AM(ptestdesign, ptestsimoptions, type, quantityfname, score);

dispstruct(ptestdesign, 20);
dispstruct(ptestsimoptions, 20);

%% lower +ve -ve test
clear ptestdesign ptestsimoptions

ptestdesign.testquant = -20;
ptestsimoptions.min_testquant = 10;
type = 'lower'
quantityfname = 'testquant';
ptestdesign.BaseScore = 1;

[score, ptestdesign, ptestsimoptions] = addpenalty_AM(ptestdesign, ptestsimoptions, type, quantityfname, score);

dispstruct(ptestdesign, 20);
dispstruct(ptestsimoptions, 20);

%% target higher
clear ptestdesign ptestsimoptions

ptestdesign.testquant = 100;
ptestsimoptions.target_testquant = 90;
type = 'target'
quantityfname = 'testquant';
ptestdesign.BaseScore = 1;

[score, ptestdesign, ptestsimoptions] = addpenalty_AM(ptestdesign, ptestsimoptions, type, quantityfname, score);

dispstruct(ptestdesign, 20);
dispstruct(ptestsimoptions, 20);

%% target lower
clear ptestdesign ptestsimoptions

ptestdesign.testquant = 80;
ptestsimoptions.target_testquant = 90;
type = 'target'
quantityfname = 'testquant';
ptestdesign.BaseScore = 1;

[score, ptestdesign, ptestsimoptions] = addpenalty_AM(ptestdesign, ptestsimoptions, type, quantityfname, score);

dispstruct(ptestdesign, 20);
dispstruct(ptestsimoptions, 20);


