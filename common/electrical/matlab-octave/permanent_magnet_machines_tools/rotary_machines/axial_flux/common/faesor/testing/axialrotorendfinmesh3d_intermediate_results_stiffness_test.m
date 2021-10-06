
 E = 207e9;

nu = 0.31;


 [feb,geom,u] = axialrotorfeb(fens, gcells, Rbi, 3 * tbi, tbi, 2*tbi, E, nu);


% Assemble the system matrix
ems = stiffness(feb, geom, u)

%%

 E = 207e9;

nu = 0.31;


 [feb,geom,u] = axialrotorfeb(finfens, fingcells, Rbi, 3 * tbi, tbi, 2*tbi, E, nu);


% Assemble the system matrix
ems = stiffness(feb, geom, u)

%%

 E = 207e9;

nu = 0.31;


 [feb,geom,u] = axialrotorfeb(shaftfens, shaftgcells, Rbi, 3 * tbi, tbi, 2*tbi, E, nu);

 %%

 E = 207e9;

nu = 0.31;


 [feb,geom,u] = axialrotorfeb(suppfens, suppgcells, Rbi, 3 * tbi, tbi, 2*tbi, E, nu);


% Assemble the system matrix
ems = stiffness(feb, geom, u)