This folder contains functions for the calculation of mechanical stress and 
strain based on the formulas presented in 'Roark's Stress and Strain 6th 
Edition'. 

The basic formulas are grouped according to their type and output. E.g. 
functions relating to beams are in the beam_theory folder. Within this folder, 
functions for calculating the second moment of area of a given beam cross-
section are found in the corresponding folder, while functions for calculating
the deflection of beams are in a separate folder.

Anyone is free to add to the function library, but is requested to stick to 
some conventions detailed below. Functions not following these conventions 
will certainly not be corrected, and may well be deleted without notice or 
any form of consultation. There may be further specific advice within each 
folder on further conventions and reccomendations regarding the functions 
contained withing it. Read these before proceeding.


**--------------------------------------------------------------------------**
		       	     NAMING CONVENTIONS
**--------------------------------------------------------------------------**

FUNCTIONS

To add a function using a formula from a table in Roark you name the function
according to the table number and case/row number for the formula you are 
adding e.g.

If you are adding the formula for calculating the 1st moment of area of a beam with 
a solid square cross section, as described on page 62, you would name this 
function:

	Table1r1_I1

If, however, you are adding the formula for calculating the 2nd moment of area 
for the same beam you would name this function:

	Table1r1_I2

If you are adding the formula for calculating the deflection of a beam with a 
simply supported left and right end undergoing a distributed load, as described 
in Table 3 on page 104 in row 2e you would call this:

	Table3r2e_Def

See existing similar functions for naming conventions. If you are adding the 
formula for something which has no precedent, feel free to make up your own,
but continue to follow this once it is established.

FOLDERS

If you are adding a new folder for a group of functions, do not use spaces in 
the name, use underscores (_) instead as this makes the path easier to reference 
in all operating systems.

**--------------------------------------------------------------------------**
		       	     CODING CONVENTIONS
**--------------------------------------------------------------------------**

First and foremost, all function must have a help section describing the inputs 
and outputs, such as the example shown below for Table3r2e_Def:


% function: Table3r2eI_1
% 
% Calculates the deflection of a beam with its left end simply supported
% and its right end simply supported, as calculated in 'Roark's Formulas
% Stress & Strain 6th edition' in table 3, page 104 row 2e.
%
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%          circular cross-section:
%          Yvars(:,1) - wa, unit load at 'a'
%          Yvars(:,2) - wl, unit load at M_B, the end of the beam
%          Yvars(:,3) - l, length of the beam
%          Yvars(:,4) - a, distance from M_A at which 'wa' is applied 
%
%   E - Young's modulus of the beam material
%
%   I - second moment of inertia of the beam cross-section
%
%   x - row vector of position values at which the deflection is to be calculated 
%
% Output:
%
%   Def - (n x 1) column vector of values of the deflection at the
%   corresponding x position
% 
% Author: Joe Blogs 
% Email: j.blogs@ed.ac.uk


Your code will be more helpful to others if it is also well commented throughout.

Where possible, your code should be vectorised, and follow a sensible and easily 
readible layout. You should also take some time to add some error checking on the 
input. 






