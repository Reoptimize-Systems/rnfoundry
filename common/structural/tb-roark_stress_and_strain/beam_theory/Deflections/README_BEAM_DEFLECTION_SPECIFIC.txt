The following is some information on adding new beam deflection cases to 
this folder.

You add beam deflection functions as per the instructions for all files. 
However, If you add a new function, please add a reference to it in the 
BeamDeflection function by adding its method in the switch case statement 
in this function.
 
This function takes a string (beamMethod) which tells it which method to 
use in calculating the deflection of the beam in question. This function 
is called by BeamDeflectionSuper which superimposes cases, so adding a 
refernce to your new function will give you acess to these higher functions, 
and allow you to superimpose multiple cases and methods. 

The case string will be the table and row for the method you are adding 
from Roark, in the format: TableNum.rowIdentifier

E.g. the function:

	Table3r2e_Def

Is used if the following string is passed in beamMethod:

	3.2e

It is not necessary to identify that this is a deflection formula as only 
references to deflection formulas should be present in the BeamDeflection 
function. How to add these should be obvious with a look at beamDeflection.