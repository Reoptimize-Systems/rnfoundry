The 'stlTools' toolbox is a collection of functions, samples and demos
to illustrate how to deal with STL files. Some of them are
contributions published in Matlab Central and properly referenced
here. 

This toolbox uses +package folders to create an 'stl' namespace, so
functions are accessed using dot notation like stl.read or stl.write. The 
following functions are available:

getFormat:  identifies the format of the STL file and returns 'binary' or 
  'ascii'. This file is inspired in the 'READ-stl' file written and 
  published by Adam H. Aitkenhead 
  (http://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation). 
  Copyright (c) 2013, Adam H. Aitkenhead.

readAscii:  reads an STL file written in ascii format. This file is 
  inspired in the 'READ-stl' file written and published by Adam H. 
  Aitkenhead (http://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation). 
  Copyright (c) 2013, Adam H. Aitkenhead

readBinary:  reads an STL file written in binary format. This file is 
  inspired in the 'READ-stl' file written and published by Adam H. 
  Aitkenhead (http://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation). 
  Copyright (c) 2013, Adam H. Aitkenhead

read:  uses 'stlGetFormat', 'stlReadAscii' and 'stlReadBinary' to make 
  STL reading independent of the format of the file

write:  writes an STL file in 'ascii' or 'binary' formats. This is 
  written and published by Sven Holcombe 
  (http://www.mathworks.com/matlabcentral/fileexchange/20922-stlwrite-filename--varargin-). 
  Copyright (c) 2012, Grant Lohsen. Copyright (c) 2015, Sven Holcombe.

slimVerts:  finds and removes duplicated vertices. This function is written 
  and published by Francis Esmonde-White as PATCHSLIM 
  (http://www.mathworks.com/matlabcentral/fileexchange/29986-patch-slim--patchslim-m-). 
  Copyright (c) 2011, Francis Esmonde-White.

getVerts:  returns a list of vertices that are 'opened' or 'closed' 
  depending on the 'mode' input parameter. An 'open' vertice is the one 
  that defines an open side. An open side is the one that only takes part 
  of one triangle

delVerts:  removes a list of vertices from STL files 

addVerts:  adds the new vertices from a list (and consequently, new faces) 
  to a STL object

plot:  is an easy way to plot the STL data

In addition a demo is available in the stl.demo script that shows how to 
use the functions, it makes use of two example stl files:

femur_binary:	is an ascii STL sample used in 'stlDemo'. It is published 
  by Eric Johnson (http://www.mathworks.com/matlabcentral/fileexchange/22409-stl-file-reader). 
  Copyright (c) 2011, Eric Johnson.

sphere_ascii:  is a binary STL sample

