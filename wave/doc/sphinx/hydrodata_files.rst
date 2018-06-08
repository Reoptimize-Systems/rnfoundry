Hydrodynamic Data Files
***********************

The case directory should contain a subdirectory 'hydroData'. The 
hydroData directory may contain various files which can be processed 
or used directly by the |TNShort|. The possible files are described 
in the following sections.
   
Nemoh
=====

Note that the Nemoh file organisation is somewhat complex. For this 
reason, the Nemoh preprocessor was developed which handles producing 
and locating all of these files as required, without intervention by 
the user. See :ref:`nemoh-interface` for more information.

+------------+---------------------------------------------------------+--------------------------------------------------------------------+
| Component  | name / extension                                        | Description                                                        |
+------------+---------------------------------------------------------+--------------------------------------------------------------------+
| Mesher     |  mesh_input_for_<body_name>                             | | directory containing Nemoh mesh processor input files for a      | 
|            |                                                         | | body, actually this directory could have any name, but the Nemoh |
|            |                                                         | | preprocessor uses this convention, where <body_name> is replaced |
|            |                                                         | | by the body name, e.g. body_id_0                                 |
+            +---------------------------------------------------------+--------------------------------------------------------------------+
|            |  mesh_input_for_<body_name>/ID.dat                      | | in mesh_input_for_<body_name> directory, needed for Nemoh mesh   |
|            |                                                         | | processor/refiner                                                |
+            +---------------------------------------------------------+--------------------------------------------------------------------+
|            |  mesh_input_for_<body_name>/Mesh.cal                    | | in mesh_input_for_<body_name> directory, contains info about the |
|            |                                                         | | Nemoh course mesh input file                                     |
+            +---------------------------------------------------------+--------------------------------------------------------------------+
|            |  mesh_input_for_body_<body_name>/mesh                   | | subdirectory of the mesh_input_for_<body_name> subdirectory,     |
|            |                                                         | | contains the Nemoh course mesh input file, and output of the     |
|            |                                                         | | mesh preprocessor                                                |
+            +---------------------------------------------------------+--------------------------------------------------------------------+
|            |  mesh_input_for_body_<body_name>/mesh/<body_name>       | | file with no extention containing the Nemoh course mesh input    |
|            |                                                         | | file                                                             |
+------------+---------------------------------------------------------+--------------------------------------------------------------------+
| Solver     |  ID.dat                                                 | needed for Nemoh solver (different from ID.dat above for the mesh) |
+            +---------------------------------------------------------+--------------------------------------------------------------------+
|            |  mesh                                                   | | copy of mesh_input_for_body_*/mesh (see above), for technical    |
|            |                                                         | | reasons this is required to be in the top level of hydroData     |
|            |                                                         | | once mesh preprocessing is complete to input to the solver       |
+            +---------------------------------------------------------+--------------------------------------------------------------------+
|            |  mesh                                                   | | copy of mesh_input_for_body_*/mesh (see above), for technical    |
|            |                                                         | | reasons this is needed in the top level of hydroData once mesh   |
|            |                                                         | | preprocessing is complete                                        |
+------------+---------------------------------------------------------+--------------------------------------------------------------------+



WAMIT
=====

+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                              | needed for Nemoh solver (different from ID.dat above for the mesh)                                                                                                            |
+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


AQUA
====

+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                              | needed for Nemoh solver (different from ID.dat above for the mesh)                                                                                                            |
+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


Processed BEM Output Files
==========================

This section refers to files which have been created by the BEMIO 
tools and are ready to be loaded directly by the ``wsim.hydroBody`` 
class.

+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  .h5                         |                                                                                                                                                                               |
+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|  .mat                        |                                                                                                                                                                               |
+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

