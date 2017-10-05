Getting Started
***************

EMST is is based in Matlab code [#f1]_, which can also run in Octave


An introduction to EMST can be provided through an example simulation setup, in
this case for a slotted radial flux permanent magnet machine. This is a very common
topology of machine, and one of the most used in EMST. Later sections in this
documentation will delve into the other machine types available and the details of
the tools implementations.

The Radial Flux Permanent Magnet Machine
========================================

A machine design is always defined using fields of a structure which we usually
named 'design'. Throughout the toolbox this is the name used. In this example
we will set up a simulation of a radial flux permanent magnet machine, and take
a look at what is actually produced by the toolbox.

The machine to be simulated is a 3-phase, twelve-pole machine with an overlapping
winding. The machine design is arbitrary and probably not a particularly good
example of a machine design, and is purely for demonstrating the tools.

Winding Design
--------------

In this case, the winding arrangement is set up first. The way this is specified
is based on the Stator slot/coil/winding terminology presented in [#f2]_.

To summarise this the following variables describe a winding, in general, it is
not necessary to specify all of these, some can be calculated from a minimal set,
and different subsets may be used depending on what is easiest for the given
task.

============  =================================================================================================
Field         Description
============  =================================================================================================
Poles         total number of magnetic Poles in the machine (the number of magnets in permanent magnet machine)
Phases        the number of electrical Phases in a machine
Qs            total number of stator slots in all Phases combined
Qc            total number of winding coils in all Phases combined
yp            Average coil pitch as defined by (Qs/Poles)
yd            Actual coil pitch as defined by round(yp) +/- k q - number of slots per pole and phase
qc            number of coils per pole and phase, i.e. the ratio coils / (Poles * Phases)
qcn           numerator of qc
qcd           denominator of qc
============  =================================================================================================

The code for seting up the winding design for our example machine is then shown below:

**code here**

Coil Design
-----------

coils can either be single layered (only one coil side per slot, not overlapping in the slots), or double layered (two coils sit on top of each other in a slot). In each case the following relationships hold:

Single layer q = 2qc Qs = 2Qc Double layer q = qc Qs = Qc

The coil design is then determined using the following variables:

CoilFillFactor
CoilTurns

.. code-block:: matlab
   :linenos:

   % we will design a 4-pole machine, so it will have 4 magnets, and 2
   % pole-pairs
   design.Poles = 4;
   % Choose the number of Phases, the conventional 3
   design.Phases = 3;
   % The desired number of layers is stored in 'CoilLayers'
   design.CoilLayers = 2;
   % The type of winding is specified as a string, it can be 'overlapping' or
   % 'nonoverlapping'
   design.WindingType = 'overlapping';
   % Next we can specify the winding by stating the ratio of coils to Poles
   % and slots. This value must be a fraction object, which is created as
   % below. To simplify things we will use one coil per pole and phase, but
   % many other ratios are possible, see the help for "fr" (run "help fr" at
   % the command prompt) for more information on using fractions objects
   design.qc = fr(36,design.Poles*design.Phases);
   % Specify the actual coil slot pitch in slots
   design.yd = 4;
   % We must also specify a fill factor for the coils, this is the wire fill
   % factor achieved
   design.CoilFillFactor = 0.6;
   % We must also specify either the number of turns, of the diameter of the
   % wire used in the coils. If we specify the turns, the wire diameter is
   % calculated, and vice-versa. The wire diameter is specified in the field
   % "Dc" and the number of turns in the field "CoilTurns". In this case we
   % will specify the number of turns
   design.CoilTurns = 200;
   % The number of series coils, or parallel branches of coils in a phase is
   % controlled with the fields 'CoilsPerBranch' or 'Branches'. We must set
   % one or both of these fields. Here we will use all coils in series by
   % setting the number of branches to one.
   design.Branches = 1;


In addition, coils can either be single layered (only one coil side per slot, not overlapping in the slots), or double layered (two coils sit on top of each other in a slot). In each case the following relationships hold:

Single layer q = 2qc Qs = 2Qc Double layer q = qc Qs = Qc

With this new knowledge we can set up the winding design we will use for our example machine.








.. [#f1] Some parts of EMST are written in C++ or C. The non-Matlab components made accessible through mex files and most users will not need to be aware of the details of the implementation of these parts.

.. [#f2] J. J. Germishuizen and M. J. Kamper, "Classification of symmetrical non-overlapping three-phase windings," in The XIX International Conference on Electrical Machines - ICEM 2010, 2010, pp. 1-6.
