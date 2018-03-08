.. _required-knowledge:

******************
Required Knowledge
******************

This section details background knowledge of Matlab (or Octave)
which is required to use the software, and common coding idioms
which occur frequently and which novice Matlab users may not have
encountered.

.. _required-knowledge-pv-pairs:

Parameter-Value Pairs
=====================

The functions and classes in this software make extensive use of the
concept of Parameter-Value pairs for passing optional arguments.
This means constructs like the following will be seen frequently::

    x = examplefcn (y, 'SomeOption', some_variable)

In the above the example function has one required input (``y``) and
can then take at least one other optional argument. The optional argument
is supplied, by giving its name in one input argument, as a string,
followed by the value for that option in the following input argument.
In the example above, the ``SomeOption`` option is given the value ``some_variable``.

A slightly less contrived example may make this clearer. Imagine we
define a function ``checkIsScalar``. This function takes one required
input ``val`` and returns true if it is a scalar and false if
it is not. e.g.::

    >> val = 1;
    >> result = checkIsScalar (val)

    ans =

      logical

       1

    >> val = [1, 2, 3];
    >> result = checkIsScalar (val)

    ans =

      logical

       0

Now imagine checkIsScalar also takes an optional argument ``'PositiveOnly'``,
which is a logical value which changes the behaviour of checkIsScalar. If
``'PositiveOnly'`` is set to true, checkIsScalar will only return true
if the input is scalar and greater than zero. By default, if the ``'PositiveOnly'``
option is not supplied, it is false (so values can be positive or negative).
Then we can use this option as follows::

    >> val = -1;
    >> result = checkIsScalar (val)

    ans =

      logical

       1

    >> val = -1;
    >> result = checkIsScalar (val, 'PositiveOnly', true)

    ans =

      logical

       0

The default behaviour can also be had with the use of the option::

    >> val = -1;
    >> result = checkIsScalar (val, 'PositiveOnly', false)

    ans =

      logical

       1

Functions may take multiple optional arguments like so::

    x = examplefcn (y, 'SomeOption', some_variable, 'OtherOption', other_variable)

And the optional arguments may be supplied in any order::

    x = examplefcn (y, 'OtherOption', other_variable, 'SomeOption', some_variable)


The Matlab Path
===============

The Matlab path is a list of directories where Matlab knows to look
for code files, so they don't have to be in the current working
directory when you want to run them. See the official Matlab
documentation for this `here`__.

.. __: https://uk.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html

Namespaces (+package directories)
=================================

The software makes extensive use of namespaces. Namespace are created by
placing code files in a directory starting with a '+'. As an example of
what a namespace is, and how it works, imagine we create a folder named
``+mypackage``. We then put two function files in this folder named
``times2`` and ``add2``. When we want to call the functions in a Matlab
script, or on the command line, instead of just using their name, we
must use the following syntax::

    x = mypackage.times2(4)

and::

    x = mypackage.add2(10)

With the folder name, but without the `+` followed by a dot, followed
by the function name. Scripts, functions and classes can all be put in
namespaces in this way. In Matlab jargon, these special folders are
referred to as packages.

More detailed information on the concept can be found in The Mathworks
official documentation `here`__.

.. __: https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html


Classes (classdef)
==================

Another programming method used extensively is object-oriented programming
using the ``classdef`` syntax. See the official Matlab documentation for this
style of programming `here`__.

.. __: https://uk.mathworks.com/help/matlab/object-oriented-programming.html
