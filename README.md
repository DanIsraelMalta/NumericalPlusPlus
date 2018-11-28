# NumericalPlusPlus

A header only collection of usefull, accurate & efficient numercial objects.

This repository holds the following numerical objects:
* FloatingPointTraits: a set of dedciated operations for floating point variables.
* Constants: varios constants.
* static_for: a simple compile time expanded for loop.
* Common: a set of general numerical utilities
* VectorN: fixed size numerical vector (designed for small vectors), include lots of numerical utilities and specialized operations.
* MatrixNM: fixed size rectangular (ROXxCOL) matrix (designed for small matrix) with a defined layout (row/colum major), include lots of usefull decompositions, equation solvers and a broad filled of constructors and modifiers.
* Quaternion: a quaternion object.
* Interval: a numerical interval object (not a 'standard' interval, but the correct numerical one; see  https://en.wikipedia.org/wiki/Interval_arithmetic).
* 3D: various 2D/3D related operations and decompositions which are not part of the VectorN/MatrixNM objects due to their specialized field of operation.

example usage and tests of the various modules is found in test.cpp.
