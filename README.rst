ESL Demonstrator - A basic DFT Code
===================================

.. image:: https://gitlab.e-cam2020.eu:10443/esl/esl-demo/badges/master/build.svg
   :alt: Build status
   :target: https://gitlab.e-cam2020.eu:10443/esl/esl-demo/commits/master

``esl-demo`` demonstrates how to use Electronic Structure Library
components to build and run a basic Density Functional Theory code.

The demonstrator implements both a plane-wave code and a atomic orbital
code.

The ESL demonstrator documentation may be found `here <esl-demo-doc>`_.

.. _esl-demo-doc: http://esl.e-cam2020.io/esl-demo//public/


Minimum requirements
--------------------

The following software elements are required to play with this project:

- A C compiler (e.g. gcc, icc, pgcc, xlc, ...)
- A Fortran compiler (e.g. gfortran, ifort, pgif90, xlf, ...)
- A recent Python interpreter (Python 2 >= 2.7.13 or Python 3 >= 3.5.1)
- GNU Make >= 2.80 (other versions may work but have not been tested)

To build the documentation, you will need a recent version of `Sphinx`_.

.. _Sphinx: http://sphinx-doc.org/


Developer requirements
----------------------

If you clone the Git repository of ``esl-demo``, you will need the following
packages to build the demonstrator:

- Autoconf >= 2.69
- Automake >= 1.15
- Libtool >= 2.4.2
- M4 >= 1.4.17
- cmake > 3.0.2

ESL Components
--------------

This ESL Demonstrator makes use of the following ESL components:

- ESCDF
- FDF
- LibXC
- OMM Bundle
- Pspio
- SQARE
- ESLW_Drivers
- flook

They will be downloaded and built automatically the first time you try to build
``esl-demo``.
