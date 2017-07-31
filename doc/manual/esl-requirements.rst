Requirements, Coding Standards, and Best Practices
==================================================

Introduction
------------

The ESL aims at providing general, standardised and efficient modules and
libraries for the development of electronic structure codes. To do so, the
modules and libraries included in the ESL must follow a minimum set of
requirements, coding standards, and best practices to ensure that they are
maintainable, useful, and easy to use. Therefore, there was a great deal of
discussions among the community about this issue and a consensus was reached
about such a minimum set that all ESL modules and libraries should follow.
There was however the concern that a very strict set of rules would alienate or
hinder the contributions from volunteers. It was thus decided that this set
should act only as guidelines for volunteer contributions. The situation is,
however, quite different for new modules and libraries to be developed directly
by personnel hired by CECAM or other institutions, as well as for EU-funded
projects like E-CAM, and we believe this set should be mandatory for such
modules and libraries.


Scope
-----

This document defines the requirements, coding standards, and best practices of
the new modules developed by the computer scientists hired by CECAM in the
framework of the E-CAM project. They are mandatory to them, while they are only
guidelines for the modules contributed to the ESL by volunteer developers. The
main objective of this document is to guarantee that all newly developed code
will be both attractive to the industry AND usable by researchers.


Terminology
-----------

Within this document, we use the verbs “must” and “should” in the following
sense:

- must: if not fulfilled, the corresponding specification will compromise the
  future of the ESL;
- should: the corresponding specification is necessary to maintain
  compatibility between E-CAM and volunteer contributions, and can be
  discussed.


Programming languages
---------------------

It is critical that ESL modules be easy to incorporate in codes written within
the electronic structure community. The most widespread language in this
community is Fortran, followed by C and C++. Therefore, all libraries developed
within E-CAM must provide a full C implementation of the API, along with
Fortran 2003 interfaces. Providing Python bindings is highly recommended but
should be implemented only once the C and Fortran implementations are complete.
Exception to this are allowed if the library is meant to be only used from a
specific language. The full language recommendations for ESL contributors can
be found here: esl.cecam.org/mediawiki/index.php/ESL_language_recommendations


Documentation
-------------

Source code documentation must be managed with Doxygen, in order to support
integration into the ESL automatic documentation system and provide a unified
style with respect to volunteer contributions. Other kinds of documentation are
free-style but must provide a wiki version on esl.cecam.org and include at
least one tutorial for beginners.


Licences
--------

All ESL libraries and modules must be released under an open source license
that allows them to be usable by non open source software. The two most popular
choices for this are BSD and LGPL. There are advantages and disadvantages to
both, which can be quite subtle yet important; in general, BSD is more
permissive, while LGPL ensures that contributions to the library remain open.
We encourage software authors to research and consider which option is most
suitable for them.

The full licensing recommendations for ESL contributors can be found here:
http://esl.cecam.org/mediawiki/index.php/ESL_licensing_recommendations


Source-code management
----------------------

The development of the ESL modules and libraries should be done using Git and
the repository should be set either on gitlab.com or github.com under the
Electronic Structure Library organization, or on a server provided by CECAM
(e.g. https://gitlab.e-cam2020.eu/). Access control and permissions are managed
by the ESL Curators and the E-CAM Infrastructure Managers.


APIs
----

The API of each ESL library must include routines providing information about
both the current version of the library and the version of the API.


Error handling
--------------

Libraries must provide informative error handling and always return control to
the parent code. The library should only be allowed to stop the execution in
case there is a clear programming error (e.g., the parent code passed a null
pointer to a function that expects a non-null pointer).


I/O
---

Libraries should, as much as possible, be silent, i.e., not perform I/O
operations at all, except when I/O is precisely the purpose of the library
(e.g. Libescdf, Pspio). All data should therefore be passed in from and out to
the parent code through the API. However, for the larger libraries performing
complex and multi-step tasks this might prove to be extremely impractical, and
a better solution should be sought (and, if necessary, developed within E-CAM).
This must respect the goals set out in the Programming Languages section for
full Fortran and C interoperability.


Parallelism
-----------

ESL developers should be aware that the libraries are very likely to be called
from massively parallel codes that often rely on hybrid MPI/OpenMP parallelism.
This fact must be taken into account in the library design and all the ESL
libraries must be thread-safe.


Build systems
-------------

In order for the ESL libraries and modules to be easy to compile, portable, and
well integrated with each other, a standard build system type, based on e.g.
the Autotools or CMake, must be used. We recommend the Autotools because the
volunteer developers are more familiar with this framework, but any build
system supporting Fortran is suitable.
