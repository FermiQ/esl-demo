Overview of The Electronic Structure Library
============================================

Summary and motivations
-----------------------

The Electronic Structure Library (ESL) is a community-maintained collection of
software of use for electronic structure simulations supported by CECAM (Centre
Européen de Calcul Atomique et Moléculaire). It is an extended library – in a
generic sense – that can be employed by everybody to build their own scientific
packages and projects. It features a structured wiki which contains entries
documenting functionalities, algorithms, interfaces, standards and pieces of
code, ranging from small routines for performing simple tasks, all the way up
to complete software libraries. Although wiki entries may document software
maintained elsewhere, the ESL also provides a software development
infrastructure for internally developed projects. This infrastructure also
permits the regular building of automatic documentation, for both internal and
external projects.

The ambition of the ESL is to segregate layers of functionality within
scientific software modules which are general, standardised and efficient. In
this way, new ideas and new science can be programmed by scientists without
needing to rewrite functionalities that are basic and already well-established,
and without having to know more software engineering than science. In other
words, separate the coding efforts for cutting-edge research from the
underlying software infrastructure which requires maintenance and refactoring
at every step of the hardware race. This objective is achieved by providing the
scientific community with organised and comprehensive information about the
existing concepts and pieces of code, platforms and methodologies to develop
and share reusable software, as well as regular training for both the
scientists and the developers of the low-level software layers.

The ESL started in July 2014 through an extended 6-week-long workshop organised
by CECAM, with 3 days dedicated to preliminary discussions and the rest of the
time to coordinated software development and documentation efforts, and to
which 20 persons participated in total. A “CECAM Summer of Code”, so to speak.
The purpose of this event was to create the core of the ESL and provide
sufficiently useful information to scientific software developers, so that they
would start to contribute themselves in the future. The long-term objective of
this intiative is to create a privileged and friendly environment for the
sharing and pooling of software, in order to put a definitive end to a wasteful
duplication of efforts and mistakes that has occured in the past in the
electronic structure community. Actually, the ESL could start because of
successful previous efforts and thanks to a critical mass of scientists in
favor of developing a common software basis for future projects.


Contents of the ESL
-------------------

The central component of the ESL is its wiki, available at
http://esl.cecam.org/. The contents of the library are organised in the form of
ESL entries, with one wiki page per entry. ESL entries are categorised as
Functionalities, Algorithms, Generic interfaces, Application Development
Interfaces (APIs), Data standards, and Software, listed as category tags at the
bottom of each wiki page. One entry can contain more than one item of each
category. However, the creation of single-category entries linked to one
another is encouraged, in order to allow an easily-navigable network to form
and expand. In particular, entries describing functionalities should remain
separated as much as possible, since each functionality will typically link to
several different algorithms and implementations.

The ESL grows with contributions open to a wide community, without imposing the
rigidity implied by the a priori definition of a conventional library, while
guaranteeing a sufficient level of practical organisation. It is also intended
to be more than just a software repository. The expectation is for the ESL to
be useful and used over a long period of time, and through that to be able to
progressively standardise software, APIs (Application Programming Interfaces)
and data standards of common use in the electronic structure community, as
opposed to randomly accumulating disconnected bits of code.


Functionalities
~~~~~~~~~~~~~~~

This entry describes a particular functionality typically taken care of by a
piece of software. This could be anything from a small routine up to a complete
library. It should represent a specific process or component of an electronic
structure calculation, but may be broad enough to allow different
approaches/methods/algorithms or a variety of smaller tasks within it. Examples
are:

- Molecular dynamics (many possible algorithms for different ensembles, thermostats, ...);
- Linear algebra (a collection of mathematical operations within a single framework);
- Memory management (a general software requirement allowing for various approaches and tasks).


Algorithms
~~~~~~~~~~

This entry describes an algorithm used to perform a functionality, or a
specific task within a functionality. This might be a physical theory or
mathematical algorithm (i.e., a particular eigensolver method), but can also
refer to a purely computational task.


Generic interfaces
~~~~~~~~~~~~~~~~~~

This entry describes an interface for a possible software implementation of an
algorithm in generic terms. It should be human-readable and descriptive,
listing the physical quantities involved (e.g., the Hamiltonian, the lattice
vectors, ...) rather than specifying how these are to be represented within the
code.


APIs
~~~~

This entry describes an API in strict terms, i.e. gives all the information
necessary for an implementation with no ambiguities.


Data standards
~~~~~~~~~~~~~~

This entry contains the specifications for storing particular types of data in
files. It provides all the necessary details for any program to be able to read
and write according to the standard.


Software
~~~~~~~~

This entry describes and links to a piece of code implementing the above
categories. The code can be specific to a single API/algorithm/functionality or
bundle together more than one of each. It includes a download link, details
about licenses, code authors, and the responsible person for the wiki entry.
The latter will normally be the creator of the entry, but can be transferred if
necessary. Appropriate documentation can be given either in the wiki entry
itself, or by linking to an external page.


Software project management
---------------------------

The ESL provides a framework for project management hosted in Launchpad, a
large open web-based platform used by the Ubuntu Linux distribution and
managing more than 35,000 projects, 800,000 development branches, as well as
1.2 millions of bug reports. The entry point of the framework consists in a
project group, located at https://launchpad.net/esl, which software projects
can join to show that they implement functionalities defined in the ESL wiki.
Launchpad provides code and package repositories, a bug tracking system, Q&A
forums, translation management, mailing lists, user and group management,
release planning and tracking, and secure authentication, covering all the
needs of most projects regarding the organisation of collaborative development.
It has been chosen over GitHub for 3 reasons:

- Bazaar (http://bazaar.canonical.com/), the Version Control System (VCS)
  Launchpad is based on, is more portable, with less hardware and software
  requirements, and easier to learn than Git while equally available, which is
  important when contributors are scientists with a limited culture in
  computer-related aspects and working on non-standard supercomputers;
- the projects related to the ESL are very modular and with a reduced number of
  contributors, making Git – which has been initially designed for massively
  distributed development – over-dimensioned for most of them;
- Launchpad focuses on project management exclusively and provides
  finer-grained control than GitHub, while CECAM provides more flexible website
  hosting than GitHub, which makes the latter redundant for ESL projects.

This organisational model divides the ESL projects in 3 categories:

- internal, for projects fully hosted on the ESL server and Launchpad;
- mirrored, for projects having their own repositories and websites but with
  synchronized copies hosted on the ESL server and/or Launchpad;
- external, for fully autonomous projects only referenced by hyperlinks in the
  ESL wiki.

Of course, projects can switch from one category to another at any time,
depending on the evolution of their needs and resources. Furthermore, the ESL
does not impose any particular development model and gives full autonomy to the
contributors of a project to organise their work as they prefer, while sharing
guidance and resources with those in need of them. Beyond purely technical
aspects, it is important that the projects are gathered within a bottom-up
structure, in order to guarantee their continued existence whatever happens to
the ESL.


Contributing to the ESL
-----------------------

Creating and maintaining ESL entries is open to all. Wiki users are free to
create and contribute to categories at all levels. Guidelines and templates are
provided to preserve a sufficient separation between categories in the wiki
over time. Contributions are edited a posteriori by a group of volunteers, in
order to fix language-related issues and unify the writing style among the
pages. The rationale behind this methodology is that the easier it is to browse
the wiki, the easier it will be for new contributors to link their projects to
existing pages and help the network grow as a coherent whole. When necessary,
discussions can also take place through a forum directly integrated with the
wiki, which facilitates the understanding between users and makes the wiki more
friendly towards newcomers.


Resources and sustainability
----------------------------

The infrastructure on which the ESL relies is provided by CECAM (wiki,
automatic documentation) and Canonical, Ltd. (project management), two healthy
European entities which have been highly involved in research and innovation
for many years and understand very well what is at stake with scientific
software development. They provide an extremely valuable support to ensure the
viability of the ESL, in particular when it comes to enhance its visibility.
The correct and efficient use of this infrastructure depends mainly on the
electronic structure community, which has acquired a lot of maturity regarding
the technical and social aspects of scientific software development over the
last 10 years, through a broader and broader decompartmentalization of its
research areas and projects. There has actually been a clearly expressed need
for such an initiative within this community since 2010, due in great part to
the announced end of Moore's law for microprocessors, while meeting current
challenges requires more and more of them. CECAM has made a bold and smart move
proposing a paradigm shift in scientific software development models at just
the right time.

On the longer run, the ESL will have a significant impact on the quality,
usability, standards-compliance, interoperability, availability, and
performance, of scientific software related to electronic structure, since each
of its base components will be developed for the whole community and will be
fairly up-to-date regarding the constantly-evolving computer architectures.
Larger developer and user communities than just one research group will allow
for faster and more complete debugging, as well as result in greater creativity
and diversity when it comes to propose solutions. It will promote collaboration
over competition within the whole electronic structure community, which in turn
will lead to a more efficient use of the available resources and quickly
increase the corpus of open-access data usable for innovation. Once
sufficiently visible, it will also become an asset in further closing the gap
between research and industry, thus facilitating highly-competitive innovation
in Europe, in particular for the design of new materials, for nanotechnology,
and for biochemistry.

Having taken useful lessons from the pilot event which lead to the creation of
the ESL, CECAM is now contacting other scientific communities, thus enlarging
the scope of its initiative in parallel to the development of this first
library. Following this approach, it will be possible to build progressively a
European network of scientific communities promoting higher-quality standards
for both software and collaborations in a coherent way, greatly facilitating
the sharing of information and good practices between different fields of
research all over the ERA. This multidisciplinary and multiscale environment
will bring a lot of new opportunities that could not even be imagined before
and will further open the doors to fruitful collaborations with
highly-innovative companies. A key ingredient for the success of this evolution
will be the presence, at strategical positions, of qualifed persons with hybrid
profiles and broad cultures covering engineering aspects as well as scientific
ones, in order to accelerate the channeling and refactoring of the scientific
ideas into innovative products and processes, thanks to the development of a
common vision and understanding.
