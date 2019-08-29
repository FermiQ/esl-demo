Build with cmake
=================

requirements: 
 - libfdf >= 0.1.1 ; build in serial mode

eg of building
 ```
 mkdir build 
 cd build 
 cmake ../
 make -j
 ```

 running
 ```
 bin/esl-demo.X <input file>
 ```
 by default will look to read from sample.inp

ELSI Support
------------
ELSI support is required for PW calculations.
```
cmake ../ -DWITH_ELSI=On
```

MPI Support
-----------
```
cmake ../ -DWITH_MPI=On
```

if you get undefined references to MPI_Get_library_version you have a very old MPI implementation
Solution:
  - update the MPI
  - or
```
  FFLAGS="-DOLDMPI" cmake ../ -DWITH_MPI=on
```

Shared libraries
-----------------
```
cmake ../ -DBUILD_SHARED_LIBS=on
```

to build API documentations
---------------------------

```
cmake ../ -DWITH_DOC=on
```
