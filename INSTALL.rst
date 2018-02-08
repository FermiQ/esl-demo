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

MPI Support
-----------
```
cmake ../ -DWITH_MPI=On
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
