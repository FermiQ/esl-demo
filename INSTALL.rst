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
 
