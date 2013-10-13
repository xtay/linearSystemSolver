#BiCGSTAB - a sparse linear system solver
============================
###basic idea###
**BiCGSTAB** is short for biconjugate gradient stabilized method, more details of this algorithm, see <http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method>

**CSR**, short for Compressed sparse row, is an effective way to store the sparse matrix. However, the solver currently supports only symmetric matrix storage, this solver stores each row from the diagnal element to the last nonzero element. For details, see <http://en.wikipedia.org/wiki/Sparse_matrix>

###compile and run###
change directory into this dir, to complie the project, execute the following shell command
```
    $ make
```

to run this program, first of all, you need to have a compressed matrix in file in "./data
", and a file contains the vector in directory "./data", then excute
```
    $./main ./data/$(FILE_OF_MATRIX) ./data/$(FILE_OF_VECTOR)
```

or, modify the shell script "run.sh", and then excute
```
    $. run.sh
```

###files and directory###
* **data**
    * a file named "mat" store the compressed matrix, which will be used as an input file of this solver, and a file name "vec", same as "mat".
* **util**
    * codes that do some basic operation around the Compressed matrix, include read matrix and vectors from a file, manage memory, basic matrix/vectors calculations, etc.
* **bicgstab**
    * codes that implements the BiCGSTAB algorithm.
* main.c
    * calls some functions provide by files above, and test them.
* run.sh
    * a shell that allows the user type less on the keyboard while running the program.

###Coming soon!!!###
* a general **CSR** way to store the matrix will be allowed.

* mpi involved, to make this solver a parallel one 

###Contact###
feel free to send me an email, if you have anything want to say about these codes
