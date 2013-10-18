#BiCGSTAB - a sparse linear system solver
============================
###Basic idea###
**BiCGSTAB** is short for biconjugate gradient stabilized method, more details of this algorithm, see <http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method>

**CSR**, short for Compressed sparse row, is an effective way to store the sparse matrix. However, the solver currently supports only symmetric matrix storage, this solver stores each row from the diagnal element to the last nonzero element. For details, see <http://en.wikipedia.org/wiki/Sparse_matrix>

**MPI**, short for Message Passing Interface, a language-independent communication protocol used to program parallel computers. This program is developed and debuged on "mpich2-1.1.1p1", for details about mpi, see <http://www.mpi-forum.org> and <http://www.mpich.org>

###Compile and run###
Change directory into this dir, to complie the project, execute the following shell command
```
    $ make
```

To run this program, first of all, you need to have a compressed matrix in file in "./data
", and a file contains the vector in directory "./data", then excute
```
    $./main ./data/<yourMatrixFile> ./data/<yourVectorFile>
```

Or, modify the shell script "run.sh", and then excute
```
    $. run.sh
```

###Files and directory###
* **data**
    * A file named "mat" store the compressed matrix, which will be used as an input file of this solver, and a file name "vec", same as "mat".
* **util**
    * Codes that do some basic operation around the Compressed matrix, include read matrix and vectors from a file, manage memory, basic matrix/vectors calculations, etc.
* **bicgstab**
    * Codes that implements the BiCGSTAB algorithm.
* main.c
    * Calls some functions provide by files above, and test them.
* run.sh
    * A shell that allows the user type less on the keyboard while running the program.

###Coming soon!!!###
* A general **CSR** way to store the matrix will be allowed.

* Uhh...I'm thinking of add some CUDA or OpenCL and OpenMP features into it.
    * A giant machine can be devided into nodes, and nodes then devided into cores, maybe a Graphics card is attached to each node. So, MPI is useful to manage processes on different nodes, while OpenMP is perfect choice to do the job within a node which has many cores, and CUDA or OpenCL is used to accelorate the program using the Graphics card

###Contact###
Feel free to send me an email, if you have anything want to say about these codes
