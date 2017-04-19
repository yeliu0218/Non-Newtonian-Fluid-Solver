Pelican: a Complex Fluid Solver
=====

This CFD solver is used for multi-phase flow simulations on arbitrary domain with rectangular meshes, using Volume-Of-Fluid(VOF) method.

The type fluid can be Newtonian/non-Newtonian.

The Navier-Stokes solver implements regularisation scheme and Augmented Lagrangian scheme. 

The interface is tracked by a Piecewise Linear Interface Construction(PLIC) scheme.

The code supports 2D and axisymmetrical geometry, and up to 3 phases of fluid together.

We are currently working on the implementation of surface tension term.

Installation of pelicans-3.1.0 (Linux)
-----

Pelicans is a library of tools used to develop computer software to solve partial differential equations
by focusing on the design of a numerical scheme suited to the physical model being studied. It is developed by 
<a href="http://www.irsn.fr/EN/Pages/home.aspx">IRSN</a>.

More details can be found <a href="http://www.irsn.fr/EN/Research/Scientific-tools/Computer-codes/PELICANS/Pages/PELICANS-software-platform.aspx">here</a>.

Prerequisites: (For Ubuntu 16.04)

gfortran, g++, freeglut3-dev, libblas-dev, liblapack-dev, libumfpack 5.7.1, libsuitesparse-dev, libspooles-dev, libmumps-dev(for parallel),
mpich(for parallel), graphviz, dot2tex, doxygen

All of the packages above can be installed by
 
`$ sudo apt-get install *`

There might be different versions of packages, for example, libumpack 5.7.1  is replaced by 5.6.0 for Ubuntu 14.04, just install the newest
version of each package if possible.

After all packages are installed, copy the two folders pelican-3.1.0 and Pel3.1_FluidSolver into one main directory. 
In my computer, it is

`$ ~/pelicans/pelicans-3.1.0`

and

`$ ~/pelicans/Pel3.1_FluidSolver`

We need to change something in the file, according to your computer settings

`$ cd pelicans-3.1.0`

`$ gedit Makefile`

Modify `LICENSE = "accept"`

If you want to use `icc` as the compiler instead of `gcc`, modify `CCC := gcc` as well.
Then

`$ cd etc` 

`$ gedit extra-Linux.mak`

Change 
`with_OPENGL = 1`
`with_UMFPACK = 1`
and the rest remains 0

And you should give it the correct path of OPENGL and UMFPACK as well, which you can find by command locate. For example

```Bash
 ifeq ($(WITH_UMFPACK),1)

     ifeq ($(MAKE_PEL),1)
     SRC += $(wildcard $(PELICANSHOME)/ExternalAPI/UMFPACK/src/*.cc)
     CPPFLAGS += -I$(PELICANSHOME)/ExternalAPI/UMFPACK/include
     endif

     UMFPACKPATH = /usr/lib/x86_64-linux-gnu
     CPPFLAGS += -I/usr/include/suitesparse
     LIBPATH  += $(UMFPACKPATH)
     LDLIBS   += -lumfpack -lamd
 
     WITH_BLAS = 1
     endif
```
And it is the same for all other packages if you change it from 0 to 1.

Now the correction of the files have been done, we proceed to compile pelican-3.1.0

`$ cd pelican-3.1.0`

`$ make clean`

`$ make all`

It takes a while to compile.

Then we need to modify bashrc 
`$ gedit ~/.bashrc`

At the bottom add the following lines:

```Bash
export PELICANSHOME=/home/leo/pelicans/pelicans-3.1.0
PATH=$PELICANSHOME/bin:$PATH
PATH=$HOME/bin:$PATH
```
Then open another terminal and type 
`$ source .bashrc`

It will read the modified bashrc and change it. 

Now if you type 

```Bash
echo $PELICANSHOME
```
it will give you the path of pelicans-3.1.0. In my case, it is

`$ /home/leo/pelicans/pelicans-3.1.0`

Now if you type `pel` it knows the command and it will give you all of the command possibilities with it.

Still in pelican-3.1.0 directory, type 
`$ make check`
to generate the exe file and to run a number of tests


Installation of Pel3.1_FluidSolver
------

`Pel3.1_FluidSolver` is based on the pelicans platform. Our main algorithm about solving non-Newtonian fluid is programmed in this package.
The installation is pretty simple, compared to `pelicans-3.1.0`

`$ cd Pel3.1_FluidSolver`

`$ make clean`

`$ make exe`

It will generate the execution file 

How to run a simulation
-------


In `Pel3.1_FluidSolver/RegressionTests` and `Pel3.1_FluidSolver/ExampleCases`, there are a number of example tests. 

Take `~/pelicans/Pel3.1_FluidSolver/RegressionTests/NavierStokes/AugmentedLagrangian/ElectricField` for an example,

data.pel is the script where you define the problem and solving scheme.

Go to the directory and type the command below to run the simulation

`$ pel run ~/pelicans/Pel3.1_FluidSolver/lib/Linux-gcc/opt0/exe data.pel result`

where result is the log file. 

Below is a simulation result of the spreading of viscoplastic fluid under gravity.
![](http://www.math.ubc.ca/~yliu0218/images/profi.jpg) 

More simulations can be found in my <a href="http://www.math.ubc.ca/~yliu0218/">website</a>







