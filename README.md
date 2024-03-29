# $k$-evolution

$k$-evolution (arXiv:1910.01104) is basen on gevolution N-body code (Adamek J., Daverio D., Durrer R., and Kunz M., Nature Phys. 12, 346 (2016), arXiv: 1604.06065).

$k$-essene in the effective field theory of dark energy framework is added into the gevolution code.

The free parameters to be chosen are $w$ (the equation of state for dark energy) and $c_s^2$ (the speed of sound squared). 


## Compilation and usage

Before compilation, make sure that all required external libraries are
installed:

* LATfield2 [version 1.1](https://github.com/daverio/LATfield2.git)
* FFTW version 3
* GNU Scientific Library (GSL) including CBLAS
* HDF5

Make sure that the include paths are set properly, or add them to the
makefile. Also check the compiler settings in the makefile. The code is
compiled by typing:

    make

A typical command to run a simulation looks like this:

    mpirun -np 16 ./gevolution -n 4 -m 4 -s settings.ini

For further about gevolution information, please refer to the User Manual (manual.pdf)

## Credits

If you use $k$-evolution for scientific purposes, we kindly ask you to cite
*Farbod Hassani, Julian Adamek, Martin Kunz, Filippo Vernizzi, JCAP12(2019)011*
in your publications.

For bug reports and other important feedback you can contact the authors,
farbod.hassani@gmail.com (for queries related to $k$-evolution) and to Julian Adamek (gevolution related queries)
developers@latfield.org (for queries related to the LATfield2 library).

