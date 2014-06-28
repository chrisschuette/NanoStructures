NanoStructures
==============

NanoStructures is a Dynamical Mean Field Theory (DMFT) solver for layered, strongly correlated nanostructures. The Density Matrix Numerical Renormalization Group (NRG) is used to solve the effective impurity problems. Parallelization using MPI allows to distribute the computational load among compute nodes in HPC cluster environment.

Required libraries
------------------
The following libraries are needed to successfully compile the project*
* [libconfig](http://www.hyperrealm.com/libconfig/)
* [GNU Scientific Library](http://www.gnu.org/software/gsl/)
* [GNU MPFR Library](http://www.mpfr.org/)
* [GNU Multiple Precision Arithmetic Library (GMP)](https://gmplib.org/)
* [FFTW 3](http://www.fftw.org/)
* [Boost.Program_options](http://www.boost.org/doc/libs/1_55_0/doc/html/program_options.html)

Documentation
-------------
The [documentation](http://chrisschuette.github.io/NanoStructures/documentation) can be found [here](http://chrisschuette.github.io/NanoStructures/documentation). Although the relevant section of the source code (essentially the NRG and DMFT class) are quite well documented, it is very hard for a reader not familiar with the NRG or DMFT to understand the source code. Consultation of the following references is therefore strongly advised

References
----------
The following references provide an introduction to the DMFT and the NRG

#### Dynamical Mean Field Theory ####
* [Georges, A., Kotliar, G., Krauth, W., & Rozenberg, M. J. (1996). Dynamical mean-field theory of strongly correlated fermion systems and the limit of infinite dimensions. Reviews of Modern Physics, 68(1), 13.](http://journals.aps.org/rmp/abstract/10.1103/RevModPhys.68.13)
* [Kotliar, Gabriel, and Dieter Vollhardt. "Strongly correlated materials: Insights from dynamical mean-field theory." Physics Today 57.3 (2004): 53-60.](http://www.physics.rutgers.edu/~kotliar/papers/PT-Kotliar_57_53.pdf)
* [Pavarini, Eva, et al. The LDA+ DMFT Approach to Strongly Correlated Materials–Lecture Notes of the Autumn School 2011 Hands-on LDA+ DMFT. 2011.](ttp://juwel.fz-juelich.de:8080/dspace/handle/2128/4467)

#### Numerical Renormalization Group ####
* [Bulla, Ralf, Theo A. Costi, and Thomas Pruschke. "Numerical renormalization group method for quantum impurity systems." Reviews of Modern Physics 80.2 (2008): 395.](http://journals.aps.org/rmp/abstract/10.1103/RevModPhys.80.395)
* [Krishna-Murthy, H. R., J. W. Wilkins, and K. G. Wilson. "Renormalization-group approach to the Anderson model of dilute magnetic alloys. I. Static properties for the symmetric case." Physical Review B 21.3 (1980): 1003.](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.21.1003)
* [Krishna-Murthy, H. R., J. W. Wilkins, and K. G. Wilson. "Renormalization-group approach to the Anderson model of dilute magnetic alloys. II. Static properties for the asymmetric case." Physical Review B 21.3 (1980): 1044.](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.21.1044)

