# SpinWave
This package computes the linear spin wave theory dispersion of magnons in crystalline materials, with a particular emphasis on the correct handling of anisotropic dispersions.

A ground-up rewrite of the popular SpinW package for 2D (3D coming soon!). Cross-compatibility is neither present nor planned.
The idea is to manufacture a module with a more transparent structure, and designed purely for pythonic open-source scientific computing backends.

The basic workflow is the following - 

- Define where the atoms are in space
- Draw bonds between them, giving them names in the process
- Specify the basis with respect to which the pin-spin coupling is defined
- Assign a specific set of numbers to the bond strengths and diagonalise the bogoliubov Hamiltonian

Once defined, it's possible to
- Calculate the spin-wave mode structure
- Compute the linear spin correlator response, with and without magnetic form factors

## TODO

1. Clean up the implementation: Separate computational parts from GUI dependent parts. Meant to support the workflow
(Build model prototype) -> (Export model [define some sensible format for this, JSON?]) -> Load model in headless compute cluster -> Optimise parameters
2. Fix all 2D-dependent variables for full 3D calculation (straightforward but irritating)
3. Classical ground state optimisation routines -- Some more rigorous optimisers would do a better job
4. Implement 2D plot fitting
5. Magnon damping (speculative)
