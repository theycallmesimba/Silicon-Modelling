# Silicon-Modelling
Code base for exchange CI code and electron shuttling code in quantum dot geometries.

There are currently two main folders in the repo:
- The ExchangeCI folder contains all pertinent files for a LCHO-configuration interaction calculation to find the many electron energy spectra energy of an arbitrary quantum dot system in any semiconductor heterostructure (assuming an effective mass approximation). A detailed overview of the math behind the calculation can be found at https://www.overleaf.com/read/jgyhxnfhvnzr
- The ElectronShuttling folder contains all pertinent files for an electron shuttling simulation.  It assumes you have produced or simulated the potentials for your particular gate geometry.  It uses a split operator method to evolve the wave function and can also perform Stark shift calculations assuming you have simulated 2D potentials (in order to find the electric field).  Effective Hamitonian simulations including orbital, valley and spin physics can also be done as well.  This code was used to generate the results in the paper https://arxiv.org/abs/2003.08018
