# Silicon-Modelling
Code base for exchange CI code and electron shuttling code (TODO: and nextnano)

Currently there are two main folders in the repo:
- The ExchangeCI folder contains all pertinent files for a configuration interaction calculation to find the exchange energy of a double quantum dot system in any semiconductor (assuming the effective mass approximation).  A detailed overview of the math behind the calculation can be found at [OVERLEAF DOCUMENT LINK]
- The ElectronShuttling folder contains all pertinent files for an electron shuttling simulation.  It assumes you have produced or simulated the potentials for your particular gate geometry.  It uses a split operator method to evolve the wave function and can also perform Stark shift calculations assuming you have simulated 2D potentials (in order to find the electric field).
