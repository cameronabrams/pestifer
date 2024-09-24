``domainswap``
==============

Parameters controlling a domain swap TMD simulation

Single-valued parameters:

  * ``swap_domain_def``: VMD atomselect string for the domain; will be appended with "and chain X" for each chain

  * ``anchor_domain_def``: VMD atomselect string for the anchor domains; will be appended with "and chain X" for each chain

  * ``chain_directional_swaps``: list of pairs of directional swaps of chainIDs [[A,B], [B,A]] swaps A and B

  * ``force_constant``: force constant used in biases for domain-swap CV MD simulation (200) (default: 200.0)

  * ``target_numsteps``: number of timesteps to run the biased MD simulation (default: 10000)

  * ``nsteps``: number of timesteps to run the biased MD simulation (default: 10000)

  * ``dcdfreq``: number of timesteps between dcd output (default: 100)

  * ``xstfreq``: number of time steps between cell size output to DCD file (default: 100)

  * ``temperature``: temperature of thermostat in MD run (default: 300)

  * ``ensemble``: ensemble in which MD is run [NVT] (default: NVT)



