# Settings which are necessary for UK R-matrix codes input files
#          except geometries (set in %geometry in a different file)

%model = (

  'directory',     "", # directory (relative path) which should correspond to a given model
                       # (you can use ${bs} as e.g. in "H2${bs}cc-pVTZ${bs}CAS")
                       # and which will be created in $dirs{'output'}
                       # if it is empty, it will be set later according to the model settings
                       # as H2/cc-pVTZ.CAS.0frozen.12active.10virtual.4states.r13.D2h
                       # suffixes for geometry and run will be also added if specified elsewhere

  # Molecule

  'groundStateSpinMultiplicity', 2,
  'molecule',    "H2OX",              # Used only for directory and descriptions
  'atoms',        ["O", "H", "H"],   # All atoms must be specified
  'nelectrons',     9,              # number of target electrons
  'symmetry',     "C2v",              # Options are "D2h", "C2v", "C2", "Cs", "C2h", "D2", "Ci", "C1"

  # Units

  'r_unit',        0,                 # distance unit: 0 - atomic units, 1 - Angstroms
                                      # for consistency reasons this will be actually set to the value given in $geometry{'length_unit'}
  'e_unit',        2,                 # energy unit:   0 - atomic units, 1 - Rydbergs, 2 - eV
  'x_unit',        1,                 # cross section: 1 - atomic units, 2 - Angstroms^2

  # Model - each model will have its own directory from these settings in the directory of the molecule
  # Important note: Put correct numbers of orbitals even for SE and SEP
  #
  # Choices below will determine which orbitals are generated and how many are used
  'basis',         "cc-pVQZ",         # make sure there are corresponding files in $dir_basis (see below)
  'orbitals',      "HF",              # which orbitals to use, HF = Hartree-Fock orbitals (default), natural = natural orbitals obtained by Molpro with CASSCF
  'charge_of',      "target",         # PHOTO/RMT: use orbitals for N-electron system (= "target") or (N+1)-electron system (= "scattering")
  'select_orb_by',  "energy",         # molden = use orbitals as ordered in Molden file, energy = use orbitals ordered by energy
  'model',           "CAS",          # options are "SE"    = static-exchange,
                                      #             "SEP"   = SE + polarization,
                                      #             "CAS"   = complete active space (as a default it is used the model B below),
                                      #             "CAS-A" = contracted version of CAS-B,
                                      #             "CAS-B" = standard close coupling with (core+cas)^N+1 and (core+cas)^N x (virtual)^1
                                      #             "CAS-C" = adds (core+cas)^N-1 x (virtual)^2 to CAS-B (to add more polarization)
                                      #   if you want to run FCI calculation then set "CAS" and 'nactive' = -1 and 'nvirtual' = 0
  'nfrozen',       1,                 # number of frozen target orbitals (used for SEP and CAS)
  'nactive',       6,                 # number of active target orbitals (used for SEP and CAS)
  # 'nactive',       8,                 # number of active target orbitals (used for SEP and CAS)
                                      #   if -1 then all available orbitals provided by the basis are used
  'nvirtual',      6,                 # number of virtual orbitals
  'nreference',   10,                 # number of orbitals used for searching reference orbitals

  # if the following arrays are empty (zeroes) then orbitals are chosen automatically according to their ordering,
  # ("energy" or "molden" for molpro orbitals, but currently just "energy" for psi4 orbitals),
  # otherwise the script uses chosen orbitals but ONLY IF the total number of chosen orbitals in arrays is consistent
  # with numbers of orbitals specified above in 'nfrozen', 'nactive', 'nvirtual'
  # be careful with these settings, if you choose your active space differently than orbitals are ordered
  # you should specify also which virtual orbitals you want to use
  'frozen_orbs',   [1,0,0,0], # C2v, which orbitals for each symmetry to use as frozen,
  'active_orbs',   [3,1,2,0], # C2v, which orbitals for each symmetry to use as active,
  # 'frozen_orbs',   [1,0,0,0], # Cs, which orbitals for each symmetry to use as frozen,
  # 'active_orbs',   [5,1,0,0], # Cs, which orbitals for each symmetry to use as active,
  # 'active_orbs',   [2,1,1,1], # which orbitals for each symmetry to use as active,
  'virtual_orbs',  [0,0,0,0], # which orbitals for each symmetry to use as virtual,
  'reference_orbs',[0,0,0,0], # which orbitals for each symmetry to use for searching reference orbitals,

  # ------------- The Multiple Active Spaces (MAS) approach -------------------
  #
  # Beyond the predefined models...
  #
  # It is strongly recommended to read the tutorial on the MAS approach, as not
  # all usage aspects are covered below.
  #
  # If 'use_MASSCF' is non-zero then the quantum chemistry calculation is performed using the MASSCF approach.
  # The options setting the orbital spaces above are overridden with the exception of 'nvirtual' and
  # 'virtual_orbs', MAS must then be used for the target and scattering calculation.
  # If 'use_MAS' is non-zero then the target and scattering calculation is performed using the MAS approach.
  # This option overides the setting of orbital spaces above, except for the total number of orbitals per
  # symmetry read in by scatci_integrals. So the orbitals used in the MAS approach must be less than or equal to
  # 'nfrozen' +  'nactive' + 'nvirtual', but there is no requirement to have the same number of frozen,
  # active or virtual orbitals as specified above. Any orbitals not included in the MAS description will be
  # considered virtual and contracted along with the continuum.
  'use_MASSCF', 0,  # == 0, don't use MASSCF approach. == 1, use ORMASSCF in Molpro.  == 2, use GASSCF in Molcas.
  'use_MAS',    0,  # == 0, don't use MAS approach. == 1, use ORMAS.  == 2, use GAS.

  # The format for defining the each active space is:
  # [orbs per sym] or orbs_per_subspace, [min_occ, max_occ], /(frozen|closed|active)/
  # where min_occ and max_occ are the local minimum and maximum electron occupancy for
  # the ORMAS case and the cumulative minimum and maximum occupancy in the GAS case.
  # Note: 'frozen_orbs' above actually refers to closed/inactive orbitals in the quantum
  # chemistry CASSCF context, not frozen orbitals. In the MAS approach the quantum chemistry
  # nomenclature is used.
  'qchem_MAS', undef, # Allows for different MAS for the quantum chemistry calculation.
  'qchem_constraints', undef, # Molpro allows simple constraints for ORMASSCF (see tutorial).

  'MAS', [
    1, [2,2], 'closed', # These orbitals are closed.
    6, [8,8], 'active', # Here we have 8 electrons in 6 orbitals CAS(8,6)
  # If a different ordering is desired than energy or molden, one can set the number of
  # orbitals per symmetry in each subspace as below.
  #  [1,0,0,0,0,0,0,0], [2,2], 'closed', # These orbitals are closed.
  #  [5,1,0,0,0,0,0,0], [8,8], 'active', # Here we have 8 electrons in 6 orbitals CAS(8,6).
  ],
  'constraints', undef, # e.g. sub {my $dist = shift; return !($dist->[1] > 1);}, # Additional constraints

  # If you want to define the l2 configuration space as something different to (N+1) electrons in
  # the target space (e.g. as in SEP, CAS-B, CAS-C etc.) then it can be done here. If constraints are set above then
  # the l2 space must be defined by hand, with care taken that the l2_constraints are consistent with
  # the constraints on the target space.
  # Below is equivalent to the CAS-B model:
  'l2_MAS', [
    1, [2,2], 'closed', # These orbitals are closed.
    6, [8,9], 'active', # Single excitations out of the CAS
    5, [0,1], 'active', # and into the virtual orbitals.
  ],
  'l2_constraints', undef, # Additional constraints.
  # ---------------------------------------------------------------------------

  # Number of target states used to generate orbitals  in MOLPRO CASSCF calculation
  # for closed-shell target, specify the number of singlets, triplets, ...
  # for   open-shell target,                       doublets, quartets, ...
  # 'ncasscf_states', {'doublet', [1,1,1,0],
  'ncasscf_states', {'doublet', [1,0,0,0],
                     'quartet', [0,0,0,0]},

  # Target states: for closed-shell target, specify the number of singlets, triplets, ...
  #                for   open-shell target,                       doublets, quartets, ...
  # for each symmetry (irreducible representation)
  # number of all target states which will be calculated
  'ntarget_states', {'doublet', [0,1,0,0],  # number of target states to calculate in each irreducible representation (IR)
                     'quartet', [0,0,0,0]},
  'ntarget_states_used', "1",        # number of target states which will be actually used in scattering calculations
                                      # (chosen according to their energy from states above)

  # Deletion thresholds, used in  scatci_integrals for the continuum orthogonalization
  'delthres',      ["1.0D-07", "1.0D-07", "1.0D-07", "1.0D-07", "1.0D-07", "1.0D-07", "1.0D-07", "1.0D-07"],

  # The options for PCOs below only need to be modified if you want to include pseudocontinuum orbitals in your calculation
  # PCO basis (to generate pseudostates 'model' should be set to CAS)
  'use_PCO',               0,                 # do you want to use a Gaussian type PCO basis?
  'reduce_PCO_CAS',        1,                 # == 0, off, target CSFs containing PCOs only use CAS^N-1 PCO^1 ; ==1, on, target CSFs containing PCOs only use GS^N-1 PCO^1.
  'maxl_PCO',              2,                 # the highest partial wave used in the PCO Gaussian basis
  'PCO_alpha0',   [0.15,0.15,0.15],           # Alpha0 parameter for PCO generation, one per PCO partial wave.
  'PCO_beta',     [1.3 ,1.3 ,1.3 ],           # Beta parameter for PCO generation, one per PCO partial wave.
  'num_PCOs',     [1   ,1   ,1   ],           # Number of PCO functions generated per PCO partial wave.
  'PCO_gto_thrs', [-1.0,-1.0,-1.0],           # Threshold for how close the cont gtos are allowed to the PCOs in terms of exponent, if <0 defaults to Alpha0*(Beta-1).
  # deletion thresholds, for PCOs
  'PCO_delthres',      ["1.0D-06", "1.0D-06", "1.0D-06", "1.0D-06", "1.0D-06", "1.0D-06", "1.0D-06", "1.0D-06"],

  # Settings for the (N+1)-electron calculationi, i.e. the scattering states

  # Scattering spin-symmetry: for closed-shell target, specify which of doublets, quartets, ... to use
  #                           for   open-shell target,                  singlets, triplets, ...
  'scattering_states', {'singlet', [1,1,1,1], # 1/0 to run/not to run scattering calculation for a given spin-symmetry
                        'triplet', [1,1,1,1]},

  # Continuum basis set
  'use_GTO',             1,              # do you want to use a Gaussion type basis?
  'radius_GTO',         13,              # radius for which Gaussian basis was optimized (basis sets available for radius 10, 13, 15, 18)
  'maxl_GTO',            5,              # the highest partial wave used in the continuum Gaussian basis
                                         # search the directory basis.sets to find out
                                         # which continuum bases are available (they start with swmol3.continuum...)

  'use_BTO',             0,              # do you want to use a B-sline type basis?
  'start_BTO',         8.0,              # radius where the B-splines start
  'order_BTO',           9,              # order of the B-splines
  'no_of_BTO',          12,              # number of B-splines
  'maxl_BTO',            4,              # the highest partial wave used in the continuum B-spline basis
  'maxl_legendre_1el',  70,              # maximum L to use in Legendre expansion for nuclear attraction integrals in the B-sline basis
  'maxl_legendre_2el',  55,              # maximum L to use in Legendre expansion for 2-electron integrals in the B-sline basis

  'rmatrix_radius',   13.0,              # R-matrix radius used for scattering, it is also where the BTO basis ends

 # Propagation step
  'max_multipole',       2,              # maximum multipole to be retained in expansion of long range potentials
  'raf',              70.0,              # radius at which continued fraction method can be used for R-matrix propagation

 # Energy grid
  'nescat',   "40          , 40          ,  40          , 500",              # number of input scattering energies in each subrange (input for R_SOLVE via namelist &rslvin)
  'einc',     "1e-4, 2.5e-5, 1e-3, 2.5e-4,  1e-2, 2.5e-3, 0.1, 0.005",              # scattering energies - initial energy, energy increment, there can be more subranges

 # Initial and final states for which the cross sections are calculated
  'maxi',              "1",              # the highest initial state for which cross sections are required
  'maxf',              "0",              # the highest final state for which cross sections are required (zero means all)

  # SCATTERING only setting
  'positron_flag',       0,             #0 for electron scattering, 1 for positron scattering

  # PHOTOIONIZATION only settings
  'initialsym',        "1",              # Symmetry of the (N+1)-electron initial state
  'first_Ip',          0.0,              # first ionization potential (used by dipelm only) - units set by e_unit
                                         # If Ip is set to 0. then dipelm will calculate the Ip automatically.

  # Options for printing information to logs
  'norbitals_to_print', 10,              # maximum number of orbitals to print in logs for each spin symmetry
  'nstates_to_print',   10,              # maximum number of states to print in logs for each spin symmetry


  # Other settings specific for R-matrix codes
  'lndo',          10000000,          # memory control in CONGEN
);
