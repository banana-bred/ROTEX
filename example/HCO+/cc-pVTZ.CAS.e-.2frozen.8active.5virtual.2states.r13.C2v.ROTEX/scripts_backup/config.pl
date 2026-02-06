# Running and saving options

%run = (

  'suffix',          "", # string added to model directory to distinguish different runs which have the same model
                         # it can be useful e.g. if you run only target calculations etc.

  'print_info',      "file", # "screen" or "file" or "both" or "none" - option whether info messages will be printed on the screen or into files in the subdirectory logs

  'add_files_to_backup', ["ukrmollib.pm", "MultiSpace.pm"], # all files given on the command line as "perl main.pl config.run.pl ..." are copied into directory script
                                           # if you want to backup additional files, specify them here

  # Running options
  # Which codes will you run? Choose one or other
  'molpro',          0,  # do SCF calculations using Molpro
  'psi4',            1,  # do SCF calculations using Psi4, a free alternative to Molpro
  # What type of calculation/process?
  'scattering',      1,  # run all programs, if you want to run target only, set to 0
  'photoionization', 0,  # calculate photoionization cross sections of the (N+1)-electron system instead of electron scattering cross sections
  'rmt_interface',   0,  # form the RMT molecular input file (rmt_interface and photoionization are mutually exclusive options)
  # Ocassionally, you may want to tun these parts of the calculation too
  'skip_radden',     1,  # skip calculating radial densities
  # Options to keep or delete files
  'gather_data',     1,  # gather eigenphase sums, cross sections, energies ...
  'clean',           1,  # removing fort.* etc (except moints with molecular integrals)
  'remove_moints',   1,  # remove moints with molecular integrals
  # Options to use existing input data
  'use_templates',   1,  # set this to 0, if you need to modify generated inputs manually and than rerun everything
                         # but do not change the filenames !

  # MPI launcher options
  'mpi_integrals',  "",  # MPI launcher command for scatci_integrals (e.g. "mpirun -ilp64 -n 1")
  'mpi_scatci',     "",  # MPI launcher command for mpi-scatci (e.g. "mpirun -np 32"; if empty, legacy scatci will be used)
  'mpi_rsolve',     "",  # MPI launcher command for mpi-rsolve (e.g. "mpirun -np 32"; if empty, serial rsolve will be used)

  # Run-time options for SCATCI_INTEGRALS
  'buffer_size',  5000,  # Size of temporary arrays for integral transformation in MiB
  'delta_r1',     0.25,  # Length, in Bohr, of the elementary radial quadrature needed for evaluation of the mixed BTO/GTO integrals
  'transform_alg',   0,  # Choice of the integral transformation algorithm: 0 = auto, 1 = sparse (optimal for BTO and mixed BTO/GTO continuum), other = dense
                         # The sparse transformation is not available in distributed (MPI) mode.

  # Saving options
  # These are files that can be plotted
  'save_eigenph',    1,  # set 1 to copy fort.10XX into a file like eigenph.singlet.Ag (eigenphase sums)
  'save_xsec',       1,  # set 1 to copy fort.12XX into a file like xsec.singlet.Ag (cross sections)
  # These are files that can be used to start an outer region calculation
  'save_channels',   1,  # set 1 to copy fort.10  into a file like channels.singlet.Ag
  'save_rmat_amp',   1,  # set 1 to copy fort.21  into a file like ramps.singlet.Ag
  # These are files that can be used for calculations beyond UKRmol+
  'save_Kmatrix',    1,  # set 1 to copy fort.9XX into a file like K-matrix.singlet.Ag
  'save_Tmatrix',    0,  # set 1 to copy fort.11XX into a file like T-matrix.singlet.Ag

  'keep_inputs',     1,  # set 1 to keep input files for UK R-matrix codes
  'keep_outputs',    1,  # set 1 to keep output files of UK R-matrix codes

  # Parallelization - you can specify a number of processes to run in parallel
  #                   separately for geometries and for symmetries (ireducible repr.)
  # 'parallel_geom' x 'parallel_sym' should be <= number of available CPUs
  # ForkManager package is used for parallelization within perl scripts
  'parallel_geom',   1,  # number of geometries to be invoked in the same time
  'parallel_symm',   1,  # number of symmetries to be invoked in the same time (values <= 0 trigger alternative workflow using mpi-scatci)

  # For very special run to get bound states of e + target without running scattering
  #   gaustail is skipped, thus integrals are not limited to the R-matrix sphere and scatci calculates regular eigenstates, not R-matrix poles
  # set this option to 1 if you want to calculate bound states
  'bound',           0,

  # EXPERT SETTING:
  # Force use of cdenprop for target calculation.
  'use_cdenprop', 0,

  # EXPERT SETTING:
  # Set a value for igh. ight = for target-scatci, ighs = for scattering-scatci.
  # 2 = Auto select (default), -1 = Arpack, 0 = Davidson, 1 = Givens-Householder
  'ight', 2,
  'ighs', 2,

  # EXPERT SETTING:
  # If set to 1 the amplitudes and channel data saved from a previous run (using the options 'save_channels' and 'save_rmat_amp')
  # will be used instead of running SWINTERF to generate them. This is useful in case the inner region data (e.g.fort.25) have
  # not been saved only a rerun of the outer region is needed e.g. with a different energy grid or using a different set of programs.
  'use_saved_ramps', 0,

  # EXPERT SETTING: in case of scattering calculations a number of programs can be run following the determination of K-matrices.
  # The list of the programs in the order they're supposed to be run is given below. Normally, you don't need to modify it.
  'run_eigenp',      1,  # calculate eigenphase sums for all symmetries
  'run_tmatrx',      1,  # calculate T-matrices for all symmetries
  'run_ixsecs',      1,  # calculate cross sections for all symmetries
  'run_reson',       1,  # calculate resonance fits for all symmetries
  'run_time_delay',  0,  # calculate time delays (requires a program that is not part of UKRmol+)

  # Photoionization
  'dipelm_smooth',   0,  # set: (0) for raw data, (1) for smoothed data or (2) for both.

  # For 'debugging' only (or maybe it can be useful to rerun outer region only, but don't forget to turn off cleaning)
  # if run_only is empty, then all programs will be executed
  # otherwise it is assumed that outputs from the previous run exist
  #           and only specified programs will be executed
  # 'only',           "target-scatci_integrals|target-congen|target-scatci|target-denprop|scattering-congen|scattering-scatci|scattering-swinterf|scattering-rsolve|scattering-eigenp|scattering-tmatrx|scattering-ixsecs|scattering-reson",  # e.g. "scattering-congen|scattering-scatci|scattering-outer" to run only last part of the scastering calculation

);
