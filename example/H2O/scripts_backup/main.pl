use Cwd;         # for changing the working directory
use File::Basename ();
use File::Spec ();

my $libpath;

BEGIN {
  # evaluate all scripts present on the command line before parsing the rest of this file
  foreach my $arg (@ARGV) {
    local $/ = undef;
    open FILE, "$arg" or die "Couldn't open file: $!";
    $content = <FILE>;
    close FILE;
    eval $content;
  }

  # if a library path is specified in %dirs, use it instead of the default "../lib"
  $libpath = File::Spec->catdir(File::Basename::dirname(Cwd::abs_path __FILE__), '../lib');
  if ($dirs{'libs'} ne "") { $libpath = $dirs{'libs'}; }
}

use lib $libpath;
use ForkManager; # for parallelization of geometries and spin-symmetries
use ukrmollib;   # my package of subroutines to handle geometries, inputs, outputs etc.

our $icform = "'F'";          # formatting of the channel information file (fort.10)
our $irform = "'F'";          # formatting of the boundary amplitudes file (fort.21)
our $ikform = "'F'";          # formatting of the K-matrix file (fort.9SM)

$run{'parallel_diag'} = 0; #Force use of scatci for build and diagonalization (use of hamdiag for diagonalization has not been implemented yet).

my %data = ( # this hash array together with option hash arrays above will be passed to subroutines calling programs
             # various subroutines which read outputs update data in this array

  # Information about molecular orbitals (read from scf output if not specfied otherwise):
  'orbitals',   {'energies',    {}, # energies of molecular orbitals as hash array, e.g. "2.1" -> -1.23456 is energy of the 1st state with symmetry 2
                 'target_all',  [], # from target calculations:    number of all target MOs in each irreducible representation (IR)
                 'target_sum',   0, # from target calculations:    sum of all target MOs from each irreducible representation (IR)
                 'occupied',    [], # from target calculations:    number of occupied target MOs in each IR
                 'frozen',      [], # user specified:              number of frozen target MOs in each IR
                 'active',      [], # user specified:              number of active target MOs in each IR
                 'reference',   [], # user specified:              number of MOs in each IR, which will be used for searching reference orbitals
                 'target',      [], # 'frozen' + 'active':         number of target MOs in each IR used to describe target states
                 'virtual',     [], # user specified:              number of virtual MOs in each IR, determined from energies
                 'target_used', [], # 'target' + 'virtual':        number of target MOs in each IR, which will be used in scattering calculation
                 'all',         [], # from scatci_integrals:       number of all basis functions in each IR
                 'continuum',   [], # 'all' - 'target_all'         number of continuum orbitals
                 'cont_used',   [], # from scatci_integrals:       number of continuum orbitals used for scattering after orthogonalization
                 'all_used',    [], # 'target_used' + 'cont_used'
                 'ifock',       []},# occupancy of orbitals in occ, 1 = full, 2 = half full (see swscf documentation, iocc and ifock)

  # Information about target states
  'target',     {'charge',         0,  # target charge
                 'hf_energy',      0,  # Hartree-Fock energy of target
                 'hf_symmetry',    0,  # Hartree-Fock target state symmetry (needed when it is doublet)
                 'nfull_occ',      0,  # number of fully occupied orbitals in SCF calculation
                 'nhalf_occ',      0,  # number of half  occupied orbitals in SCF calculation
                 'noccupied',      0,  # total number of occupied orbitals in SCF calculation
                 'spin',           "", # current spin multiplicity of the target state
                 'symmetry',       0,  # current symmetry of the target state
                 'states',         {}, # energies of all target states as read from target.scatci output with keys of the form "spin.sym.n" such as "1.0.n" or "3.2.n" where n numbers states of the same spin-symmetry
                 'ordered_states', [], # array contains strings "spin.sym.n" (see above) describing target states ordered by energy
                 'ground_state',   "", # contains string "spin.sym.1" as above of the target ground state
                 'spinsym_order',  [], # order in which different target state spin-symmetries were calculated (each element contains string "spin.symmetry")
                 'used_tgt_states',{}, # number of target state for each spin-symmetry which are used in scattering calculations (e.g. "1.0" -> 2, "1.1" -> 1, ...)
                 # auxiliary variables to use for scatci and denprop inputs (generated automatically according to target states requirements)
                 'ntgt',           0,  # number of target spin-symmetries used in scattering calculations (ntgt -> denprop, ntgsym -> scattering.scatci)
                 'ntgtl',          [], # number of states used in scattering calculations in each symmetry (ntgtl -> denprop, numtgt -> scattering.scatci)
                 'mcont',          [], # symmetry of the continuum orbitals associated with each target spin-symmetry (see MCONT in scattering.scatci), but changed in make_scatci_input according to symmetry of the final state if necessary)
                 # auxiliary variables to use for outer input
                 'idtarg',         [], # mapping of target states to energy order (idtarg for swinterf in outer input)
                 # auxiliary arrays for energies etc for all geometries
                 'states_all_geom',{}, # arrays of energies of all target states for all geometries as read from target.scatci output with keys of the form "spin.sym.n" such as "1.0.n" or "3.2.n" where n numbers states of the same spin-symmetry
                                       # e.g. $data{'target'}->{'states_all_geom'}->{'1.0.1'} = [E_geom1, E_geom2, ...]
                 'dipole_all_geom',{}, # arrays of dipole moments of the ground state for all geometries as read from the file borndata during read_denprop_output
                                       # with keys 'x', 'y' and 'z' e.g. $data{'target'}->{'dipole_all_geom'}->{'x'} = [dmx_geom1, dmx_geom2, ...]
                },

  # Information about inner states
  'scattering', {'basis',    "", # string such as "q0.r10.L6" determined by which_continuum_basis_file_to_use() below
                 'spin',     "", # current spin multiplicity of the scattering state
                 'symmetry', 0,  # current symmetry of the scattering state
                 'cont_csf', 0,  # number of N+1 configurations containing the continuum orbitals (LCONT)
                 'states_all_geom',{}, # arrays of energies of all R-matrix poles for all geometries as read from scattering.scatci output with keys of the form "spin.sym.n" such as "1.0.n" or "3.2.n" where n numbers states of the same spin-symmetry
                                       # e.g. $data{'scattering'}->{'states_all_geom'}->{'1.0.1'} = [E_geom1, E_geom2, ...]
                 'resonance_position',   {}, # arrays of positions of resonance for all geometries as read from scattering.outer output with keys of the form "spin.sym"
                 'resonance_width',      {}, # arrays of widths    of resonance for all geometries as read from scattering.outer output with keys of the form "spin.sym"
                },

  # determined automatically (whatever is set will be overwritten)
  'geometries',    [], # array of all geometries
  'geom_labels',   "", # will be labels for geometries used at the first line in files for potential curves etc.
  'task',          "", # will be the current regime, "target" or "scattering"
  'program',       "", # will be the current program to be executed
  'inputfile',     "", # will be the current input file
  'outputfile',    "", # will be the current output file
  'logfile',       "", # will be the current log file
  'igeom',         0,  # index of current geometry to work with
  'geom',          0,  # will be set in the geometry loop as a reference to a given geometry
  'nir',           0,  # number of irreducible representation
  'orbs_ok',       1,  # if setting of frozen and active orbitals is inconcistent this is set to 0
  'scf_ok',        0,  # set to 1, if iocc is consistent with result of SCF calculations
  'no_scatci',     0,  # set to 1, if there is no configuration generated during congen run -> skipping scatci

);

# ================== beginning of processing =================

# ----------------------- preliminaries ----------------------

# collection of hash arrays for simple passing of all variables to subroutines
my %parameters = ('model', \%model,
                  'run',   \%run,
                  'dirs',  \%dirs,
                  'data',  \%data);

# Set number of irreducible representations - used as abbreviation
$data{'nir'} = scalar @{$irred_repr{$model{'symmetry'}}};

# Add current working directory to directories specified by user
# (necessary because we change the working directory later)
$dirs{'cwd'} = getcwd();
if ($sys eq "win") { $dirs{'cwd'} =~ s/\//\\/g; } # For some reason the path on windows is returned with / instead of \

&add_cwd_to_dirs(\%dirs);

# If MASSCF is used get number of closed and active orbitals from the qchem MAS
# Note that $model{'nvitual'} is still used to specify the number of virtuals
# kept.
if ($model{'use_MASSCF'} > 0){
  my %orbs_per_space_type;
  if ($model{'qchem_MAS'}){
    %orbs_per_space_type = MultiSpace::orbs_per_subspace_type($model{'qchem_MAS'});
  } else {
    %orbs_per_space_type = MultiSpace::orbs_per_subspace_type($model{'MAS'});
  }
  $model{'nfrozen'} = $orbs_per_space_type{'closed'};
  $model{'nactive'} = $orbs_per_space_type{'active'};
}

my $nPCO=0;
my $l=0;
foreach( @{$model{'num_PCOs'}}){ $nPCO += $_ * (2*$l+1) ;$l++;}

my $PCO_name = "";
if ( $model{'use_PCO'} ) { $PCO_name =  ".".$nPCO."nPCO.".@model{'maxl_PCO'}."maxl_PCO"; }

if ( $model{'positron_flag'} eq 0 ) {$elpos = "e-";}
else {$elpos = "e+";}


# Set and make the model directory
if ($model{'directory'} eq "") {
  # if the model is "CAS", 'nactive' = -1 and 'nvirtual' = 0 then FCI is performed
  if ($model{'model'} eq "CAS" && $model{'nactive'} < 0 && $model{'nvirtual'} == 0) {
    $dirs{'model'} = "$dirs{'output'}${bs}$model{'molecule'}${bs}$model{'basis'}.FCI.$model{'nfrozen'}frozen$PCO_name.$model{'ntarget_states_used'}states.r$model{'rmatrix_radius'}.$model{'symmetry'}";
  }
  else {
    # example: output/H2/cc-pVTZ.CAS.e-.0frozen.12active.10virtual.4states.r13.D2h
    $dirs{'model'} = "$dirs{'output'}${bs}$model{'molecule'}${bs}$model{'basis'}.$model{'model'}.$elpos.$model{'nfrozen'}frozen.$model{'nactive'}active.$model{'nvirtual'}virtual$PCO_name.$model{'ntarget_states_used'}states.r$model{'rmatrix_radius'}.$model{'symmetry'}";
  }
}
else {
  $dirs{'model'} = "$dirs{'output'}${bs}$model{'directory'}";
}
if ($geometry{'suffix'} ne "") { $dirs{'model'} .= ".$geometry{'suffix'}"; }
if ($run{'suffix'}      ne "") { $dirs{'model'} .= ".$run{'suffix'}"; }
if ($run{'bound'}       == 1)  { $dirs{'model'} .= '.bound'; }
$dirs{'model'} =~ s/\.\././g; # get rid of double dots if user specified suffixes starting with a dot

&make_dir($dirs{'model'});

# Set scripts dir and copy used scripts into it
# so one can later easily find out settings for results
$dirs{'scripts'} = "$dirs{'model'}${bs}scripts_backup";
&make_dir($dirs{'scripts'});
foreach my $file ($0, @ARGV) { # $0 is the main script
  &copy_file($file, "$dirs{'scripts'}${bs}$file");
}
foreach my $file (@{$run{'add_files_to_backup'}}) { # Copy the libraries
  &copy_file($INC{$file}, "$dirs{'scripts'}${bs}$file");
}

# Set log dir and file and print basic information about run
$dirs{'logs'} = "$dirs{'model'}${bs}logs";
&make_dir($dirs{'logs'});
$data{'logfile'} = "$dirs{'logs'}${bs}main.log";

# Make directories for collected scattering data
$dirs{'data'} = "$dirs{'model'}${bs}collected_scattering_data";
&make_dir($dirs{'data'});
&make_dir("$dirs{'data'}${bs}xsec");     # directory for cross sections
&make_dir("$dirs{'data'}${bs}eigenph");  # directory for cross sections
if ($run{'save_channels'}) { &make_dir("$dirs{'data'}${bs}channels"); } # directory for channel data
if ($run{'save_rmat_amp'}) { &make_dir("$dirs{'data'}${bs}rmat_amplitudes"); } # directory for channel data
if ($run{'save_Kmatrix'}) { &make_dir("$dirs{'data'}${bs}K-matrices"); } # directory for K-matrices
if ($run{'save_Tmatrix'}) { &make_dir("$dirs{'data'}${bs}T-matrices"); } # directory for T-matrices

# Determine the target charge
$data{'target'}->{'charge'} = -$model{'nelectrons'};
foreach (@{$model{'atoms'}}) { $data{'target'}->{'charge'} += $atomic_number{$_}};

# determine which continuum basis file will be used, saved into $data{'scattering'}->{'basis'}
if ($run{'scattering'} == 1) { &which_continuum_basis_file_to_use(\%parameters); }

# print basic information about run
&print_info("Current working directory: $dirs{'cwd'}\n\n", \%parameters);
&print_info("Running UK R-matrix codes for e + $model{'molecule'}\n", \%parameters);
&print_info("Target charge q = $data{'target'}->{'charge'}\n", \%parameters);
if ($run{'scattering'} == 1) { &print_info("  the continuum basis swmol3.continuum.$data{'scattering'}->{'basis'} will be loaded.\n", \%parameters); }
&print_info("\nOutput in $dirs{'model'}\n\n", \%parameters);
&print_info("Symmetry $model{'symmetry'} with $data{'nir'} IRs for all geometries.\n\n", \%parameters);

# Check consistency of settings for orbitals and states
&check_settings(\%parameters);

# These directories will be created and used in directories for each geometry
$dirs{'inputs'}  = "inputs";      # Directory where all inputs will be stored  - e.g. as target.molpro.inp
$dirs{'outputs'} = "outputs";     # Directory where all outputs will be stored - e.g. as target.molpro.out

# Next we get all geometries using the subroutine generate_geomatries
# which also creates the file 'geometries' and directories for each geometry
# For each geometry it returns hash array { 'geometry', "string of distances and angles as in the file geometries",
#                                           'dir', directory,
#                                           'symmetry', sym,
#                                           'natype', value,
#                                           'atoms', [ ['A',x,y,z for atom 1], ..., ['B',x,y,z for last atom] ] }
# where 'sym' is one of (D2h, C2v, C2, Cs, C2h, D2, Ci) and 'A' or 'B' are symbols of atoms

push( @{$data{'geometries'}}, &generate_geometries(\%parameters, \%geometry));
print_info("\nStarting main loop over geometries.\n", \%parameters);
print_info("Logging continues into files corresponding to each geometry.\n", \%parameters);

# ================ main loop over geometries =================

if ($run{'parallel_geom'} > 1) { $pgeom = Parallel::ForkManager->new($run{'parallel_geom'}); }
foreach $r_geom (@{$data{'geometries'}}) {
  $data{'igeom'} += 1;
  if ($data{'igeom'} >= $geometry{'start_at_geometry'} &&
      ($geometry{'stop_at_geometry'} == 0 || $data{'igeom'} <= $geometry{'stop_at_geometry'})) {

    if ($run{'parallel_geom'} > 1) { $pgeom->start and next; }
    $data{'geom'} = $r_geom;
    $dirs{'geom'} = $r_geom->{'dir'};
    chdir($dirs{'geom'});

    # Create directories for inputs and outputs
    &make_dir($dirs{'inputs'});
    &make_dir($dirs{'outputs'});

    $data{'logfile'} = "$dirs{'logs'}${bs}geom$data{'igeom'}.log";
    &print_info("\nGeometry \#$data{'igeom'}:\n===========\n", \%parameters);

    # =================== start of target ======================
    $data{'task'} = "target";
    &print_info("\nTarget calculations:\n\n", \%parameters);

    $data{'scf_ok'} = 0; # indicator whether SCF calculation converged

    my $qchem = $run{'psi4'} ? "psi4" : "molpro";

    if (exists($model{'use_MASSCF'}) && $model{'use_MASSCF'} == 1){
      # ORMASSCF
      $qchem = "molpro"
    } elsif (exists($model{'use_MASSCF'}) && $model{'use_MASSCF'} == 2){
      # GASSCF
      $qchem = "molcas"
    }

    # Run pre-defined models using the MAS approach if they have been implemented.
    # TODO: Implement PCO, CHF models and positron scattering using the MAS approach.
    if (!exists($model{'models_using_mas'}) || !defined($model{'models_using_mas'})){
      $model{'models_using_mas'} = {
      'SE',1, 'SEP',1, 'CAS',1, 'CAS-A',1, 'CAS-B',1, 'CAS-C',1, 'CHF',0, 'CHF-A',0, 'CHF-B',0
      };
    }
    if ($model{'positron_flag'} == 0 && $model{'use_PCO'} == 0) {
      if (!exists($model{'use_MAS'}) || !defined($model{'use_MAS'}) || $model{'use_MAS'} == 0){
        if($model{'models_using_mas'}->{$model{'model'}}){
          $model{'use_MAS'} = 1;
          $model{'MAS'} = undef;  # Ensure model is chosen by model type
          $model{'l2_MAS'} = undef; # Ensure model is chosen by model type
        }
      }
    }

    # if number of frozen and active orbitals in each IRs are not given
    # we have to run MOLPRO or MOLCAS twice to get number of orbitals from HF calculations first
    if ($data{'orbs_ok'} == 0 && $model{'model'} !~ m/^(SE|FCI) && $model{'orbitals'} eq "natural"/) {
      &print_info("We have to run ".uc($qchem)." twice, first for HF only, then CASSCF...\n", \%parameters);
      &run_code($qchem, \%parameters);
    }
    $data{'scf_ok'} = 1;
    if (exists($model{'use_MASSCF'}) && $model{'use_MASSCF'} > 0){
      &setup_mas("qchem", \%parameters);
      &print_info("Quantum chemistry MAS:\n", \%parameters)
      &print_info($data{'MAS'}->{'qchem'}->string_mas(), \%parameters)
    }
    &run_code($qchem, \%parameters);
    &run_code("scatci_integrals", \%parameters);
    &run_sub("make_symlink_ALWAYS", \%parameters, "moints", "fort.16");

    # from here we run congen and scatci for each spin state (multiplicity) and symmetry (IR) separately

      if (exists($model{'use_MAS'}) && $model{'use_MAS'} > 0){

        # Use predefined model names to select the model instead of setting
        # MAS by hand. Does not work for CHF models currently!
        if(!$model{'MAS'}){
          choose_model_by_name_mas(\%parameters);
        }

        # Set up the MAS and write the details to the log files.
        &setup_mas("rmat", \%parameters);
        &print_info("Target MAS:\n", \%parameters);
        &print_info($data{'MAS'}->{'target'}->string_mas(), \%parameters);
        if (exists($model{'l2_MAS'}) && defined($model{'l2_MAS'})){
          &print_info("L^2 MAS:\n", \%parameters);
          &print_info($data{'MAS'}->{'l2'}->string_mas(), \%parameters);
        }
      }

    if ($model{'model'} =~ /(^CHF)/) {
      &initialize_CHF_data(%parameters);
      &print_info("\nCoupled HF model -- overriding target states:\n", \%parameters);
      &print_info("ntarget_states_used = $model{'ntarget_states_used'}\n", \%parameters);
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'ntarget_states'}}) {
        &print_info("ntarget_states $statespin = [".join(",", @{$model{'ntarget_states'}->{$statespin}})."]\n", \%parameters);
      }
    }

    if ($model{'model'} eq 'CHF-A') {

      my @states = @{$data{'CHF'}->{'states'}};
      my $no_states = scalar @states;

      for $state (@states){
        my $id = $state->{'id'};
        my $sym = $state->{'sym'};
        my $mult = $state->{'multiplet'};
        $data{'CHF'}->{'current_state'} = $id;
        $data{'target'}->{'spin'} = $mult;
        $data{'target'}->{'symmetry'} = $sym;

        &run_code("congen", \%parameters);
        if ($data{'no_scatci'} == 1) {
          $data{'no_scatci'} = 0;
        }
        elsif ($run{'parallel_symm'} >= 1) {
          &run_code("scatci", \%parameters);
        }
      }

      if ($run{'parallel_symm'} <= 0) {
        &run_code("mpi-scatci", \%parameters);
      }

      for (my $i = 0; $i < $data{'nir'}; $i++) {
        foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'ntarget_states'}}) {
          $model{'ntarget_states'}->{$statespin}->[$i] = 0;
        }
      }
      $model{'ntarget_states_used'} = 0;

      for $state (@states){
        my $sym = $state->{'sym'};
        my $mult = $state->{'multiplet'};
        my $ncsfs = $state->{'no_csfs'};
        $model{'ntarget_states'}->{$mult}->[$sym] += $ncsfs;
        $model{'ntarget_states_used'} += $ncsfs;
      }

    } else {

      for (my $i = 0; $i < $data{'nir'}; $i++) {
        foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'ntarget_states'}}) {
          $data{'target'}->{'spin'} = $statespin;
          if ($model{'ntarget_states'}->{$statespin}->[$i] > 0) { # skip if target states of a given spinsymmetry are not required
            $data{'target'}->{'symmetry'} = $i;
            &run_code("congen", \%parameters);
            if ($data{'no_scatci'} == 1) {
              $data{'no_scatci'} = 0;
            }
            elsif ($run{'parallel_symm'} >= 1) {
              &run_code("scatci", \%parameters);
            }
          }
        } # foreach spin state (multiplicity)
      } # for each IRs

      if ($run{'parallel_symm'} <= 0) {
        &run_code("mpi-scatci", \%parameters);
      }

      if ($model{'model'} eq "CHF-B") {
        my $ntarg = 0;
        foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'ntarget_states'}}) {
          for (my $i = 0; $i < $data{'nir'}; $i++) {
            if ($model{'ntarget_states'}->{$statespin}->[$i] > 0){
              $ntarg += $model{'ntarget_states'}->{$statespin}->[$i];
            }
          }
        }
        $model{'ntarget_states_used'} = $ntarg;
      }
    }

    if ($model{'model'} =~ /(^CHF)/) {
      &print_info("\nCoupled HF model -- overriding target states:\n", \%parameters);
      &print_info("ntarget_states_used = $model{'ntarget_states_used'}\n", \%parameters);
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'ntarget_states'}}) {
        &print_info("ntarget_states $statespin = [".join(",", @{$model{'ntarget_states'}->{$statespin}})."]\n", \%parameters);
      }
    }
    if ($model{'use_MAS'} && $qchem eq 'molpro'){
      # Analysis of MAS target state energies.
      if (!($data{'scf_ok'} == 0 || $model{'model'} =~ m/^SE/ || $model{'orbitals'} eq "HF")) {
        print_mas_analysis(\%parameters);
      }
    }

    if ($run{'parallel_symm'} >= 1) {
      &run_code("denprop", \%parameters);
    }

    &run_system("$cp_cmd fort.$luprop prop.out", \%parameters);

    # ===================== end of target =======================

    # ================== start of scattering ====================
    if ($run{'scattering'} == 1 or $run{'photoionization'} == 1 or $run{'rmt_interface'} == 1) {
      $data{'task'} = "scattering";
      &print_info("\nScattering calculations with $model{'model'} model:\n\n", \%parameters);

      if ($run{'parallel_symm'} >= 1) {

        # in the case of photoionization we need a specific order, so turn parallel mode off
        if ($run{'photoionization'} == 1 or $run{'rmt_interface'} == 1) { $run{'parallel_symm'} = 1; }

        # diagonalization order; must start with the initial state in photoionization runs
        my @symm_order = 0 .. $data{'nir'}-1;
        if ($run{'photoionization'} == 1 and 1 <= $model{'initialsym'} and $model{'initialsym'} <= $data{'nir'}) {
          $symm_order[0] = $symm_order[$model{'initialsym'} - 1];
          $symm_order[$model{'initialsym'} - 1] = 0;
        }

        # from here we run congen and scatci for each spin state (multiplicity) and symmetry (IR) separately
        if ($run{'parallel_symm'} > 1) { $psymm = Parallel::ForkManager->new($run{'parallel_symm'}); }
        for (my $i = 0; $i < $data{'nir'}; $i++) {
          foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'scattering_states'}}) {
            if ($run{'parallel_symm'} > 1) { $psymm->start and next; }
            $data{'scattering'}->{'spin'} = $statespin;
            # skip spin-symmetries which are not required
            if ($model{'scattering_states'}->{$statespin}->[$symm_order[$i]] == 1) {
              $data{'scattering'}->{'symmetry'} = $symm_order[$i];
              my $spin_sym = "$statespin.$irred_repr{$model{'symmetry'}}->[$symm_order[$i]]";
              my $spin_sym_suffix = "$spin_multiplicity{$statespin}$symm_order[$i]";
              # each spin-symmetry calculation is done in a local directory
              &make_dir($spin_sym);
              chdir($spin_sym);
              &run_sub("make_symlink_ALWAYS", \%parameters, "..${bs}fort.$luamp", "fort.$luamp");
              &run_sub("make_symlink_ALWAYS", \%parameters, "..${bs}fort.$luprop", "fort.$luprop");
              if ($run{'rmt_interface'} == 1) {
                &run_sub("make_symlink_ALWAYS", \%parameters, "..${bs}fort.$lucsf_base$spin_sym_suffix", "fort.$lucsf_base$spin_sym_suffix");  # N+1 CSFs (to be written)
                &run_sub("make_symlink_ALWAYS", \%parameters, "..${bs}fort.$luci", "fort.$luci");  # scattering eigenstates (to be written)
              }
              &run_code("congen", \%parameters);
              if ($data{'no_scatci'} == 1) { # skip scatci if there is no configuration generated by congen
                $data{'no_scatci'} = 0;      # reset the flag
              }
              else {
                &run_sub("make_symlink_ALWAYS", \%parameters, "..${bs}fort.$lucitgt", "fort.$lucitgt");  # target states (ready)
                &run_code("scatci", \%parameters);
                if ($run{'parallel_diag'} == 1) { &run_code("hamdiag", \%parameters); }
                if ($run{'bound'} == 0) { # if bound states are required we do not run outer region codes
                  if ($run{'photoionization'} == 1) {
                    if ($i == 0) {
                      &run_system("$cp_cmd fort.$lucsf_base$spin_sym_suffix ..${bs}fort.${lucsf_base}00", \%parameters);
                      &run_system("$cp_cmd fort.$luci ..${bs}fort.${luci}0", \%parameters);
                    }
                    &run_sub("make_symlink_ALWAYS", \%parameters, "..${bs}fort.${luci}0", "fort.${luci}0");
                    &run_sub("make_symlink_ALWAYS", \%parameters, "..${bs}fort.${lucsf_base}00", "fort.${lucsf_base}00");
                    &run_code("cdenprop", \%parameters);
                    &run_code("swinterf", \%parameters);
                    &run_code("rsolve", \%parameters);
                    for (my $k = 0; $k < 3; $k++) {
                      if ($group_table[$symm_order[$i]]->[$symm_order[0]] == $coor_repr{$model{'symmetry'}}->[$k]) {
                        my $unit = $coor_repr{$model{'symmetry'}}->[$k];
                        &run_system("$cp_cmd fort.${lupwd_base}1$unit ..${bs}fort.${lupwd_base}1$unit", \%parameters);
                      }
                    }
                  } elsif ($run{'rmt_interface'} == 0) {
                    &run_code("swinterf", \%parameters);
                    &run_code("rsolve", \%parameters);
                    if ($run{"run_eigenp"}) { &run_code("eigenp", \%parameters); }
                    if ($run{"run_tmatrx"}) {
                      &run_code("tmatrx", \%parameters);
                      if ($run{"run_ixsecs"}) {
                        &run_code("ixsecs", \%parameters);
                      }
                    }
                    if ($run{"run_reson"} == 1) {
                      &run_code("reson", \%parameters);
                    }
                    if ($run{"run_time_delay"} == 1) {
                      &run_code("time-delay", \%parameters);
                    }
                  }
                  if ($run{"save_channels"}) { &run_system("$mv_cmd fort.$luchan $dirs{'data'}${bs}channels${bs}channels.geom$data{'igeom'}.$spin_sym", \%parameters); }
                  if ($run{"save_rmat_amp"}) { &run_system("$mv_cmd fort.$lurmt $dirs{'data'}${bs}rmat_amplitudes${bs}ramps.geom$data{'igeom'}.$spin_sym", \%parameters); }
                  if ($run{"save_Kmatrix"})  { &run_system("$mv_cmd fort.$lukmt_base$spin_sym_suffix $dirs{'data'}${bs}K-matrices${bs}K-matrix.geom$data{'igeom'}.$spin_sym", \%parameters); }
                  if ($run{"run_eigenp"} and $run{"save_eigenph"}) { &run_system("$cp_cmd fort.$lueig_base$spin_sym_suffix ..${bs}eigenph.$spin_sym", \%parameters); }
                  if ($run{"run_ixsecs"} and $run{"save_xsec"})    { &run_system("$cp_cmd fort.$luxsn_base$spin_sym_suffix ..${bs}xsec.$spin_sym", \%parameters); }
                  if ($run{"run_tmatrx"} and $run{"save_Tmatrix"}) { &run_system("$mv_cmd fort.$lutmt_base$spin_sym_suffix $dirs{'data'}${bs}T-matrices${bs}T-matrix.geom$data{'igeom'}.$spin_sym", \%parameters); }
                } # bound
              } # no_scatci
              if ($run{'clean'}) {
                system("$rm_cmd fort.*");
                if ($run{'bound'} == 0) { system("$rm_cmd reson_message"); }
                system("$rm_cmd log_file.*");
              }
              chdir($dirs{'geom'});
              if ($run{'clean'}) { rmdir($spin_sym); }
            } # run given spin-symmetry
            if ($run{'parallel_symm'} > 1) { $psymm->finish; }
          } # foreach spin state (multiplicity)
        } # for each IRs
        if ($run{'parallel_symm'} > 1) { $psymm->wait_all_children; }

        # --------------- photoionization cross section -------------------
        if ($run{'photoionization'} == 1) {
          if (%parameters{'run'}->{'dipelm_smooth'} == 0 or %parameters{'run'}->{'dipelm_smooth'} == 1){
            &run_code("dipelm", \%parameters);
          }
          else {
            # run both smooth and unsmoothed
            # start with raw/unsmoothed
            %parameters{'run'}->{'dipelm_smooth'} = 0;
            &run_code("dipelm", \%parameters);

            # rename unsmoothed data (prepend with us_ for unsmoothed)
            &run_system("$cp_cmd photo_beta_1c us_photo_beta_1c", \%parameters);
            &run_system("$cp_cmd photo_beta_2c us_photo_beta_2c", \%parameters);
            &run_system("$cp_cmd photo_beta_2l us_photo_beta_2l", \%parameters);
            &run_system("$cp_cmd photo_total_xsec us_photo_total_xsec", \%parameters);
            &run_system("$cp_cmd photo_xsec us_photo_xsec", \%parameters);

            # now run smoothed
            %parameters{'run'}->{'dipelm_smooth'} = 1;
            &run_code("dipelm", \%parameters);
          }
          my @direction = ('x', 'y', 'z');
          for (my $k = 0; $k < 3; $k++) {
            my $id = $coor_repr{$model{'symmetry'}}->[$k];
            &run_system("$cp_cmd fort.${lupwd_base}1$id pwdips-$direction[$k]", \%parameters);
          }
        }

        # --------------- rmt interface -------------------
        elsif ($run{'rmt_interface'} == 1) {
          &run_code("cdenprop_all", \%parameters);
          &run_code("rmt_interface", \%parameters);
        }

      } else {

        # alternative workflow using MPI-SCATCI: do all CONGEN runs first
        for (my $i = 0; $i < $data{'nir'}; $i++) {
          foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'scattering_states'}}) {
            $data{'scattering'}->{'spin'} = $statespin;
            if ($model{'scattering_states'}->{$statespin}->[$i] == 1) {
              $data{'scattering'}->{'symmetry'} = $i;
              my $spin_sym = "$statespin.$irred_repr{$model{'symmetry'}}->[$i]";
              &run_code("congen", \%parameters);
            }
          }
        }

        &run_code("mpi-scatci", \%parameters);

        if ($run{'save_channels'}) {
          &run_system("$mv_cmd fort.$luchan $dirs{'data'}${bs}channels${bs}channels.geom$data{'igeom'}.$spin_sym", \%parameters);
        }
        if ($run{'save_rmat_amp'}) {
          &run_system("$mv_cmd fort.$lurmt $dirs{'data'}${bs}rmat_amplitudes${bs}ramps.geom$data{'igeom'}.$spin_sym", \%parameters);
        }

        for (my $i = 0; $i < $data{'nir'}; $i++) {
          foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'scattering_states'}}) {
            $data{'scattering'}->{'spin'} = $statespin;
            if ($model{'scattering_states'}->{$statespin}->[$i] == 1) {
              $data{'scattering'}->{'symmetry'} = $i;
              my $spin_sym = "$statespin.$irred_repr{$model{'symmetry'}}->[$i]";
              my $spin_sym_suffix = "$spin_multiplicity{$statespin}$i";

              # if bound states are required we do not run outer region codes
              if ($run{'bound'} == 0) {
                if ($run{'photoionization'} == 1) {
                  &run_code("cdenprop", \%parameters);
                  &run_system("$mv_cmd dyson_orbitals.molden0 dyson_orbitals.$spin_sym.molden", \%parameters);
                  &run_system("$mv_cmd dyson_orbitals.ukrmolp0 dyson_orbitals.$spin_sym.ukrmolp", \%parameters);
                  if ($run{'clean'} == 1) {
                    &run_system("$rm_cmd target.phases.data");
                  }
                }
                if ($run{'scattering'} == 1 or $run{'photoionization'} == 1) {
                  &run_code("rsolve", \%parameters);
                  if ($run{'save_Kmatrix'}) {
                    &run_system("$cp_cmd fort.$lukmt_base$spin_sym_suffix $dirs{'data'}${bs}K-matrices${bs}K-matrix.geom$data{'igeom'}.$spin_sym", \%parameters);
                  }
                }
                if ($run{'scattering'} == 1 and $run{'run_eigenp'} == 1) {
                  &run_code("eigenp", \%parameters);
                  if ($run{'save_eigenph'}) {
                    &run_system("$cp_cmd fort.$lueig_base$spin_sym_suffix eigenph.$spin_sym", \%parameters);
                  }
                }
                if ($run{'scattering'} == 1 and $run{'run_tmatrx'} == 1) {
                  &run_code("tmatrx", \%parameters);
                  if ($run{'save_Tmatrix'}) {
                    &run_system("$mv_cmd fort.$lutmt_base$spin_sym_suffix $dirs{'data'}${bs}T-matrices${bs}T-matrix.geom$data{'igeom'}.$spin_sym", \%parameters);
                  }
                }
                if ($run{'scattering'} == 1 and $run{'run_ixsecs'} == 1) {
                  &run_code("ixsecs", \%parameters);
                  if ($run{'save_xsec'}) {
                    &run_system("$cp_cmd fort.$luxsn_base$spin_sym_suffix xsec.$spin_sym", \%parameters);
                  }
                }
                if ($run{'scattering'} == 1 and $run{'run_reson'} == 1) {
                  &run_code("reson", \%parameters);
                  if ($run{'clean'} == 1) {
                    &run_system("$rm_cmd reson_message", \%parameters);
                  }
                }
                if ($run{"run_time_delay"} == 1) {
                  &run_code("time-delay", \%parameters);
                }
              } # bound

            } # run given spin-symmetry
          } # foreach spin state (multiplicity)
        } # for each IRs

        # --------------- photoionization cross section -------------------
        if ($run{'photoionization'} == 1) {
          if (%parameters{'run'}->{'dipelm_smooth'} == 0 or %parameters{'run'}->{'dipelm_smooth'} == 1) {
            &run_code("dipelm", \%parameters);
          } else {
            # run both smooth and unsmoothed
            # start with raw/unsmoothed
            %parameters{'run'}->{'dipelm_smooth'} = 0;
            &run_code("dipelm", \%parameters);

            # rename unsmoothed data (prepend with us_ for unsmoothed)
            &run_system("$cp_cmd photo_beta_1c us_photo_beta_1c", \%parameters);
            &run_system("$cp_cmd photo_beta_2c us_photo_beta_2c", \%parameters);
            &run_system("$cp_cmd photo_beta_2l us_photo_beta_2l", \%parameters);
            &run_system("$cp_cmd photo_total_xsec us_photo_total_xsec", \%parameters);
            &run_system("$cp_cmd photo_xsec us_photo_xsec", \%parameters);

            # now run smoothed
            %parameters{'run'}->{'dipelm_smooth'} = 1;
            &run_code("dipelm", \%parameters);
          }
          my @direction = ('x', 'y', 'z');
          for (my $k = 0; $k < 3; $k++) {
            my $id = $coor_repr{$model{'symmetry'}}->[$k];
            &run_system("$cp_cmd fort.${lupwd_base}1$id pwdips-$direction[$k]", \%parameters);
          }
        }

      } # alternative workflow using MPI-SCATCI

      # --------------- gathering scattering data -------------------
      if ($run{'gather_data'} == 1) {
        if ($run{'bound'} == 0) {
          if ($run{'run_eigenp'} == 1) {
            &gather_eigenphases("eigenph", \%parameters);
          }
          if ($run{'run_ixsecs'} == 1) {
            &convert_and_gather_cross_sections("xsec", \%parameters);
          }
        }
      }
    }
    # =================== end of scattering =====================


    # ------------- cleaning before next geometry ---------------
    foreach my $key (%{$data{'orbitals'}}) {
      $data{'orbitals'}->{$key} = [];
    }
    $data{'orbitals'}->{'energies'} = {};
    $data{'target'}->{'states'}          = {};
    $data{'target'}->{'ordered_states'}  = [];
    $data{'target'}->{'spinsym_order'}   = [];
    $data{'target'}->{'used_tgt_states'} = {};
    $data{'target'}->{'ntgt'}            = 0;
    $data{'target'}->{'ntgtl'}           = [];
    $data{'target'}->{'mcont'}           = [];
    $data{'target'}->{'idtarg'}          = [];

    if ($run{'clean'} == 1) {
        system("$rm_cmd fort.*");
        system("$rm_cmd log_file.*");
        if ($run{'psi4'} == 1) {
          system("$rm_cmd timer.dat");
        }
    }
    if ($run{'remove_moints'} == 1) {
        system("$rm_cmd moints");
    }
    chdir($dirs{'cwd'});
    if ($run{'parallel_geom'} > 1) { $pgeom->finish; }
  } # if ($data{'igeom'} >= $geometry{'start_at_geometry'})

} # foreach $r_geom
if ($run{'parallel_geom'} > 1) { $pgeom->wait_all_children; }

# ================ collecting data for all geometries =================

$data{'logfile'} = "$dirs{'logs'}${bs}main.log";
&print_info("\nStart of collecting data for all geometries ...\n", \%parameters);

# if the script runs in parallel mode we have to collect data separately
# because everything stored in %data during parallel run was lost
# and because of possible cleaning of inputs and output it is done always
$data{'igeom'} = 0;
$run{'only'} = "nothing"; # to skip running codes when collecting data
foreach $r_geom (@{$data{'geometries'}}) {
  $data{'igeom'} += 1;
  $data{'geom'} = $r_geom;
  $dirs{'geom'} = $r_geom->{'dir'};
  chdir($dirs{'geom'});

  $data{'task'} = "target";
  # from here we read scatci output for each spin state (multiplicity) and symmetry (IR) separately
  for (my $i = 0; $i < $data{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'ntarget_states'}}) {
      $data{'target'}->{'spin'} = $statespin;
      if ($model{'ntarget_states'}->{$statespin}->[$i] > 0) { # skip if target states of a given spinsymmetry are not required
        $data{'target'}->{'symmetry'} = $i;
        if ($run{'parallel_symm'} >= 1) {
          &run_code("scatci", \%parameters);
        }
      }
    } # foreach spin state (multiplicity)
  } # for each IRs

  if ($run{'parallel_symm'} <= 0) {
    if ($run{'clean'}) { &run_system("$cp_cmd prop.out fort.$luprop"); }
    &run_code("mpi-scatci", \%parameters);
    if ($run{'clean'}) { &run_system("$rm_cmd fort.$luprop"); }
  }

  if ($run{'scattering'} == 1) {
    $data{'task'} = "scattering";
    if ($run{'parallel_symm'} <= 0) {
      &run_code("mpi-scatci", \%parameters);
    }
    # from here we read scatci output for each spin state (multiplicity) and symmetry (IR) separately
    for (my $i = 0; $i < $data{'nir'}; $i++) {
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$model{'scattering_states'}}) {
        $data{'scattering'}->{'spin'} = $statespin;
      # skip spin-symmetries which were not required
        if ($model{'scattering_states'}->{$statespin}->[$i] == 1) {
          $data{'scattering'}->{'symmetry'} = $i;
          my $spin_sym = "$statespin.$irred_repr{$model{'symmetry'}}->[$i]";
          if ($run{'parallel_symm'} >= 1) {
            &run_code("scatci", \%parameters);
          }
          if ($run{'bound'} == 0) { # if bound states were required we did not run outer region codes
            if ($run{'scattering'} == 1 and $run{'run_reson'} == 1) {
              &run_code("reson", \%parameters);
            }
          }
        }
      } # foreach spin state (multiplicity)
    } # for each IRs
  }

  # ------------- cleaning of inputs and outputs ---------------
  if ($run{'keep_inputs'} == 0) {
    system("$rm_cmd $dirs{'inputs'}${bs}*.inp");
    rmdir($dirs{'inputs'});
  }
  if ($run{'keep_outputs'} == 0) {
    system("$rm_cmd $dirs{'outputs'}${bs}*.out");
    rmdir($dirs{'outputs'});
  }
  chdir($dirs{'cwd'});

} # foreach $r_geom

if ($run{'gather_data'} == 1) {
  &save_target_energies("target.energies", \%parameters);
  &save_dipole_moments("target.dipole.moments", \%parameters);
  if ($run{'scattering'} == 1) {
    &save_rmatrix_energies("Rmatrix.energies", \%parameters, 5);
    if ($run{'bound'} == 0) {
      &save_resonance_positions_and_widths("resonance.positions.and.widths", \%parameters);
      if ($run{'run_eigenp'}) {
        &make_gnuplot_files_for_eigenphase_sums("eigenph", \%parameters);
      }
      if ($run{'run_ixsecs'}) {
        &make_gnuplot_files_for_cross_sections("xsec", \%parameters);
      }
    }
  }
}
