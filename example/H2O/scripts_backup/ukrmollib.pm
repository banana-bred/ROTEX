use dirfile;
use Storable 'dclone';
use List::Util qw(max min sum0);
use MultiSpace;

# Settings of Molpro input for different symmetries
our %symmetries = (
  'D2h',   [ 3, "'X', 'Y', 'Z'" ],
  'C2v',   [ 2, "'X', 'Y'" ],
  'C2h',   [ 2, "'XY', 'Z'" ],
  'D2',    [ 2, "'XZ', 'YZ'" ],
  'C2',    [ 1, "'XY'" ],
  'Cs',    [ 1, "'X'" ],
  'Ci',    [ 1, "'XYZ'" ],
  'C1',    [ 0, "''" ],
);

# Order of irreducible representations in Molpro
our %irred_repr = (
  'D2h',   [ "Ag", "B3u", "B2u", "B1g", "B1u", "B2g", "B3g", "Au" ],
  'C2v',   [ "A1", "B1", "B2", "A2" ],
  'C2h',   [ "Ag", "Au", "Bu", "Bg" ],
  'D2',    [ "A", "B3", "B2", "B1" ],
  'C2',    [ "A", "B" ],
  'Cs',    [ "Ap", "App" ], # p stands for prime ', there is a problem on linux machines with '
  'Ci',    [ "Ag", "Au" ],
  'C1',    [ "A" ],
);

# Order of irreducible representations in Psi4 in terms of the Molpro-ordered irreducible representations
our %irred_psi4 = (
  'D2h',   [ 0, 3, 5, 6, 7, 4, 2, 1 ],
  'C2v',   [ 0, 3, 1, 2 ],
  'C2h',   [ 0, 3, 1, 2 ],
  'D2',    [ 0, 3, 2, 1 ],
  'C2',    [ 0, 1 ],
  'Cs',    [ 0, 1 ],
  'Ci',    [ 0, 1 ],
  'C1',    [ 0 ],
);

# Group multiplication table corresponding to Molpro ordering
# For smaller groups than D2h only part of this table (e.g. 4x4 for C2v or 2x2 for C2) is used
our @group_table = (
  [ 0, 1, 2, 3, 4, 5, 6, 7 ],
  [ 1, 0, 3, 2, 5, 4, 7, 6 ],
  [ 2, 3, 0, 1, 6, 7, 4, 5 ],
  [ 3, 2, 1, 0, 7, 6, 5, 4 ],
  [ 4, 5, 6, 7, 0, 1, 2, 3 ],
  [ 5, 4, 7, 6, 1, 0, 3, 2 ],
  [ 6, 7, 4, 5, 2, 3, 0, 1 ],
  [ 7, 6, 5, 4, 3, 2, 1, 0 ]
);

# Irreducible representations for X, Y, Z coordinates
our %coor_repr = (
  'D2h',   [ 1, 2, 4 ],  # X = B3u, Y = B2u, Z = B1u
  'C2v',   [ 1, 2, 0 ],  # X = B1,  Y = B2,  Z = A1
  'C2h',   [ 1, 1, 2 ],  # X = Bu,  Y = Bu,  Z = Au
  'D2',    [ 1, 2, 3 ],  # X = B3,  Y = B2,  Z = B1
  'C2',    [ 1, 1, 0 ],  # X = B,   Y = B,   Z = A
  'Cs',    [ 0, 0, 1 ],  # X = Ap,  Y = Ap,  Z = App
  'Ci',    [ 1, 1, 1 ],  # X = Au,  Y = AU,  Z = Au
  'C1',    [ 0, 0, 0 ]   # X = A,   Y = A,   Z = A
);

our %spin_multiplicity = (
  'singlet', 1,
  'doublet', 2,
  'triplet', 3,
  'quartet', 4,
  'quintet', 5,
  'sextet', 6,
  'septet', 7,
  'octet', 8,
  'nonet', 9,
);

# Relative atomic masses - used to determine the center of mass
# from http://www.nist.gov/pml/data/comp.cfm
# only the most abundant isotopes
our %mass = (
  'H' ,   1.007825,
  'D' ,   2.014102,
  'T' ,   3.016049,
  'He',   4.002603, # He4 100.0% (almost)
  'Li',   7.016005, # Li7  92.4%
  'Be',   9.012182,
  'B' ,  11.009305, # B11  80.1%
#  'C' ,  12.000000, # C12  98.9%
  'C' ,  12.0110,   # Molpro
  'N' ,  14.003074, # N14  99.6%
  'O' ,  15.9994,   # Molpro
#  'O' ,  15.994915, # O16  99.8%
  'F' ,  18.998403,
  'Ne',  19.992440, # Ne20 90.5%
  'Na',  22.989769,
  'Mg',  23.985042, # Mg24 79.0%
  'Al',  26.981539,
  'Si',  27.976927, # Si28 92.2%
  'P' ,  30.973762,
  'S' ,  31.972071, # S32  95.0%
  'Cl',  34.968853  # Cl35 75.8%
);

# Atomic numbers - used to automatic generation of the swmol3 bases
our %atomic_number = (
  'H' ,   1,
  'H@2.014101779' ,   1,
  'T' ,   1,
  'He',   2,
  'Li',   3,
  'Be',   4,
  'B' ,   5,
  'C' ,   6,
  'N' ,   7,
  'O' ,   8,
  'F' ,   9,
  'Ne',  10,
  'Na',  11,
  'Mg',  12,
  'Al',  13,
  'Si',  14,
  'P' ,  15,
  'S' ,  16,
  'Cl',  17,
  'Ar',  18,
  'K' ,  19,
  'Ca',  20,
  'Sc',  21,
  'Ti',  22,
  'V' ,  23,
  'Cr',  24,
  'Mn',  25,
  'Fe',  26,
  'Co',  27,
  'Ni',  28,
  'Cu',  29,
  'Zn',  30,
  'Ga',  31,
  'Ge',  32,
  'As',  33,
  'Se',  34,
  'Br',  35,
  'Kr',  36,
  'Rb',  37,
  'Sr',  38,
  'Y' ,  39,
  'Zr',  40,
  'Nb',  41,
  'Mo',  42,
  'Tc',  43,
  'Ru',  44,
  'Rh',  45,
  'Pd',  46,
  'Ag',  47,
  'Cd',  48,
  'In',  49,
  'Sn',  50,
  'Sb',  51,
  'Te',  52,
  'I' ,  53,
  'Xe',  54,
  'Cs',  55,
  'Ba',  56,
  'La',  57,
  'Ce',  58,
  'Pr',  59,
  'Nd',  60,
  'Pm',  61,
  'Sm',  62,
  'Eu',  63,
  'Gd',  64,
  'Tb',  65,
  'Dy',  66,
  'Ho',  67,
  'Er',  68,
  'Tm',  69,
  'Yb',  70,
  'Lu',  71,
  'Hf',  72,
  'Ta',  73,
  'W' ,  74,
  'Re',  75,
  'Os',  76,
  'Ir',  77,
  'Pt',  78,
  'Au',  79,
  'Hg',  80,
  'Tl',  81,
  'Pb',  82,
  'Bi',  83,
  'Po',  84,
  'At',  85,
  'Rn',  86,
  'Fr',  87,
  'Ra',  88,
  'Ac',  89,
  'Th',  90,
  'Pa',  91,
  'U' ,  92,
  'Np',  93,
  'Pu',  94,
  'Am',  95,
  'Cm',  96,
  'Bk',  97,
  'Cf',  98,
  'Es',  99,
  'Fm', 100,
  'Md', 101,
  'No', 102,
  'Lr', 103,
  'Rf', 104,
  'Db', 105,
  'Sg', 106,
  'Bh', 107,
  'Hs', 108,
  'Mt', 109,
  'Ds', 110,
  'Rg', 111
);

our $atomic_orbitals = "spdfghiklmno";

# --------------- hard-coded unit numbers and formats ----------------

# Fortran scratch unit numbers used in the scripts
our $luchan = 10;           # channel information from SWINTERF or MPI-SCATCI (fort.10)
our $luamp = 16;            # molecular integrals and one-electron properties from SCATCI_INTEGRALS (fort.16 -> moints)
our $lurmt = 21;            # boundary amplitudes from SWINTERF or MPI-SCATCI (fort.21)
our $luprop = 24;           # target state properties (fort.24)
our $luci = 25;             # scattering eigenstates (fort.25)
our $lucitgt = 26;          # target eigenstates (fort.26)
our $luitdip = 667;         # inner region eigenstates transition dipoles from CDENPROP_ALL (fort.667)

# Fortran scratch unit number prefixes (S = spin multiplicity, M = irreducible representation)
our $lupwd_base = 4;        # prefix for partial wave dipoles from RSOLVE or MPI_RSOLVE (fort.4SM)
our $luidip_base = 6;       # prefix for inner region dipoles from CDENPROP or MPI_RSOLVE (fort.6SM)
our $lucsf_base = 7;        # prefix for configuration state functions from CONGEN (fort.7SM)
our $luham_base = 8;        # prefix for Hamiltonian matrices from SCATCI and MPI-SCATCI (fort.8SM)
our $lukmt_base = 9;        # prefix for K-matrices from RSOLVE or MPI_RSOLVE (fort.9SM)
our $lueig_base = 10;       # prefix for eigenphases from EIGENP (fort.10SM)
our $lutmt_base = 11;       # prefix for T-matrices from TMATRX (fort.11SM)
our $luxsn_base = 12;       # prefix for cross sections from IXSECS (fort.12SM)
our $lures_base = 13;       # prefix for resonance fits (fort.13SM)

# Fortran scratch unit formatting
our $icform = "'U'";          # formatting of the channel information file (fort.10)
our $irform = "'U'";          # formatting of the boundary amplitudes file (fort.21)
our $ikform = "'U'";          # formatting of the K-matrix file (fort.9SM)
our $itform = "'U'";          # formatting of the T-matrix file (fort.11SM)

# Overide MultiSpace versions of these variables.
@MultiSpace::group_table = @group_table;
%MultiSpace::irred_repr = %irred_repr;

# --------------- subroutines for geometries ----------------

# One can use this subroutine to add more geometries by stretching particular bond(s)
# Several atoms which form a group can be moved at the same time in the direction of the specific bond
# An example of use can be a symmetric stretch of H2O. Assuming the equilibrium geometry was set manually
# in geometry.pl, one can call (also in geometry.pl after %geometry = (...);)
# &set_geometries_by_stretching_bonds(\%geometry,
#   {'ranges', [ [-0.4, -0.2, 0.2],
#                [ 0.1,  0.4, 0.1] ],
#    'shift', { '1-2', [2],
#               '1-3', [3] },
#    'gnuplot_prefix', "Delta O-H = "
#   }
# );
# to stretch both O-H bonds simultaneously.
# Antisymmetric stretch can be obtained by setting '3-1', [3] insted of '1-3', [3]
# Note: if only one range is specified there still must be two pairs of [] like this
#   'ranges', [ [-0.3, 0.3, 0.1] ],
sub set_geometries_by_stretching_bonds {
  my ($r_geometry, $r_param) = @_;

  # Collect values to use for stretching, add_range also sorts and deletes duplicates
  my @stretching_values = ();
  foreach $range (@{$r_param->{'ranges'}}) {
    &add_range(\@stretching_values, $range->[0], $range->[1], $range->[2]);
  }

  # Add another suffix to the directory with stretching parameters
  $r_geometry->{'suffix'} .= ".stretched.by".sprintf("%.2f", $stretching_values[0])."-".sprintf("%.2f", $stretching_values[-1]);

  # Copy the equilibrium geometry to a temporary array (it must be done using dclone from Storable)
  # and delete equilibrium geometry
  my $equilibrium_geometry = dclone($r_geometry->{'geometries'}->[0]);
  $r_geometry->{'geometries'} = [];

  # Loop over all stretching values
  foreach my $dr (@stretching_values) {

    # First, copy the equilibrium geometry to a new geometry and change description
    my $new_geometry = dclone($equilibrium_geometry);
    $new_geometry->{'description'} = sprintf("%8.5f", $dr);
    $new_geometry->{'gnuplot_desc'} = $r_param->{'gnuplot_prefix'}.sprintf("%6.2f", $dr);

    # Then shift given atoms in the direction of bonds specified in parameters as 'shift'
    foreach my $bond (keys %{$r_param->{'shift'}}) {
      my ($from, $to) = ($bond =~ /(\d+)\-(\d+)/);

      # Calculate the shift vector, we have to subtract 1 from $from and $to to get a correct array element
      my @n = map { $equilibrium_geometry->{'atoms'}->[$to - 1]->[$_] - $equilibrium_geometry->{'atoms'}->[$from - 1]->[$_] } 1..3; # direction
      my $norm = 0.0; $norm += $_*$_ for @n;
      @n = map { $dr * $_ / sqrt($norm) } @n; # adjust magnitude

      # Move all specified atoms for the given bond
      foreach my $atom (@{$r_param->{'shift'}->{$bond}}) {
        foreach my $i (1..3) {
          $new_geometry->{'atoms'}->[$atom - 1]->[$i] += $n[$i - 1];
        }
      }
    }
    push(@{$r_geometry->{'geometries'}}, $new_geometry);
  }
  return 1;
}

# One can use this subroutine to read geometries from a given file or several files
# where the geometries are given in the molden format.
# [Atoms] (Angs|AU)
# element_name number atomic_number x y z
# ...
# Here is an example of the equilibrium geometry of the water molecule given in Angstroms
# [Atoms] Angs
# O     1    8         0.0000000000        0.0000000000        0.0656240182
# H     2    1         0.0000000000        0.7572286165       -0.5208964435
# H     3    1         0.0000000000       -0.7572286165       -0.5208964435
# Several geometries can be specified in one file, each starting with [Atoms], all of them are read,
# or if $filename equals "*string*" then all files containg 'string' in their name are read
# and searched for geometries. ('string' would usually be 'molden', I guess.)
sub read_geometries_from_files_in_molden_format {
  my ($filename, $r_geometry) = @_;

  # Prepare @filelist with all files to be read
  my @filelist = ();
  if ($filename =~ m/^\*([^\*]+)\*/) {
    my $str = $1;
    foreach my $file (&dirfilelist(".")) {
      if ($file =~ m/$str/i) { push(@filelist, $file); }
    }
  }
  else {
    push(@filelist, $filename);
  }

  # Read each file and search for geometries
  $r_geometry->{'geometries'} = [];
  my $ngeom = 0;
  foreach my $file (@filelist) {
    if (open(INPUT, "$file")) {
      my $in_atoms = 0;
      while (my $line = <INPUT>) {
        chomp($line);
        if ($line =~ m/^\[Atoms\]\s+(Angs|AU)/i) {
          $in_atoms = 1;
          $r_geometry->{'length_unit'} = ($1 eq "Angs" ? 1 : 0);
          $ngeom++;
          push(@{$r_geometry->{'geometries'}},
            { 'description', "see geometries.molden",
              'gnuplot_desc', "Geometry $ngeom", # used in gnuplot files for keys
              'atoms', [ ]
            }
          );
        }
        elsif ($in_atoms == 1 && $line =~ m/^\s*([A-Za-z]+)\s+(\d+)\s+(\d+)/) {
          $line =~ s/^\s*//;
          $line =~ s/\s*$//;
          my @values = split(/\s+/, $line);
          push(@{$r_geometry->{'geometries'}->[$ngeom-1]->{'atoms'}},
            [ $values[0], $values[3], $values[4], $values[5] ]
          );
        }
        else {
          $in_atoms = 0;
        }
      }
      close(INPUT);
    }
    else {
      die "Error in read_geometries_from_files_in_molden_format():\nCould not open the file $file with geometries!\n";
    }
  }
  return 1;
}

sub construct_geometries_from_normal_modes {
  my ($filenamesRef, $q_valuesRef, $save_geoms_to_disk, $r_geometry) = @_;
     my @filenames = @{$filenamesRef};
     my @q_values = @{$q_valuesRef};

     $n_modes = scalar @filenames;
     $n_modes_chck = scalar @q_values;

     #Error checking and info printing
     if ($n_modes != $n_modes_chck) {
        die "Error in construct_geometries_from_normal_modes: on input the number of normal modes is not consistent!\n"
     }
     print "Number of normal modes: $n_modes\n";

     print "Starting geometry:\n";
     $eq_geom = $r_geometry->{'geometries'}->[0];

     $n_at = scalar @{$eq_geom->{'atoms'}};
     for (my $i = 0; $i < scalar @{$eq_geom->{'atoms'}}; $i++) {
       my $r_atom = $eq_geom->{'atoms'}->[$i];
       print " $r_atom->[0]";
       print " $r_atom->[1]";
       print " $r_atom->[2]";
       print " $r_atom->[3]\n";
     }

     # Read-in the normal mode displacements into the 'geometries' hash array.
     $i = -1;
     @n_qs = ();
     @modes = ();
     foreach(@q_values) {
        $n = @$_;
        push @n_qs, $n;
        $i = $i + 1;
        print "For mode $filenames[$i] we're using $n displacements\n";

        if (open(NORMALMODE, "< $filenames[$i]")) {
           $mode="";
           $nm = -1;
           push @modes, %eq_geom; #copy the contents of the equilibrium geometry
           while (my $line = <NORMALMODE>) {
              if ($line =~ m/^\s*([+-]?[0-9]+\.[0-9]+)\s+([+-]?[0-9]+\.[0-9]+)\s+([+-]?[0-9]+\.[0-9]+)/i) {
                 $nm = $nm + 1;
                 $modes[$i]->{'atoms'}->[$nm]->[1] = $1;
                 $modes[$i]->{'atoms'}->[$nm]->[2] = $2;
                 $modes[$i]->{'atoms'}->[$nm]->[3] = $3;
                 my $r_atom = $modes[$i]->{'atoms'}->[$nm];
                 $mode = $mode."$r_atom->[1] $r_atom->[2] $r_atom->[3]\n";
              }
           }
           print "NORMAL MODE (q):\n";
           print "$mode";
           if ($nm+1 != $n_at) {die "$nm; $n_at Normal mode incompatible with starting geometry! \n";}
           close NORMALMODE;
        } else {
           die "Error in construct_geometries_from_normal_modes: couldn't open the file $filenames[$i] !\n";
        }
     }

     # Erase all geometries input by the user: they will be replaced by the displaced geometries.
     $r_geometry->{'geometries'} = [];

     if ($n_modes == 1) {

        $r_geometry->{'geometry_labels'} = " q1 ";

        #Implement the single mode by a fake 2-mode set-up
        $n_qs[1] = 1; #only one displacement of the 2nd mode
        $q_values[1]->[0] = 0.0; #zero value of the q-displacement of the 2nd mode
        $mode2_ind = 0; #2nd mode index points in fact to the first mode

     } elsif ($n_modes == 2) {

       $r_geometry->{'geometry_labels'} = " q1    q2 ";
       $mode2_ind = 1; #point to the second mode

     } else {
       die "construct_geometries_from_normal_modes: Not yet implemented for more than two modes!";
     }

     $ngeom = 0;
     for($q1 = 0; $q1 < $n_qs[0]; $q1++) {

        $q1_displ = $q_values[0]->[$q1];

        for($q2 = 0; $q2 < $n_qs[1]; $q2++) {

           $q2_displ = $q_values[1]->[$q2];
           $ngeom++;

           print "DISPLACED GEOMETRY $ngeom: q1 = $q1_displ, q2 = $q2_displ\n";

           push(@{$r_geometry->{'geometries'}},
            { 'description', " $ngeom  $q1_displ  $q2_displ",
              'gnuplot_desc', "Geometry $ngeom: q1 = $q1_displ, q2 = $q2_displ", # used in gnuplot files for keys
              'atoms', [ ]
            }
           );

           if ($save_geoms_to_disk == 1) {open(DISPLACED, "> geom_$ngeom") or die $!;}

           for (my $i = 0; $i < $n_at; $i++) {
              my $q1_atom = $modes[0]->{'atoms'}->[$i];
              my $q2_atom = $modes[$mode2_ind]->{'atoms'}->[$i];
              my $eq_atom = $eq_geom->{'atoms'}->[$i];

              my $x = $eq_atom->[1] + $q1_displ*$q1_atom->[1] + $q2_displ*$q2_atom->[1];
              my $y = $eq_atom->[2] + $q1_displ*$q1_atom->[2] + $q2_displ*$q2_atom->[2];
              my $z = $eq_atom->[3] + $q1_displ*$q1_atom->[3] + $q2_displ*$q2_atom->[3];

              if ($save_geoms_to_disk == 1) {print DISPLACED "$eq_geom->{'atoms'}->[$i]->[0] $x $y $z\n";}
              print "$eq_geom->{'atoms'}->[$i]->[0] $x $y $z\n";

              push(@{$r_geometry->{'geometries'}->[$ngeom-1]->{'atoms'}},
                [ $eq_geom->{'atoms'}->[$i]->[0], $x, $y, $z ]
              );
           }

           if ($save_geoms_to_disk == 1) {close DISPLACED;}

        }
     }

     return 1;
}

sub generate_geometries {
  my ($r_par, $r_geometry) = @_;

  my @atoms = @{$r_par->{'model'}->{'atoms'}};
  my $na = scalar @atoms;                # number of atoms
  my $dir = $r_par->{'dirs'}->{'model'}; # model directory

  my $r_geom;           # auxiliary reference to one geometry
  my $dir_geom;         # geometry directory
  my @geometries = ();  # array of all geometries - returned at the end
  my $ng = 0;           # number of geometries

  open(GEOM, ">$dir${bs}geometries"); # File with geometries
  open(MOLDEN, ">$dir${bs}geometries.molden"); # File with geometries in molden format

  # For consistency reasons here we set the same length unit as given in %geometry
  $r_par->{'model'}->{'r_unit'} = $r_geometry->{'length_unit'};
  $r_par->{'data'}->{'geom_labels'} = "  $r_geometry->{'geometry_labels'}";
  print GEOM '#geom', $r_par->{'data'}->{'geom_labels'}, "\n";
  foreach my $r_given_geom ( @{$r_geometry->{'geometries'}} ) {

    # directory
    $ng++;
    $dir_geom = "$dir${bs}geom$ng";
    &make_dir($dir_geom);

    # general settings for geometry
    $r_geom = {
      'geometry', "  $r_given_geom->{'description'}",  # string to use in output files
      'gnuplot_desc', $r_given_geom->{'gnuplot_desc'}, # gnuplot description                   \
      'dir', $dir_geom,                                # directory                              > of a given geometry
      'symmetry', $r_par->{'model'}->{'symmetry'},     # symmetry (this is not used right now) /
      'atoms', $r_given_geom->{'atoms'}
    };
    if ($r_geometry->{'correct_cm'} == 1) { &make_cm_correction($r_geom); }
    &coord_to_strings($r_geom);
    push @geometries, $r_geom;

    # Write information about geometry into the string (used also later in files for potential curves etc.)
    # and in the file 'geometries' and 'geometries.molden'
    print GEOM sprintf("%5d", $ng), $r_geom->{'geometry'}, "\n";
    print MOLDEN "! Geometry ".sprintf("%5d", $ng)."\n";
    print MOLDEN "[Atoms] ".sprintf("%s", $r_geometry->{'length_unit'} == 0 ? "AU" : "Angs")."\n";
    for (my $i = 0; $i < scalar @{$r_geom->{'atoms'}}; $i++) {
      my $r_atom = $r_geom->{'atoms'}->[$i];
      print MOLDEN sprintf("%-4s", $r_atom->[0]);
      print MOLDEN sprintf("%3d", $i + 1);
      print MOLDEN sprintf("%5d", $atomic_number{$r_atom->[0]});
      print MOLDEN sprintf("%20.10f", $r_atom->[1]);
      print MOLDEN sprintf("%20.10f", $r_atom->[2]);
      print MOLDEN sprintf("%20.10f", $r_atom->[3])."\n";
    }
  }
  &print_info("Number of geometries to run: $ng\n", $r_par);

  close(GEOM);
  return @geometries;
}

sub deg2rad {
  my ($angle) = @_;
  return $angle * 3.14159265358979323846 / 180;
}

sub make_cm_correction {
  my ($r_geom) = @_;
  my $n = scalar @{$r_geom->{'atoms'}};

  # determine cm
  my $mT = 0.0;             # total mass
  my @xT = (0.0, 0.0, 0.0); # coordinates of CM
  for (my $i = 0; $i < $n; $i++) { # read all atoms
     my $ma = $mass{$r_geom->{'atoms'}->[$i][0]};
     $mT += $ma;
     for (my $ir = 1; $ir <= 3; $ir++) { $xT[$ir - 1] += $ma * $r_geom->{'atoms'}->[$i][$ir]; }
  }
  foreach (@xT) { $_ = $_/$mT; }

  # shift all atoms
  for (my $i = 0; $i < $n; $i++) {
     for (my $ir = 1; $ir <= 3; $ir++) { $r_geom->{'atoms'}->[$i][$ir] -= $xT[$ir - 1]; }
  }
  return 1;
}

sub coord_to_strings {
  my ($r_geom) = @_;
  my $n = scalar @{$r_geom->{'atoms'}};
  for (my $i = 0; $i < $n; $i++) {
     foreach (@{$r_geom->{'atoms'}->[$i]}[1..3]) { $_ = sprintf("%13.9f", $_); }
  }
}

# ============ end of subroutines for GEOMETRIES ============

# ================= subroutines for INPUTS ==================

# according to model settings this subroutine searches for the continuum basis file
# if a specific file with the given L is not available, a file with the basis for higher L is used
sub which_continuum_basis_file_to_use {
  my ($r_par) = @_;
  my $basis = "q$r_par->{'data'}->{'target'}->{'charge'}.r".sprintf("%.0f", $r_par->{'model'}->{'radius_GTO'}).".L";
  if (-e "$r_par->{'dirs'}->{'basis'}${bs}swmol3.continuum.$basis$r_par->{'model'}->{'maxl_GTO'}") {
    $basis .= $r_par->{'model'}->{'maxl_GTO'};
  } else {
    foreach my $file (grep(/continuum/, &dirfilelist($r_par->{'dirs'}->{'basis'}))) {
      if ($file =~ m/$basis(\d+)/) {
        if ($1 > $r_par->{'model'}->{'maxl_GTO'}) { $basis .= $1; }
      }
    }
  }
  if ($basis =~ m/L$/) {
    &print_info("No suitable continuum basis for $basis$r_par->{'model'}->{'maxl_GTO'} exists !\n", $r_par);
    die;
  }
  $r_par->{'data'}->{'scattering'}->{'basis'} = $basis;
  return 1;
}

sub get_basis_for_swmol3_input {
  my ($dir, $basis, $atom, $xxx, $yyy, $zzz, $r_par) = @_;
  foreach ($xxx, $yyy, $zzz) { $_ = sprintf("%13.9f", $_); }
  my $str = "";
  if (!(&read_file("$dir${bs}swmol3.$atom.$basis", \$str))) {
    my %molpro_basis = ();
    &read_molpro_basis_to_hash($dir, $basis, $atom, \%molpro_basis);
    $str = molpro_basis_to_swmol3_string(\%molpro_basis);
  }
  &replace_in_template(\$str, "XXX", $xxx);
  &replace_in_template(\$str, "YYY", $yyy);
  &replace_in_template(\$str, "ZZZ", $zzz);
  if ($atom eq "continuum") {
    my $L = $r_par->{'model'}->{'maxl_GTO'} + 1;
    $str =~ s/iqm\s+=\s+\d+/iqm = $L/sg;
  }
  return $str;
}

# Next subroutine creates a hash array of the following structrure
# containing all information about basis from Molpro format
# e.g. for oxygen in the cc-pVTZ basis we should get
# %data = (
#   'atom',         'O',
#   'comments',     '! OXYGEN       (10s,5p,2d,1f) -> [4s,3p,2d,1f]\n',
#   'exponents',    { 's', [15330.0000000, 2299.0000000, ...],
#                     'p', [34.4600000, 7.7490000, ...],
#                      ... }
#   'contractions', { 's', [ [ 1, 8, 0.0005080, 0.0039290, ...],
#                            [ 1, 8, -0.0001150, -0.0008950, ...],
#                            [ 9, 9, 1]
#                            ...],
#                     'p', ... }
# );
sub read_molpro_basis_to_hash {
  my ($dir, $basis, $atom, $r_data) = @_;
  %$r_data = ();
  $r_data->{'exponents'} = {};
  $r_data->{'contractions'} = {};
  my $basis_comment = "";
  my $atom_ok = 0;
  my $current_orb = "";
  if (!(open(BASIS, "$dir${bs}molpro.$atom.$basis"))) {
    if (!(open(BASIS, "$dir${bs}molpro.$basis"))) {
      if (!(open(BASIS, "$dir${bs}$basis"))) {
        die "A file with the basis $basis does not exist !!!\n";
      }
    }
  }
  LINE: while (my $line=<BASIS>) {
    if ($line =~ /^!/) {
      last LINE if $atom_ok == 1;
      $basis_comment .= $line;
    }
    else {
      chomp($line);
      $line =~ s/\s*$//;
      if ($line =~ /^\s*([$atomic_orbitals])\s*,\s*(\w+)\s*,\s*(.*)$/i) {
        $current_orb  = $1;
        my $current_atom = $2;
        my $current_exp  = $3;
        # check whether we have the right atom
        if (lc($current_atom) eq lc($atom)) {
          $atom_ok = 1;
          # write the atom and comment
          if (lc($current_orb) eq "s") {
            $r_data->{'atom'} = $atom;
            $r_data->{'comment'} = $basis_comment;
          }
          # write exponents and create an empty array of contractions
          $r_data->{'exponents'}->{$current_orb} = [split(/\s*,\s*/, $current_exp)];
          $r_data->{'contractions'}->{$current_orb} = [];
        } # if (lc($current_atom) eq lc($atom))
        $basis_comment = "";
      }
      elsif ($line =~ /^\s*c\s*,\s*(\d+)\.(\d+)\s*,\s*(.*)$/i) {
        if ($atom_ok == 1) {
          # save another contraction
          push(@{$r_data->{'contractions'}->{$current_orb}}, [$1, $2, split(/\s*,\s*/, $3)]);
        }
      }
    } # if ($line =~ /^!/)
  } # while (my $line=<BASIS>)
  close(BASIS);
  if ($atom_ok == 0) { die "A file with the basis $basis for $atom does not exist !!!\n";}
  else {
    # Ulozit posledni kontrakce
  }
  return 1;
}

# Subroutine which writes Molpro basis (read by read_molpro_basis_to_hash
# into a hash array $r_data) to a string
sub molpro_basis_to_molpro_string {
  my ($r_data) = @_;
  my $str = "";
  $str .= $r_data->{'comment'};
  ORBS: foreach my $orb (split('', $atomic_orbitals)) {
    if ($r_data->{'exponents'}->{$orb}) {
      $str .= "$orb, $r_data->{'atom'}, ".join(', ', @{$r_data->{'exponents'}->{$orb}})."\n";
      foreach my $r_contr (@{$r_data->{'contractions'}->{$orb}}) {
        $str .= "c, $r_contr->[0].$r_contr->[1], ".join(', ', @{$r_contr}[2..(scalar @$r_contr - 1)])."\n";
      }
    }
    else { last ORBS; }
  }
  return $str;
}

# Subroutine which converts Molpro basis format (read by read_molpro_basis_to_hash
# into a hash array $r_data) to swmol3 basis format and returns a string containing
# &atom namelist, geometry is given as
#   atcord = >>>XXX<<<, >>>YYY<<<, >>>ZZZ<<<
# and must be replaced elsewhere
sub molpro_basis_to_swmol3_string {
  my ($r_data) = @_;
  my $iqm = scalar keys %{$r_data->{'exponents'}}; # number of atomic symmetries
  my $str = "";
  $str .= '&atom'."\n";
  $str .= '  atcord = >>>XXX<<<, >>>YYY<<<, >>>ZZZ<<<,'."\n";
  $str .= "  atnam = '$r_data->{'atom'}',\n";
  $str .= "  q = $atomic_number{$r_data->{'atom'}},\n";
  $str .= "  iqm = $iqm,\n";
  $str .= "  jco = ".join(',', (1) x $iqm).",\n";
  $str .= $r_data->{'comment'};
  $str .= "/\n";
  ORBS: foreach my $orb (split('', $atomic_orbitals)) {
    if ($r_data->{'exponents'}->{$orb}) {
      my $n_exp = scalar @{$r_data->{'exponents'}->{$orb}};
      my $r_contr = $r_data->{'contractions'}->{$orb};
      my $n_contr = scalar @$r_contr;
      # line with number of exponents and number of contractions
      $str .= sprintf("%5d%5d%5d", $n_exp, $n_contr, 0)."\n";
      # for each exponent write contration coeeficients
      for (my $i = 0; $i < $n_exp; $i++) {
        $str .= sprintf("%15.7f", $r_data->{'exponents'}->{$orb}->[$i]);
        for (my $j = 0; $j < $n_contr; $j++) {
          # check whether the coefficient for the current exponent is non-zero
          if ($i + 1 >= $r_contr->[$j]->[0] && $i + 1 <= $r_contr->[$j]->[1]) {
            $str .= sprintf(" %15.7E", $r_contr->[$j]->[$i + 3 - $r_contr->[$j]->[0]]);
          }
          else {
            $str .= sprintf(" %10.7f     ", 0.0);
          }
        }
        $str .= "\n";
      }
    }
    else { last ORBS; }
  }
  $str .= "/\n";
  return $str;
}


# ------------------------- MAS setup ------------------------

sub setup_mas {
  my ($mas_type, $r_par) = @_;

  my $prefix = "";
  my $suffix = "";

  if ($mas_type eq 'qchem' && $r_par->{'model'}{'qchem_MAS'}){  # MAS for quantum chemistry calculation
    $prefix = "qchem_";
    $suffix = "SCF";
  }

  my $MAS = MultiSpace::set_orbs_per_irrep_per_subspace_automatically(
    $r_par->{'model'}{$prefix.'MAS'}, $r_par->{'data'}{'MAS'}{'sorted_orb'}
  );

  my $target_space = MultiSpace->new({
    symmetry => $r_par->{'model'}->{'symmetry'},
    MAS => $MAS,
    constraints => $r_par->{'model'}->{$prefix.'constraints'},
    type => $r_par->{'model'}->{'use_MAS'.$suffix}}
  );

  $r_par->{'data'}->{'MAS'}->{$mas_type eq 'qchem' ? 'qchem' : 'target'} = $target_space;

  if ($mas_type eq 'rmat' && $r_par->{'model'}->{'l2_MAS'}){
      my $l2_MAS = MultiSpace::set_orbs_per_irrep_per_subspace_automatically(
        $r_par->{'model'}{'l2_MAS'}, $r_par->{'data'}{'MAS'}{'sorted_orb'}
      );
      my $l2_space = MultiSpace->new({
        symmetry => $r_par->{'model'}->{'symmetry'},
        MAS => $l2_MAS,
        constraints => $r_par->{'model'}->{'l2_constraints'},
        type => $r_par->{'model'}->{'use_MAS'}}
      );
      $r_par->{'data'}->{'MAS'}->{'l2'} = $l2_space;
  }
  return 1;
}

sub choose_model_by_name_mas {
  my ($r_par) = @_;

  my $model_type = $r_par->{'model'}->{'model'};

  if (!exists($r_par->{'model'}->{'models_using_mas'}->{$model_type})) {
    die "$model_type is not a recognized model.";
  } elsif (!$r_par->{'model'}->{'models_using_mas'}->{$model_type}){
    die "$model_type model not yet implemented in the MAS approach";
  }

  my $nel_closed = min($r_par->{'model'}->{'nfrozen'}*2, $r_par->{'model'}->{'nelectrons'});
  my $nel_active = $r_par->{'model'}->{'nelectrons'} - $nel_closed; # ORMAS style entry

  my @frozen_orbs = (0) x 8;
  my @active_orbs = (0) x 8;
  my @virtual_orbs = (0) x 8;
  my @pco_orbs = (0) x 8;
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    if ($r_par->{'model'}{'nfrozen'})   {$frozen_orbs[$i] += $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i];}
    if ($r_par->{'model'}{'nactive'})   {$active_orbs[$i] += $r_par->{'data'}->{'orbitals'}->{'active'}->[$i];}
    if ($r_par->{'model'}{'nvirtual'})  {$virtual_orbs[$i] += $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i]};
    if ($r_par->{'model'}->{'use_PCO'}) {$pco_orbs[$i] += $r_par->{'data'}->{'orbitals'}->{'PCO'}->[$i]};
  }

  my @MAS;

  # Add closed orbitals if any.
  if ($r_par->{'model'}{'nfrozen'} > 0){
    if ($r_par->{'model'}->{'use_PCO'} && $r_par->{'model'}->{'reduce_PCO_CAS'}) {
      # We allow single excitations from the closed space to PCOs in the 'reduce_PCO_CAS' model
      push(@MAS, (\@frozen_orbs, [$nel_closed-1, $nel_closed], "closed") );
    } else {
      push(@MAS, (\@frozen_orbs, [$nel_closed, $nel_closed], "closed") );
    }
  }

  # Add active orbitals if any.
  my @CAS;
  if ($r_par->{'model'}{'nactive'} > 0){
    if ($r_par->{'model'}->{'use_PCO'} && $r_par->{'model'}->{'reduce_PCO_CAS'}) {
      # We need to split the CAS into 2 subspaces: active part of GS and
      # the rest of the CAS active orbitals. All occupancies are allowed so it
      # is equivalent to the single space description of the CAS.
      @CAS = split_cas($r_par);
    } else {
      push(@CAS,  (\@active_orbs, [$nel_active, $nel_active], "active"));
    }
    push(@MAS, @CAS);
  }

  # Add PCO if we have them.
  if ($r_par->{'model'}->{'use_PCO'}) {
    if ($r_par->{'model'}{'nvirtual'} > 0) {
      # Virtuals must come before PCO so must be added even though they are not
      # used in the target run.
      if ($model_type eq "CAS-A") {
        die "Error: 'nvirtual' must be zero for CAS-A with PCO";
      }
      push(@MAS,  (\@virtual_orbs, [0, 0], "active"));
    }
    add_excitation(\@MAS);
    push(@MAS,  (\@pco_orbs, [0, 1], "active"));
  }
  $r_par->{'model'}->{'MAS'} = \@MAS;

  # We need to set further constraints for the reduce_PCO_CAS model
  if ($r_par->{'model'}->{'use_PCO'} && $r_par->{'model'}->{'reduce_PCO_CAS'}) {
    $r_par->{'model'}->{'constraints'} = sub {
      my $occ = shift;
      my $val = 1;
      # Only allow PCO for single excitations from GS
      if ($r_par->{'model'}{'nfrozen'} > 0){
        $val = !($occ->[-1] > 0 && $occ->[0] + $occ->[1] != $r_par->{'model'}->{'nelectrons'} - 1);
        $val = $val && !($occ->[-1] == 0 && $occ->[0] != $nel_closed);
      } else {
        $val = !($occ->[-1] > 0 && $occ->[0] != $r_par->{'model'}->{'nelectrons'}-1);
      }
      return $val;
    }
  }

  my @l2_MAS = @{dclone(\@MAS)};
  add_electron(\@l2_MAS);
  $polarization = 0;

  if ($r_par->{'model'}->{'nvirtual'} > 0) {
    if ($model_type eq "CAS" || $model_type eq "CAS-B") { # Spingle excitations to the virtuals
      $polarization = 1;
      if ($r_par->{'model'}->{'use_PCO'}) {
        $l2_MAS[-5]->[1] = 1;
      } else {
        push(@l2_MAS, (\@virtual_orbs, [0, 1], "active"));
      }
    } elsif ($model_type eq "SEP" || $model_type eq "CAS-C"){ # Spingle and double excitations to the virtuals
      $polarization = 2;
      add_excitation(\@l2_MAS);
      if ($r_par->{'model'}->{'use_PCO'}) {
        $l2_MAS[-5]->[1] = 2 ;
      } else {
        push(@l2_MAS, (\@virtual_orbs, [0, 2], "active"));
      }
    }
  }
  $r_par->{'model'}->{'l2_MAS'} = \@l2_MAS;

  # We need to set further constraints if we have PCO
  if ($r_par->{'model'}->{'use_PCO'}) {
    $r_par->{'model'}->{'l2_constraints'} = sub {
      my $occ = shift;
      my $nelectrons = $r_par->{'model'}->{'nelectrons'};
      my $val = 1;

      if ($polarization) {
        # Virtuals and PCO should not be simultaneously populated.
        $val = !($occ->[-1] > 0 && $occ->[-2] > 0);
      }
      if ($r_par->{'model'}->{'reduce_PCO_CAS'}){
        # Only PCO CSF of the type  GS^N PCO^1 and GS^N-1 PCO^2
        if ($r_par->{'model'}{'nfrozen'} > 0){
          if ( ($occ->[-1] == 0 && $occ->[0] == $nel_closed)
                ||  ($occ->[-1] == 2 && $occ->[0] + $occ->[1] == $nelectrons-1)
                ||  ($occ->[-1] == 1 && $occ->[0] + $occ->[1] == $nelectrons) ) {
            $val = $val && 1;
          } else {
            $val = 0;
          }
        } else {
          if ($occ->[-1] == 0 ||  ($occ->[-1] == 2 && $occ->[0] == $nelectrons-1)
              || ($occ->[-1] == 1 && $occ->[0] == $nelectrons)) {
            $val = $val && 1;
          } else {
            $val = 0;
          }
        }
      }
      return $val;
    }
  }
  return 1;
}

sub split_cas {
  my ($r_par) = @_;
  # For use with ORMAS only.

  my @frozen_orbs = (0) x 8;
  my @active_orbs = (0) x 8;
  my @HartreeFockGS_orbs = (0) x 8;

  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    if ($r_par->{'model'}{'nfrozen'})  {$frozen_orbs[$i] += $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i];}
    $active_orbs[$i] += $r_par->{'data'}->{'orbitals'}->{'active'}->[$i];
    $HartreeFockGS_orbs[$i] += $r_par->{'data'}->{'orbitals'}->{'GS'}->[$i];
  }

  my @CAS;
  my @active_GS_orbs = (0) x 8;  # GS orbitals that are not closed.
  my @active_virtual_orbs = (0) x 8;  # Virtual to GS but still in CAS.

  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      $active_GS_orbs[$i] = $HartreeFockGS_orbs[$i] - $frozen_orbs[$i];
      $active_virtual_orbs[$i] = $active_orbs[$i] - $active_GS_orbs[$i];
  }

  # No restrictions on the number electrons in the two spaces - the alg. will ensure
  # that the correct number of electrons are placed
  push(@CAS, (\@active_GS_orbs, [0, 2*sum0(@active_GS_orbs)], "active"));
  push(@CAS, (\@active_virtual_orbs, [0, 2*sum0(@active_virtual_orbs)], "active"));

  return @CAS;
}

sub add_electron{
  my ($MAS) = @_;

  # ORMAS with subspace labels only
  for( my $i = 0; $i < @{$MAS}; $i += 3 ) {
    if ($MAS->[$i+1][1] < 2*sum0(@{$MAS->[$i]})) {
      $MAS->[$i+1][1]++;
    }
  }
}

sub add_excitation{
  my ($MAS) = @_;

  # ORMAS with subspace labels only
  for( my $i = 0; $i < @{$MAS}; $i += 3 ) {
    if ($MAS->[$i+1][0] > 0 && $MAS->[$i+2] eq "active") {
      $MAS->[$i+1][0]--;
    }
  }
}

# ------------------------- molpro --------------------------

sub make_molpro_input {
  my ($r_par, $r_str) = @_;
  my $r_geom = $r_par->{'data'}->{'geom'};
  my $natoms = scalar @{$r_geom->{'atoms'}};
  my $sym_op = "SYMMETRY,".$symmetries{$r_geom->{'symmetry'}}->[1]; # Symmetry operations
  if ($r_geom->{'symmetry'} eq "C1") { $sym_op = "NOSYM"; }
  $sym_op =~ s/\'//g; # In MOLPRO input there are no ''

  # Substitutions to input
  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "BASISNAME",$r_par->{'model'}->{'basis'});
  if ($r_par->{'model'}->{'r_unit'} == 1) { &replace_in_template($r_str, "ANGSTROM", "ANGSTROM"); }
  else                                       { &replace_in_template($r_str, "ANGSTROM", ""); }
  &replace_in_template($r_str, "SYMOP", $sym_op);

  # Adding geometries and basis sets for each atom
  my $str_basis = "";
  my $str_geometry = "";
  my %used_basis = ();
  for (my $i = 0; $i < $natoms; $i++) {
    my ($atom, $xxx, $yyy, $zzz) = @{$r_geom->{'atoms'}->[$i]};
    if (!$used_basis{$atom}) {
      my %molpro_basis = ();
      &read_molpro_basis_to_hash($r_par->{'dirs'}->{'basis'}, $r_par->{'model'}->{'basis'}, $atom, \%molpro_basis);
      $str_basis .= molpro_basis_to_molpro_string(\%molpro_basis);
      $used_basis{$atom} = 1;
    }
    foreach ($xxx, $yyy, $zzz) { $_ = sprintf("%13.9f", $_); }
    $str_geometry .= "$atom,, $xxx, $yyy, $zzz\n";
  }
  &replace_in_template($r_str, "BASIS", $str_basis);
  &replace_in_template($r_str, "GEOMETRY", $str_geometry);

  if ($r_par->{'model'}->{'model'} =~ /^SE/ || $r_par->{'data'}->{'scf_ok'} == 0 ||
        $r_par->{'model'}->{'orbitals'} eq "HF") {
    &replace_in_template($r_str, "METHOD", "HF");
  }
  else {
    &replace_in_template($r_str, "METHOD", "CASSCF");
  }

  # Hartree-Fock wave function (must be set if target has charge)
  my $nelec = $r_par->{'model'}->{'nelectrons'};
  if ($r_par->{'model'}->{'charge_of'} eq "scattering") { $nelec++; }
  &replace_in_template($r_str, "HFWF", "wf,nelec=$nelec;");

  # Adding closed and occupied orbitals and wf
  # if MOLPRO has not run yet and CASSCF has to be set
  # we have to use information provided by user
  if (!(@{$r_par->{'data'}->{'orbitals'}->{'target'}})) {
    # Target orbitals = frozen + active
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i] = $r_par->{'model'}->{'frozen_orbs'}->[$i];
      $r_par->{'data'}->{'orbitals'}->{'target'}->[$i] = $r_par->{'model'}->{'frozen_orbs'}->[$i]
                                                       + $r_par->{'model'}->{'active_orbs'}->[$i];
    }
  }

  my $target_space;
  my %molpro_input;
  if ($r_par->{'data'}->{'scf_ok'} == 1 && exists($r_par->{'model'}->{'use_MASSCF'}) && $r_par->{'model'}->{'use_MASSCF'} == 1){
    $target_space = $r_par->{'data'}{'MAS'}{'qchem'};
    if ($target_space->{'type'} != 1){
      die "ORMAS approach must be used for Molpro ORMASSCF run!";
    }
    %molpro_input = $target_space->molpro_ormas_input();
    # TODO: AH: Perhaps I should reset $r_par->{'data'}->{'orbitals'}->{'frozen'} etc here?
    &replace_in_template($r_str, "CLOSED", join(",", @{$molpro_input{'CLOSED'}}));
    &replace_in_template($r_str, "OCC", join(",", @{$molpro_input{'OCC'}}));
  } else {
    &replace_in_template($r_str, "CLOSED", join(",", @{$r_par->{'data'}->{'orbitals'}->{'frozen'}}));
    &replace_in_template($r_str, "OCC", join(",", @{$r_par->{'data'}->{'orbitals'}->{'target'}}));
  }

  # Adding all 'wf' to MOLPRO input for which we want to optimise orbitals
  # and also for all target states needed
  my $casscf_states = "";
  my $all_states = "";
  if (exists($r_par->{'model'}->{'use_MASSCF'}) && $r_par->{'model'}->{'use_MASSCF'} == 1){
    if ($molpro_input{'ORMAS'}){
      $casscf_states .= $molpro_input{'ORMAS'};
      my $commented_lines = "";
      foreach my $line (split(/\n/, $molpro_input{'ORMAS'})){
        if ($line ne "config;"){
          $line = "!".$line;
        }
        $commented_lines .= $line."\n";
      }
      $all_states .= $commented_lines."!";
    }
  }
  for (my $i = 1; $i <= $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'ncasscf_states'}}) {
      my $spin = $spin_multiplicity{$statespin} - 1;
      my $nstates = $r_par->{'model'}->{'ncasscf_states'}->{$statespin}->[$i - 1];
      if ($nstates > 0) {
        $casscf_states .= "wf,nelec=$nelec,sym=$i,spin=$spin; state,$nstates; ";
      }
      my $nstates = $r_par->{'model'}->{'ntarget_states'}->{$statespin}->[$i - 1];
      if ($nstates > 0) {
        $all_states .= "wf,nelec=$nelec,sym=$i,spin=$spin; state,$nstates; ";
      }
    }
  }
  &replace_in_template($r_str, "CASWF", $casscf_states);
  &replace_in_template($r_str, "ALLWF", $all_states);
  return 1;
}

# ------------------------- psi4 --------------------------

sub make_psi4_input {
  my ($r_par, $r_str) = @_;
  my $r_geom = $r_par->{'data'}->{'geom'};
  my $nir = $r_par->{'data'}->{'nir'};
  my $natoms = scalar @{$r_geom->{'atoms'}};
  my $nelectrons = $r_par->{'model'}->{'nelectrons'};
  my $nprotons = 0;
  my @frozen, @active;
  my $groundStateSpinMultiplicity = $r_par->{'model'}->{'groundStateSpinMultiplicity'};

  # Compose the geometry section and calculate number of +1 charges
  $str_geometry = "";
  for (my $i = 0; $i < $natoms; $i++) {
    my ($atom, $xxx, $yyy, $zzz) = @{$r_geom->{'atoms'}->[$i]};
    foreach ($xxx, $yyy, $zzz) { $_ = sprintf("%13.9f", $_); }
    $str_geometry .= "    $atom $xxx $yyy $zzz\n";
    $nprotons += $atomic_number{$atom};
  }

  # Determine which electron count (target/scattering) to use for orbitals
  if ($r_par->{'model'}->{'charge_of'} eq "scattering") { $nelectrons++; }

  # Set the QC method
  if ($r_par->{'model'}->{'model'} =~ /^SE/ || $r_par->{'data'}->{'scf_ok'} == 0 || $r_par->{'model'}->{'orbitals'} eq "HF") {
    &replace_in_template($r_str, "METHOD", "scf");
  }
  else {
    &replace_in_template($r_str, "METHOD", "casscf");
  }

  # Target orbitals = frozen + active
  if (!(@{$r_par->{'data'}->{'orbitals'}->{'target'}})) {
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i] = $r_par->{'model'}->{'frozen_orbs'}->[$i];
      $r_par->{'data'}->{'orbitals'}->{'target'}->[$i] = $r_par->{'model'}->{'frozen_orbs'}->[$i]
                                                       + $r_par->{'model'}->{'active_orbs'}->[$i];
    }
  }

  # Reorder the irreducible representations from the Molpro order to the Psi4 order
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    $frozen->[$i] = $r_par->{'model'}->{'frozen_orbs'}->[$irred_psi4{$r_par->{'model'}->{'symmetry'}}->[$i]];
    $active->[$i] = $r_par->{'model'}->{'active_orbs'}->[$irred_psi4{$r_par->{'model'}->{'symmetry'}}->[$i]];
  }

  # Average energies of all required targets, but only with single IRR; as of version 1.3.2, Psi4 still
  # can't average energies of states with different IRRs.
  my $numroots = $r_par->{'model'}->{'ntarget_states_used'};
  if ($numroots < 1 || $nir > 1) { $numroots = 1 }

  &replace_in_template($r_str, "NUMROOTS", $numroots);
  &replace_in_template($r_str, "AVGSTATES", join(",", map { $_ } 0 .. $numroots-1));
  &replace_in_template($r_str, "AVGWEIGHTS", join(",", map { 1 } 0 .. $numroots-1));

  &replace_in_template($r_str, "CLOSED", join(",", @$frozen));
  &replace_in_template($r_str, "ACTIVE", join(",", @$active));

  &replace_in_template($r_str, "CHARGE",    $nprotons - $nelectrons);
  &replace_in_template($r_str, "SPIN",      $groundStateSpinMultiplicity);
  &replace_in_template($r_str, "REFERENCE", $groundStateSpinMultiplicity  == 1 ? "rhf" : "rohf");

  &replace_in_template($r_str, "MOLECULE",  lc($r_par->{'model'}->{'molecule'}));
  &replace_in_template($r_str, "BASISNAME", $r_par->{'model'}->{'basis'});
  &replace_in_template($r_str, "GEOMETRY",  $str_geometry);
  &replace_in_template($r_str, "SYMMETRY",  $r_geom->{'symmetry'});
  &replace_in_template($r_str, "UNIT",      $r_par->{'model'}->{'r_unit'} == 1 ? "ang" : "au");
}

# ------------------------- molcas --------------------------

sub make_molcas_input {
  my ($r_par, $r_str) = @_;
  my $r_geom = $r_par->{'data'}->{'geom'};
  my $nir = $r_par->{'data'}->{'nir'};
  my $natoms = scalar @{$r_geom->{'atoms'}};
  my $nelectrons = $r_par->{'model'}->{'nelectrons'};
  my $sym_op = $symmetries{$r_geom->{'symmetry'}}->[1]; # Symmetry operations
  my $point_group = $r_par->{'model'}->{'symmetry'};

  $sym_op =~ s/'/$1/g;
  $sym_op =~ s/,/$1/g;

  my $nprotons = 0;

  # Compose the geometry section and calculate number of +1 charges
  $str_geometry = "";
  for (my $i = 0; $i < $natoms; $i++) {
    my ($atom, $xxx, $yyy, $zzz) = @{$r_geom->{'atoms'}->[$i]};
    foreach ($xxx, $yyy, $zzz) { $_ = sprintf("%13.9f", $_); }
    $str_geometry .= "    $atom $xxx $yyy $zzz\n";
    $nprotons += $atomic_number{$atom};
  }

  # Determine which electron count (target/scattering) to use for orbitals
  if ($r_par->{'model'}->{'charge_of'} eq "scattering") { $nelectrons++; }

  &replace_in_template($r_str, "NATOMS",  $natoms);
  &replace_in_template($r_str, "UNIT",      $r_par->{'model'}->{'r_unit'} == 1 ? "Angstrom" : "a.u.");
  &replace_in_template($r_str, "GEOMETRY",  $str_geometry);
  &replace_in_template($r_str, "SYMOP", $sym_op);
  &replace_in_template($r_str, "BASISNAME", $r_par->{'model'}->{'basis'});
  &replace_in_template($r_str, "CHARGE",    $nprotons - $nelectrons);

  if ($r_par->{'model'}->{'model'} =~ /^SE/ || $r_par->{'data'}->{'scf_ok'} == 0 ||
        $r_par->{'model'}->{'orbitals'} eq "HF") {
    &replace_in_template($r_str, "METHOD", "HF");
  }
  else {
    &replace_in_template($r_str, "METHOD", "GASSCF");
  }
  #  Choose spin, irrep and number of states to average over Note: Molcas
  #  cannot average over states of different symmetry (unlike Molpro). Defaults
  #  to the lowest spin totally symmetric state (as in Molcas).
  my $spin = ($nprotons - $nelectrons) % 2 + 1;
  my $irrep = 1;
  my $ciroot = [1, 1, 1];

  %spin_labels = reverse %spin_multiplicity;
  my $found_states = 0;
  for my $i (1..9){
    if ($found_states) {last;}
    if (defined $r_par->{'model'}->{'ncasscf_states'}->{$spin_labels{$i}}){
      for my $isym (0..$nir-1){
        if ($r_par->{'model'}->{'ncasscf_states'}->{$spin_labels{$i}}->[$isym] > 0){
          my $nstates = $r_par->{'model'}->{'ncasscf_states'}->{$spin_labels{$i}}->[$isym];
          $spin = $i;
          $irrep = MultiSpace::map_to_qchem_irrep_order($isym, $point_group, "molcas") + 1;
          $ciroot = [$nstates, $nstates, 1];
          $found_states = 1;
          last;
        }
      }
    }
  }

  my %molcas_input;
  if ($r_par->{'model'}->{'model'} =~ /^SE/ || $r_par->{'data'}->{'scf_ok'} == 0 || $r_par->{'model'}->{'orbitals'} eq "HF") {
    $molcas_input{'NACTEL'} = [0,0,0];
    $molcas_input{'FROZEN'} = [0,0,0,0,0,0,0,0];
    $molcas_input{'INACTIVE'} = [0,0,0,0,0,0,0,0];
    $molcas_input{'GASSCF'} = [[0], [0,0,0,0,0,0,0,0], [0,0]];
  } else {
    my $target_space = $r_par->{'data'}{'MAS'}{'qchem'};
    if ($target_space->{'type'} != 2){
      die "GAS approach must be used for Molcas GASSCF run!";
    }
    %molcas_input = $target_space->molcas_gasscf_input();
  }
  &replace_in_template($r_str, "IRREP",  $irrep);
  &replace_in_template($r_str, "SPIN",  $spin);
  &replace_in_template_with_arrays($r_str, "CIROOT",  $ciroot, ' ');
  &replace_in_template_with_arrays($r_str, "NACTEL",  $molcas_input{'NACTEL'}, ' ');
  &replace_in_template_with_arrays($r_str, "FROZEN",  $molcas_input{'FROZEN'}, ' ');
  &replace_in_template_with_arrays($r_str, "INACTIVE",  $molcas_input{'INACTIVE'}, ' ');
  &replace_in_template_with_arrays($r_str, "GASSCF",  $molcas_input{'GASSCF'}, ' ');

  return 1;
}

# ------------------------- integrals -------------------------

sub make_scatci_integrals_input {
  my ($r_par, $r_str) = @_;

### TODO: for bound run set RMATR < 0 !!!
  if ($r_par->{'run'}->{'bound'} == 1) {
    &replace_in_template($r_str, "RMATR", "-1.0");
    &replace_in_template($r_str, "RUN_FREE_SCATTERING", ".false.");
  }
  else {
    &replace_in_template($r_str, "RMATR",    $r_par->{'model'}->{'rmatrix_radius'});
    &replace_in_template($r_str, "RUN_FREE_SCATTERING", ".true.");
  }
  &replace_in_template($r_str, "NSYMOP",   $symmetries{$r_geom->{'symmetry'}}->[0]);
  &replace_in_template($r_str, "SYMOP",    $symmetries{$r_geom->{'symmetry'}}->[1]);
  &replace_in_template($r_str, "MOLECULE", lc($r_par->{'model'}->{'molecule'}));
  &replace_in_template($r_str, "NOB", join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}}));
  if ($r_par->{'model'}->{'select_orb_by'} eq "molden") {
    &replace_in_template($r_str, "SELECT_ORB_BY", 1);
  }
  else {
    &replace_in_template($r_str, "SELECT_ORB_BY", 2);
  }

  my $str_PCO = "";
  if ($r_par->{'model'}->{'use_PCO'}) { #Adding PCOs
    my $nir = $r_par->{'data'}->{'nir'};
    my $l = $r_par->{'model'}->{'maxl_PCO'};
    $str_PCO = "&pco_data \n";
    $str_PCO .= "  min_PCO_l = 0, max_PCO_l = $l,\n";
    $str_PCO .= "  PCO_alpha0(0:$l) = ".join(",", @{$r_par->{'model'}->{'PCO_alpha0'}}[0..$l]).",\n";
    $str_PCO .= "  PCO_beta(0:$l) = ".join(",", @{$r_par->{'model'}->{'PCO_beta'}}[0..$l]).",\n";
    $str_PCO .= "  PCO_gto_thrs(0:$l) = ".join(",", @{$r_par->{'model'}->{'PCO_gto_thrs'}}[0..$l]).",\n";
    $str_PCO .= "  num_PCOs(0:$l) = ".join(",", @{$r_par->{'model'}->{'num_PCOs'}}[0..$l]).",\n";
    $str_PCO .= "  PCO_del_thrs = ".join(",", @{$r_par->{'model'}->{'PCO_delthres'}}[0..$nir-1]).",\n";
    $str_PCO .= "&end";
  }
  &replace_in_template($r_str, "PCO_BASIS", $str_PCO);
  my $nir = $r_par->{'data'}->{'nir'};
  &replace_in_template($r_str, "DELTHRES", join(",", @{$r_par->{'model'}->{'delthres'}}[0..$nir-1]));

  # add Gaussian basis if wanted
  my $str_GTO = "";
  if ($r_par->{'run'}->{'scattering'} == 1 and $r_par->{'model'}->{'use_GTO'}) {
    $str_GTO .= "min_l = 0, max_l = ".sprintf("%d", $r_par->{'model'}->{'maxl_GTO'}).",\n";
    $str_GTO .= "  ".&get_exponents_from_swmol3_continuum_basis($r_par).",\n\n";
  }
  &replace_in_template($r_str, "GTO_BASIS", $str_GTO);

  # add B-spline basis if wanted
  my $str_BTO = "";
  if ($r_par->{'run'}->{'scattering'} == 1 and $r_par->{'model'}->{'use_BTO'}) {
    $str_BTO .=   "bspline_grid_start = ".sprintf("%.3f", $r_par->{'model'}->{'start_BTO'}).",\n";
    $str_BTO .= "  bspline_order = ".sprintf("%d", $r_par->{'model'}->{'order_BTO'}).",\n";
    $str_BTO .= "  no_bsplines = ".sprintf("%d", $r_par->{'model'}->{'no_of_BTO'}).",\n\n";
    $str_BTO .= "  min_bspline_l = 0, max_bspline_l = ".sprintf("%d", $r_par->{'model'}->{'maxl_BTO'}).",\n";
    my $first_bspline = 3;
    if ($r_par->{'model'}->{'start_BTO'} <= 0.001) { $first_bspline = 2; }
    for (my $l = 0; $l <= $r_par->{'model'}->{'maxl_BTO'}; $l++) {
      $str_BTO .= "  bspline_indices(1,$l) = $first_bspline,\n";
      $str_BTO .= "  bspline_indices(2,$l) = $r_par->{'model'}->{'no_of_BTO'},\n\n";
    }
  }
  &replace_in_template($r_str, "BTO_BASIS", $str_BTO);

  &replace_in_template($r_str, "BUFFER_SIZE", $r_par->{'run'}->{'buffer_size'});
  &replace_in_template($r_str, "INTEGRAL_TRANS_ALG", $r_par->{'run'}->{'transform_alg'});
  &replace_in_template($r_str, "DELTAR1", $r_par->{'run'}->{'delta_r1'});
  &replace_in_template($r_str, "MAXL_LEGENDRE_1EL", $r_par->{'model'}->{'maxl_legendre_1el'});
  &replace_in_template($r_str, "MAXL_LEGENDRE_2EL", $r_par->{'model'}->{'maxl_legendre_2el'});

  return 1;
}


sub get_exponents_from_swmol3_continuum_basis {
  my ($r_par) = @_;
  my $str_exp = "";
  my $str_basis = &get_basis_for_swmol3_input($r_par->{'dirs'}->{'basis'}, $r_par->{'data'}->{'scattering'}->{'basis'}, "continuum",  0.0, 0.0, 0.0, $r_par);
  my @nexp = ();
  if ($str_basis =~ m/jco = ((\d+,)+)\s/s) {
    @nexp = split(',', $1);
  }
  else { die "Error: failed to read 'jco' from continuum basis file !\n"; }
  my @all_exp = split(/\s+1\s+1\s+/s, $str_basis);  # the first element is &atom namelist
  foreach (@all_exp) { $_ =~ s/\s+1\.\s*$// ; }     # but other elements are only exponents
  my $ne = 0;
  for (my $n = 0; $n < scalar @nexp; $n++) {
    $str_exp .= "exponents(:,$n) = ";
    for (my $i = 1; $i <= $nexp[$n]; $i++) {
      $ne++;
      $str_exp .= "$all_exp[$ne], ";
    }
    if ($n < (scalar @nexp) - 1) { $str_exp .= "\n  "; }
  }
  return $str_exp;
}


# ------------- congen - multiple active spaces -------------

sub make_congen_mas_input {
  my ($r_par, $r_str) = @_;

  my $task        = $r_par->{'data'}->{'task'};
  my $statesym    = $r_par->{'data'}->{$task}->{'symmetry'};                      # number 0,1,2,...
  my $stateir     = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym];    # string A1,B1,...
  my $statespin   = $r_par->{'data'}->{$task}->{'spin'};                          # string singlet,doublet,...
  my $statespinno = $spin_multiplicity{lc($statespin)};                           # number 1,2,3,...
  my $nelectrons  = $r_par->{'model'}->{'nelectrons'};                            # number of electrons, first target
  my $n_irreps = scalar @{$irred_repr{$r_par->{'model'}->{'symmetry'}}};


  # Get the target space (and l2 space if given)
  # ==============================================

  my @nob0 = ();

  my $target_space = $r_par->{'data'}{'MAS'}{'target'};
  my @total_target_used = $target_space->total_orbitals();

  my $l2_space;
  my @total_l2_used = (0) x 8;
  if (@{$r_par->{'model'}->{'l2_MAS'}}){
    $l2_space = $r_par->{'data'}{'MAS'}{'l2'};
    @total_l2_used = $l2_space->total_orbitals();
  };

  # If not used in the target or l2 space then virtual orbitals are considered to
  # be part of the continuum, and will be contracted with them.
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    $nob0[$i] = max($total_target_used[$i], $total_l2_used[$i]);
  }

  # We need to add the virtuals to the continuum
  $r_par->{'data'}->{'orbitals'}->{'virt_cont_used'} = [];
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] - $nob0[$i];
    $r_par->{'data'}->{'orbitals'}->{'virt_cont_used'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$i] + $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i];
  }


  # Start - Taken directly from make_congen_input
  #----------------------------------------------
  # &state namelist
  #================

  my %state_namelist = (
    'MOLECULE', $r_par->{'model'}->{'molecule'},           # For scattering: "e + " will be added
    'SPIN',     $statespin,
    'SYMMETRY', $stateir,
    'MEGUL',    "$lucsf_base$statespinno$statesym",           # For target: e.g. 713, 1 for singlet, 3 for symmetry
    'ISCAT',    1,                                            # For scattering: changed to 2
    'SYMTYP',   2,
    'QNTOT',    "$statespinno,$statesym,0",
    'NELECT',   $nelectrons,                                  # For scattering: +1 will be added
    'IPOSIT',   0, # "", # flag for positron scattering # TODO: Change when positrons implemented in MAS
    'NOB',      "",
    'NOBE',     "",
    'NOBP',     "",
    'NOBV',     "",
    'NOB0',     "",
    'NREFO',    0,
    'REFORB',   "",
    'LNDO',     $r_par->{'model'}->{'lndo'},
    'GFLAG',    1,
    'POSITRON_FLAG',  "!",
  );

  &print_info("\nSearching for reference orbitals for $statespin $stateir state.\n", $r_par);
  my $nrefo  = 0;
  my @reforb = ();

  # Special settings for target
  if ($task eq "target") {
    if ($r_par->{'model'}->{'model'} eq "CHF-A") {
      my $megul = sprintf('%03s', $r_par->{'data'}->{'CHF'}->{'current_state'});
      $state_namelist{'MEGUL'} = $lucsf_base.$megul;
    }
    else{
      $state_namelist{'MEGUL'} = $lucsf_base.$statespinno.$statesym;
    }
    $state_namelist{'NOB'}   = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all'}});
    $state_namelist{'NOB0'}  = join(",", @nob0);

    # reference orbitals depend on target state spin and symmetry
    # only frozen + active orbitals are used to find reference orbital
    &print_info("Putting $nelectrons electrons into (".join(',', @{$r_par->{'data'}->{'orbitals'}->{'reference'}}).") orbitals\n", $r_par);
    ($nrefo, @reforb) = &get_reference_orbitals($r_par, $nelectrons, $statesym, $statespinno);
    if ($nrefo == 0) {
      die   "Increase your active space or change the target settings:\n
             put 0 for target state: $statespin $stateir($statesym) !\n";
    }
  }

  # Special settings for scattering calculation
  else {

    $nelectrons++;
    $state_namelist{'MOLECULE'} = "e + $state_namelist{'MOLECULE'}";
    $state_namelist{'ISCAT'} = 2;
    $state_namelist{'NOB'}  = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}});
    # $state_namelist{'NOB0'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}});
    $state_namelist{'NOB0'}  = join(",", @nob0);
    if ($r_par->{'model'}->{'model'} eq "CHF-A") {
      $state_namelist{'GFLAG'} = 0;
    }

    # reference orbitals depend on target state spin and symmetry
    # now all available orbitals are used
    &print_info("Putting $nelectrons electrons into (".join(',', @{$r_par->{'data'}->{'orbitals'}->{'reference'}}).") orbitals\n", $r_par);
    ($nrefo, @reforb) = &get_reference_orbitals($r_par, $nelectrons, $statesym, $statespinno);
    if ($nrefo == 0) {
      die   "Increase the number of orbitals for searching reference orbitals !\n";
    }
  }

  $state_namelist{'NELECT'} = $nelectrons;
  $state_namelist{'NREFO'} = $nrefo;
  $state_namelist{'REFORB'} = join(",", @reforb);
  $state_namelist{'REFORB'} =~ s/(\d+,\d+,\d+,\d+,\d+,)/\1 /g; # Add spaces for better reading
  &print_info("Reference orbitals found:\n$state_namelist{'REFORB'}\n\n", $r_par);

  &replace_all_in_template($r_str, \%state_namelist);

  # End - Taken directly from make_congen_input
  #----------------------------------------------

  if ($task ne "target") {
    $nelectrons--;
  }

  # &wfngrp namelist(s)
  #====================

  my %wfngrp_namelist = (
    'GNAME',    "",
    'QNTAR',    "-1,0,0", # Default for target, no constraints
    'NELECG',   $r_par->{'model'}->{'nelectrons'},         # For scattering: +1
    'NREFOG',   $state_namelist{'NREFO'},
    'REFORG',   $state_namelist{'REFORB'},
    'NDPROD',   0,
    'NELECP',   "",
    'NSHLP',    "",
    'PQN',      "",
    'MSHL',     "",
    'LSQUARE',   1,
    'L2_POSITRON_FLAG',  "!",  #, "",   TODO: This should be changed when positron
                               #        scattering is implemented in MAS.
  );

  # Read &wfngrp namelist template
  my $wfngrp_template = "";
  if(!&read_file("$r_par->{'dirs'}->{'templates'}${bs}congen.wfngrp.inp", \$wfngrp_template)) {
    &print_info("Warning: no template file for congen.wfngrp.inp !\n", $r_par);
  }

  if ($task eq "target") {

    my @wfngrps = $target_space->congen_target($nelectrons);

    foreach my $wg (@wfngrps){
      my $wg_namelist = dclone(\%wfngrp_namelist);
      my %updated_wfngrp = (%{$wg_namelist}, %{$wg});
      $wfngrp = $wfngrp_template;
      &replace_all_in_template_with_arrays(\$wfngrp, \%updated_wfngrp, ',');
      $$r_str .= $wfngrp;
    };

  } else {
      $wfngrp_namelist{'NELECG'} += 1;
      my $total_spinsym  = [$statespinno, $statesym];
      # my @nob0 = @{$r_par->{'data'}->{'orbitals'}->{'target_used'}};

      my @all_target_spinsyms = ();
      for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
        foreach my $targetstatespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
          if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$spin_multiplicity{lc($targetstatespin)}.".$i"}) {
              push(@all_target_spinsyms, [$spin_multiplicity{lc($targetstatespin)}, $i]);
          }
        }
      }

      my @wfngrps = ();
      if (@{$r_par->{'model'}->{'l2_MAS'}}){
        @wfngrps = $target_space->congen_scattering($nelectrons, \@nob0, $total_spinsym, \@all_target_spinsyms, $l2_space);
      } else {
        @wfngrps = $target_space->congen_scattering($nelectrons, \@nob0, $total_spinsym, \@all_target_spinsyms);
      }

      foreach my $wg (@wfngrps){
        my $wg_namelist = dclone(\%wfngrp_namelist);
        my %updated_wfngrp = (%{$wg_namelist}, %{$wg});
        $wfngrp = $wfngrp_template;
        &replace_all_in_template_with_arrays(\$wfngrp, \%updated_wfngrp, ',');
        $$r_str .= $wfngrp;
      };
  }
  return 1;
}
# ------------------------ congen ---------------------------

sub make_congen_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym    = $r_par->{'data'}->{$task}->{'symmetry'};                      # number 0,1,2,...
  my $stateir     = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym];    # string A1,B1,...
  my $statespin   = $r_par->{'data'}->{$task}->{'spin'};                          # string singlet,doublet,...
  my $statespinno = $spin_multiplicity{lc($statespin)};                           # number 1,2,3,...
  my $nelectrons  = $r_par->{'model'}->{'nelectrons'};                            # number of electrons, first target
  my $nactive_electrons = $nelectrons - 2 * $r_par->{'model'}->{'nfrozen'};       # number of target active electrons
  my $iposit = -$r_par->{'model'}->{'positron_flag'};
  my $pos_scattering = "";

  # &state namelist
  #================

  my %state_namelist = (
    'MOLECULE', $r_par->{'model'}->{'molecule'},           # For scattering: "e + " will be added
    'SPIN',     $statespin,
    'SYMMETRY', $stateir,
    'ISCAT',    1,                                            # For scattering: changed to 2
    'SYMTYP',   2,
    'QNTOT',    "$statespinno,$statesym,0",
    'NELECT',   $nelectrons,                                  # For scattering: +1 will be added
    'IPOSIT',   "",                               # flag for positron scattering
    'NOB',      "",
    'NOBE',     "",
    'NOBP',     "",
    'NOBV',     "",
    'NOB0',     "",
    'NREFO',    0,
    'REFORB',   "",
    'LNDO',     $r_par->{'model'}->{'lndo'},
    'GFLAG',    1,
    'POSITRON_FLAG',  "!",
  );

  &print_info("\nSearching for reference orbitals for $statespin $stateir state.\n", $r_par);
  my $nrefo  = 0;
  my @reforb = ();

  # Special settings for target
  if ($task eq "target") {
    $iposit = 0;
    if ($r_par->{'model'}->{'model'} eq "CHF-A") {
      my $megul = sprintf('%03s', $r_par->{'data'}->{'CHF'}->{'current_state'});
      $state_namelist{'MEGUL'} = $lucsf_base.$megul;
    }
    else{
      $state_namelist{'MEGUL'} = $lucsf_base.$statespinno.$statesym;
    }
    $state_namelist{'NOBE'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all'}});
    $state_namelist{'NOBP'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'empty'}});
    $state_namelist{'NOB'}  = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all'}});
    $state_namelist{'NOBV'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}});
    $state_namelist{'NOB0'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}});

    # reference orbitals depend on target state spin and symmetry
    # only frozen + active orbitals are used to find reference orbital
    &print_info("Putting $nelectrons electrons into (".join(',', @{$r_par->{'data'}->{'orbitals'}->{'reference'}}).") orbitals\n", $r_par);
    ($nrefo, @reforb) = &get_reference_orbitals($r_par, $nelectrons, $statesym, $statespinno);
    if ($nrefo == 0) {
      die   "Increase your active space or change the target settings:\n
             put 0 for target state: $statespin $stateir($statesym) !\n";
    }
  }

  # Special settings for scattering calculation
  else {

    $nelectrons++;
    $state_namelist{'MOLECULE'} = "e + $state_namelist{'MOLECULE'}";
    $iposit = -$r_par->{'model'}->{'positron_flag'};
    $state_namelist{'IPOSIT'} = $iposit;
    if ($iposit == -1) {$state_namelist{'POSITRON_FLAG'} = $pos_scattering;}
    $state_namelist{'ISCAT'} = 2;
    $state_namelist{'MEGUL'} = $lucsf_base.$statespinno.$statesym;
    if ($iposit == 0) {
      $state_namelist{'NOB'}  = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}});
      $state_namelist{'NOBE'}  = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}});
      $state_namelist{'NOBP'}   = join(",", @{$r_par->{'data'}->{'orbitals'}->{'empty'}});
      $state_namelist{'NOB0'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}});
      $state_namelist{'NOBV'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}});
    }
    else {
      $state_namelist{'NOBE'}  = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}});
      $state_namelist{'NOBP'}  = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}});
      $state_namelist{'NOB'}  = join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used_positron'}});
      $state_namelist{'NOB0'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}});
      $state_namelist{'NOBV'} = join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}});
    }
    if ($r_par->{'model'}->{'model'} eq "CHF-A") {
      $state_namelist{'GFLAG'} = 0;
    }

    # reference orbitals depend on target state spin and symmetry
    # now all available orbitals are used
    &print_info("Putting $nelectrons electrons into (".join(',', @{$r_par->{'data'}->{'orbitals'}->{'reference'}}).") orbitals\n", $r_par);
    ($nrefo, @reforb) = &get_reference_orbitals($r_par, $nelectrons, $statesym, $statespinno);
    if ($nrefo == 0) {
      die   "Increase the number of orbitals for searching reference orbitals !\n";
    }
  }

  $state_namelist{'NELECT'} = $nelectrons;
  $state_namelist{'IPOSIT'} = $iposit;
  $state_namelist{'NREFO'} = $nrefo;
  $state_namelist{'REFORB'} = join(",", @reforb);
  $state_namelist{'REFORB'} =~ s/(\d+,\d+,\d+,\d+,\d+,)/\1 /g; # Add spaces for better reading
  &print_info("Reference orbitals found:\n$state_namelist{'REFORB'}\n\n", $r_par);

  &replace_all_in_template($r_str, \%state_namelist);

  # &wfngrp namelist(s)
  #====================

  my %wfngrp_namelist = (
    'GNAME',    "",
    'QNTAR',    "-1,0,0", # Default for target, no constraints
    'NELECG',   $r_par->{'model'}->{'nelectrons'},         # For scattering: +1
    'NREFOG',   $state_namelist{'NREFO'},
    'REFORG',   $state_namelist{'REFORB'},
    'NDPROD',   0,
    'NELECP',   "",
    'NSHLP',    "",
    'PQN',      "",
    'MSHL',     "",
    'LSQUARE',   1,
    'L2_POSITRON_FLAG',  "",
  );

  # Read &wfngrp namelist template
  my $wfngrp_template = "";
  if(!&read_file("$r_par->{'dirs'}->{'templates'}${bs}congen.wfngrp.inp", \$wfngrp_template)) {
    &print_info("Warning: no template file for congen.wfngrp.inp !\n", $r_par);
  }
  my $wfngrp = $wfngrp_template;

  # For target we add one &wfngrp namelist
  # it is assumed that information about frozen orbitals is
  #   in the array $r_par->{'data'}->{'orbitals'}->{'frozen'}
  # and about active (= valence for SE and SEP) orbitals are
  #   in the array $r_par->{'data'}->{'orbitals'}->{'active'}
  if ($task eq "target") {
    $wfngrp_namelist{'L2_POSITRON_FLAG'} = "!";

    if ($r_par->{'model'}->{'model'} =~ /(^CHF)/) {

      my @states = @{$r_par->{'data'}->{'CHF'}->{'states'}};
      my $nactive = $r_par->{'model'}->{'nactive'};

      for $state (@states){
        my $current_id = $r_par->{'data'}->{'CHF'}->{'current_state'};
        my $id = $state->{'id'};
        my $sym = $state->{'sym'};
        my $mult = $state->{'multiplet'};
        my @config = @{$state->{'config'}};
        if (($sym eq $statesym) && ($mult eq $statespin)) {

          next if ($r_par->{'model'}->{'model'} eq "CHF-A" && $id ne $current_id);

          # photoionization calculation (HF-like targets with single excitation)
          # electron-scattering calculation (ground state HF)
          if (scalar @config <= $nactive) {
            $wfngrp_namelist{'GNAME'} = "HF $id";
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            for my $shell (@config){
              my ($shell_sym, $shell_num, $ne) = split /,/, $shell;
              &add_active_shell_for_CHF($r_par, \%wfngrp_namelist, $shell_sym, $shell_num, $ne);
            }
          }
          # -----------CHF configurations-----------
          else {
            $wfngrp_namelist{'GNAME'} = "CHF $id";
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            for(my $i = 0; $i < $nactive; $i++){
              my ($shell_sym, $shell_num, $ne) = split /,/, @config[$i];
              &add_active_shell_for_CHF($r_par, \%wfngrp_namelist, $shell_sym, $shell_num, $ne);
            }
            my ($shell_sym, $shell_num, $ne) = split /,/, @config[$nactive];
            &add_virtual_shell_for_CHF($r_par, \%wfngrp_namelist, $shell_sym, $shell_num);
          }
          $wfngrp = $wfngrp_template;
          &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
          $$r_str .= $wfngrp;
        }
      }
    }
    else {
      $wfngrp_namelist{'GNAME'} = "$state_namelist{'MOLECULE'} - $state_namelist{'SPIN'} $state_namelist{'SYMMETRY'}";
      &add_frozen_shells($r_par, \%wfngrp_namelist);
      &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
      &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
      $$r_str .= $wfngrp;

      if ( $r_par->{'model'}->{'use_PCO'} && $r_par->{'model'}->{'model'} eq "CAS") {
        if ( $r_par->{'model'}->{'reduce_PCO_CAS'} ) {
          $wfngrp_namelist{'GNAME'} = "$state_namelist{'MOLECULE'} - $state_namelist{'SPIN'} $state_namelist{'SYMMETRY'} GS^N-1 PCO^1";
          &add_GS_shells($r_par, \%wfngrp_namelist, $r_par->{'model'}->{'nelectrons'}-1 );
          &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
          $wfngrp = $wfngrp_template;
          &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
          $$r_str .= $wfngrp;
        }else{
          $wfngrp_namelist{'GNAME'} = "$state_namelist{'MOLECULE'} - $state_namelist{'SPIN'} $state_namelist{'SYMMETRY'} Core^Nc CAS^N-Nc-1 PCO^1";
          &add_frozen_shells($r_par, \%wfngrp_namelist);
          &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons-1);
          &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
          $wfngrp = $wfngrp_template;
          &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
          $$r_str .= $wfngrp;
        }
      } # For PCOs
    }

  }
  # For scattering calculation the number of &wfngrp namelists
  # depends on model
  else {
    $wfngrp_namelist{'NELECG'}++;
    if ($iposit == 0) {
      $wfngrp_namelist{'L2_POSITRON_FLAG'} = "!";
    }

    # SE and SEP models
    if ($r_par->{'model'}->{'model'} =~ /(SE|SEP)/) {

      # Determine spin and symmetry of the ground state
      my ($ground_state_spin, $ground_state_sym) = split(/\./, $r_par->{'data'}->{'target'}->{'ground_state'});

      # In both SE and SEP we use
      # ground state x continuum^1
      $wfngrp_namelist{'GNAME'} = "Ground state x continuum^1";
      $wfngrp_namelist{'L2_POSITRON_FLAG'} = '!';
      $wfngrp_namelist{'LSQUARE'} = 0;
      $wfngrp_namelist{'QNTAR'} = "$ground_state_spin,$ground_state_sym,0";
      &add_frozen_shells($r_par, \%wfngrp_namelist);
      &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
      &add_continuum_shells($r_par, \%wfngrp_namelist, $ground_state_sym, $iposit);
      $wfngrp = $wfngrp_template;
      &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
      $$r_str .= $wfngrp;

      #IF DOING PCCHF - MANY TARGET STATES
      print($r_par->{'model'}->{'ntarget_states_used'});
      if ($r_par->{'model'}->{'ntarget_states_used'} != 1) {
        for (my $i =0; $i < $r_par->{'data'}->{'nir'}; $i++) {
          foreach my $targetstatespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
            if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$spin_multiplicity{lc($targetstatespin)}.".$i"}) { # skip if target states of a given spin-symmetry are not used
              $found = 0;
              while ($found == 0) {
                for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
                  if ($group_table[$i]->[$j] == $statesym) {
                    $found = 1;
                    #if using virtual orbitals then:
                    if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$j] > 0) {
                      my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i];
                      $wfngrp_namelist{'GNAME'} = "$str_ir (GS^N-1 virtual^1) x continuum^1";
                      $wfngrp_namelist{'L2_POSITRON_FLAG'} = '!';
                      $wfngrp_namelist{'LSQUARE'} = 0;
                      $wfngrp_namelist{'QNTAR'} = "$spin_multiplicity{lc($targetstatespin)},$i,0";
                      &add_frozen_shells($r_par, \%wfngrp_namelist);
                      &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
                      &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
                      &add_continuum_shells($r_par, \%wfngrp_namelist, $j, $iposit);
                      $wfngrp = $wfngrp_template;
                      &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                      $$r_str .= $wfngrp;
                    }
                  }
                }
              }
            }
          }
        }
      }

      if ( $r_par->{'model'}->{'use_PCO'} ) {
        if ($r_par->{'model'}->{'model'} eq "SEP") {
          $wfngrp_namelist{'GNAME'} = "( GS^N-1 PCO^1 ) x continuum^1";
          $wfngrp_namelist{'L2_POSITRON_FLAG'} = '!';
          $wfngrp_namelist{'LSQUARE'} = 0;
          $wfngrp_namelist{'QNTAR'} = "$ground_state_spin,$ground_state_sym,0";
          &add_frozen_shells($r_par, \%wfngrp_namelist);
          &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons-1);
          &add_PCO_shell_for_given_target_symmetry($r_par, \%wfngrp_namelist, $ground_state_sym, 1, 0);
          &add_continuum_shells($r_par, \%wfngrp_namelist, $ground_state_sym, $iposit);
          $wfngrp = $wfngrp_template;
          &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);

          $$r_str .= $wfngrp;
          $wfngrp = $wfngrp_template;
          if ($iposit == 0) {
            $wfngrp_namelist{'GNAME'} = "GS^N-1 PCO^2";
            $wfngrp_namelist{'QNTAR'} = "-1,0,0";
            $wfngrp_namelist{'LSQUARE'} = 1;
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons-1);
            &add_PCO_shells($r_par, \%wfngrp_namelist, 2);
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
          else {
            for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
              $found = 0;
              while ($found == 0) {
                for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
                  if ($group_table[$i]->[$j] == $statesym) {
                    $found = 1;
                    next if ($r_par->{'data'}->{'orbitals'}->{'PCO'}->[$j] <= 0);
                    $wfngrp_namelist{'GNAME'} = "GS^N-1 PCO^1 positron-PCO^1";
                    $wfngrp_namelist{'QNTAR'} = "$ground_state_spin,$i,0";
                    $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
                    $wfngrp_namelist{'LSQUARE'} = 1;
                    &add_frozen_shells($r_par, \%wfngrp_namelist);
                    &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons-1);
                    &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
                    &add_PCO_shell_for_given_target_symmetry($r_par, \%wfngrp_namelist, $j, 1, $iposit);
                    $wfngrp = $wfngrp_template;
                    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                    $$r_str .= $wfngrp;
                  }
                }
              }
            }
          }
        }
        # GS^N PCO^1

        $wfngrp = $wfngrp_template;
        if ($iposit == 0) {
          $wfngrp_namelist{'GNAME'} = "GS^N PCO^1";
          $wfngrp_namelist{'QNTAR'} = "-1,0,0";
          $wfngrp_namelist{'LSQUARE'} = 1;
          &add_frozen_shells($r_par, \%wfngrp_namelist);
          &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
          &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
        }
        else {
          $wfngrp_namelist{'GNAME'} = "GS^N positron-TGT-to-PCO^1";
          $wfngrp_namelist{'QNTAR'} = "$ground_state_spin,$ground_state_sym,0";
          $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
          $wfngrp_namelist{'LSQUARE'} = 1;
          &add_frozen_shells($r_par, \%wfngrp_namelist);
          &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
          &add_positron_TGTPCO_shell_for_given_target_symmetry($r_par, \%wfngrp_namelist, $statesym);
        }
        &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
        $$r_str .= $wfngrp;
        $wfngrp = $wfngrp_template;

      } # For PCOs
      # ground state x virtual^1 (if any)
      if ($r_par->{'data'}->{'orbitals'}->{'virtual'}->[$statesym] > 0) {
        if ($iposit == 0) {
          $wfngrp_namelist{'GNAME'}  = "Ground state x virtual^1";
          $wfngrp_namelist{'QNTAR'} = "-1,0,0";
          $wfngrp_namelist{'LSQUARE'} = 1;
          &add_frozen_shells($r_par, \%wfngrp_namelist);
          &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
          &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
          $wfngrp = $wfngrp_template;
          &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
          $$r_str .= $wfngrp;
        }
      }
      if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$statesym] > 0) {
        if ($iposit == -1) {
          if ($r_par->{'model'}->{'use_PCO'} == 0) { #if using PCOs then these configurations are included in GS^N positron-TGT-to-PCO^1, therefore skip
            $wfngrp_namelist{'GNAME'}  = "Ground state x positron-TGT-to-virtual^1";
            $wfngrp_namelist{'QNTAR'} = "1,0,0";
            $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
            $wfngrp_namelist{'LSQUARE'} = 1;
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
            &add_positron_frozact_shells($r_par, \%wfngrp_namelist, $statesym);
            $wfngrp = $wfngrp_template;
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
        }
      }

      # and for SEP only
      if ($r_par->{'model'}->{'model'} eq "SEP") {
        # core valence^-1 x virtual^2
        if ($r_par->{'data'}->{'orbitals'}->{'virtual'}->[$statesym] > 0) {
          if ($iposit == 0) {
            $wfngrp_namelist{'GNAME'}  = "Core valence^-1 x virtual^2";
            $wfngrp_namelist{'LSQUARE'} = 1;
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
            &add_virtual_shells($r_par, \%wfngrp_namelist, 2);
            $wfngrp = $wfngrp_template;
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
        }
        if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$statesym] > 0) {
          if ($iposit == -1) {
            for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
              $found = 0;
              while ($found == 0) {
                for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
                  if ($group_table[$i]->[$j] == $statesym) {
                    $found = 1;
                    next if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$j] <=0);
                    $wfngrp_namelist{'GNAME'}  = "Core valence^-1 x virtual^1 x positron-TGT-to-virtual^1";
                    $wfngrp_namelist{'QNTAR'} = "1,$i,0";
                    $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
                    $wfngrp_namelist{'LSQUARE'} = 1;
                    &add_frozen_shells($r_par, \%wfngrp_namelist);
                    &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
                    &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
                    &add_positron_frozact_shells($r_par, \%wfngrp_namelist, $j);
                    $wfngrp = $wfngrp_template;
                    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                    $$r_str .= $wfngrp;
                  }
                }
              }
            }
          }
        }
      }
      # for SE, SEP (Photoionization only)
      if ($r_par->{'run'}->{'photoionization'} == 1) {
        $wfngrp_namelist{'GNAME'}  = "Ground state";
        $wfngrp_namelist{'QNTAR'} = "-1,0,0";
        &add_frozen_shells($r_par, \%wfngrp_namelist);
        &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons + 1);
        $wfngrp = $wfngrp_template;
        &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
        $$r_str .= $wfngrp;
      }

    }
    elsif ($r_par->{'model'}->{'model'} =~ /CHF/) {

      my @states = @{$r_par->{'data'}->{'CHF'}->{'states'}};
      my $nactive = $r_par->{'model'}->{'nactive'};

      for $state (@states){
        my $id = $state->{'id'};
        my $sym = $state->{'sym'};
        my $mult = $spin_multiplicity{lc($state->{'multiplet'})};
        my @config = @{$state->{'config'}};
        &add_frozen_shells($r_par, \%wfngrp_namelist);

        # -----------SE-like configurations (photoionization)-----------
        if (scalar @config <= $nactive) {
          $wfngrp_namelist{'GNAME'} = "HF $id x cont^1";
          for my $shell (@config){
            my ($shell_sym, $shell_num, $ne) = split /,/, $shell;
            &add_active_shell_for_CHF($r_par, \%wfngrp_namelist, $shell_sym, $shell_num, $ne);
          }
        }
        # -----------CHF configurations-----------
        else {
          $wfngrp_namelist{'GNAME'} = "CHF $id x cont^1";
          for(my $i = 0; $i < $nactive; $i++){
            my ($shell_sym, $shell_num, $ne) = split /,/, @config[$i];
            &add_active_shell_for_CHF($r_par, \%wfngrp_namelist, $shell_sym, $shell_num, $ne);
          }
          my ($shell_sym, $shell_num, $ne) = split /,/, @config[$nactive];
          &add_virtual_shell_for_CHF($r_par, \%wfngrp_namelist, $shell_sym, $shell_num)
        }
        # couple all target states to appropriate symmetry
        $wfngrp_namelist{'L2_POSITRON_FLAG'} = '!';
        $wfngrp_namelist{'LSQUARE'} = 0;
        $wfngrp_namelist{'QNTAR'} = "$mult,$sym,0";
        &add_continuum_shells($r_par, \%wfngrp_namelist, $sym, $iposit);
        $wfngrp = $wfngrp_template;
        &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
        $$r_str .= $wfngrp;
      }

      # ground state x virtual^1 (if any)
      if ($iposit == 0) {
        $wfngrp_namelist{'GNAME'}  = "Ground state x virtual^1";
        $wfngrp_namelist{'QNTAR'} = "-1,0,0";
        &add_frozen_shells($r_par, \%wfngrp_namelist);
        &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
        &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
        $wfngrp_namelist{'L2_POSITRON_FLAG'} = '!';
        $wfngrp_namelist{'LSQUARE'} = 1;
        $wfngrp = $wfngrp_template;
        &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
        $$r_str .= $wfngrp;

        # core valence^-1 x virtual^2
        $wfngrp_namelist{'GNAME'}  = "Core valence^-1 x virtual^2";
        &add_frozen_shells($r_par, \%wfngrp_namelist);
        &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
        &add_virtual_shells($r_par, \%wfngrp_namelist, 2);
        $wfngrp = $wfngrp_template;
        &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
        $$r_str .= $wfngrp;
      }
      else {
        if ($r_par->{'model'}->{'use_PCO'} == 0) { #if using PCOs then these configurations are included in GS^N positron-TGT-to-PCO^1, therefore skip
          $wfngrp_namelist{'GNAME'}  = "Ground state x positron-TGT-to-virtual^1";
          $wfngrp_namelist{'QNTAR'} = "1,0,0";
          $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
          $wfngrp_namelist{'LSQUARE'} = 1;
          &add_frozen_shells($r_par, \%wfngrp_namelist);
          &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
          &add_positron_frozact_shells($r_par, \%wfngrp_namelist, $statesym);
          $wfngrp = $wfngrp_template;
          &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
          $$r_str .= $wfngrp;
        }
        for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
          $found = 0;
          while ($found == 0) {
            for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
              if ($group_table[$i]->[$j] == $statesym) {
                $found = 1;
                next if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$j] <=0);
                $wfngrp_namelist{'GNAME'}  = "Core valence^-1 x virtual^1 x positron-TGT-to-virtual^1";
                $wfngrp_namelist{'QNTAR'} = "1,$i,0";
                $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
                $wfngrp_namelist{'LSQUARE'} = 1;
                &add_frozen_shells($r_par, \%wfngrp_namelist);
                &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
                &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
                &add_positron_frozact_shells($r_par, \%wfngrp_namelist, $j);
                $wfngrp = $wfngrp_template;
                &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                $$r_str .= $wfngrp;
              }
            }
          }
        }
      }

      # Photoionization only
      if ($r_par->{'run'}->{'photoionization'} == 1) {
        $wfngrp_namelist{'GNAME'}  = "Ground state";
        $wfngrp_namelist{'QNTAR'} = "-1,0,0";
        &add_frozen_shells($r_par, \%wfngrp_namelist);
        &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons + 1);
        $wfngrp = $wfngrp_template;
        &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
        $$r_str .= $wfngrp;
      }

    }
    # CAS models
    else {
      # First loop over all target states spin symmetries
      # and get Target state x continuum^1 wfngrp
      for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
        foreach my $targetstatespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
          if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$spin_multiplicity{lc($targetstatespin)}.".$i"}) { # skip if target states of a given spin-symmetry are not used
            my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i];
            $wfngrp_namelist{'GNAME'} = "$targetstatespin $str_ir x continuum^1";
            $wfngrp_namelist{'L2_POSITRON_FLAG'} = '!';
            $wfngrp_namelist{'LSQUARE'} = 0;
            $wfngrp_namelist{'QNTAR'} = "$spin_multiplicity{lc($targetstatespin)},$i,0";
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
            &add_continuum_shells($r_par, \%wfngrp_namelist, $i, $iposit);
            $wfngrp = $wfngrp_template;
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
            if ( $r_par->{'model'}->{'use_PCO'} ) {
              if ( $r_par->{'model'}->{'reduce_PCO_CAS'} ) {
                $wfngrp_namelist{'GNAME'} = "$targetstatespin $str_ir (GS^N-1 PCO^1) x continuum^1";
                $wfngrp_namelist{'L2_POSITRON_FLAG'} = '!';
                $wfngrp_namelist{'LSQUARE'} = 0;
                $wfngrp_namelist{'QNTAR'} = "$spin_multiplicity{lc($targetstatespin)},$i,0";
                &add_GS_shells($r_par, \%wfngrp_namelist, $r_par->{'model'}->{'nelectrons'}-1 );
                &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
                &add_continuum_shells($r_par, \%wfngrp_namelist, $i, $iposit);
                $wfngrp = $wfngrp_template;
                &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                $$r_str .= $wfngrp;
              }else{
                $wfngrp_namelist{'GNAME'} = "$targetstatespin $str_ir (core^Nc cas^N-Nc-1 PCO^1) x continuum^1";
                $wfngrp_namelist{'L2_POSITRON_FLAG'} = '!';
                $wfngrp_namelist{'LSQUARE'} = 0;
                $wfngrp_namelist{'QNTAR'} = "$spin_multiplicity{lc($targetstatespin)},$i,0";
                &add_frozen_shells($r_par, \%wfngrp_namelist);
                &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons-1);
                &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
                &add_continuum_shells($r_par, \%wfngrp_namelist, $i, $iposit);
                $wfngrp = $wfngrp_template;
                &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                $$r_str .= $wfngrp;
              }
            } # For PCOs
          }
        } # foreach spin state (multiplicity)
      } # for each IRs

      if ( $r_par->{'model'}->{'use_PCO'} ) {
        if ( $r_par->{'model'}->{'reduce_PCO_CAS'} ) {
          if ($iposit == 0) {
            $wfngrp = $wfngrp_template;
            $wfngrp_namelist{'LSQUARE'} = 1;
            $wfngrp_namelist{'GNAME'}  = "GS^N PCO^1";
            $wfngrp_namelist{'QNTAR'} = "-1,0,0";
            &add_GS_shells($r_par, \%wfngrp_namelist, $r_par->{'model'}->{'nelectrons'} );
            &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
	    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
          else {
            $wfngrp = $wfngrp_template;
            $wfngrp_namelist{'LSQUARE'} = 1;
            $wfngrp_namelist{'GNAME'}  = "GS^N positron-TGT-to-PCO^1";
            $wfngrp_namelist{'QNTAR'} = "1,$statesym,0";
            $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
            &add_GS_shells($r_par, \%wfngrp_namelist, $r_par->{'model'}->{'nelectrons'} );
            &add_positron_TGTPCO_shell_for_given_target_symmetry($r_par, \%wfngrp_namelist, $statesym);
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }

          $wfngrp = $wfngrp_template;
          if ($iposit == 0) {
            $wfngrp_namelist{'GNAME'}  = "GS^N-1 PCO^2";
            $wfngrp_namelist{'QNTAR'} = "-1,0,0";
            $wfngrp_namelist{'LSQUARE'} = 1;
            &add_GS_shells($r_par, \%wfngrp_namelist, $r_par->{'model'}->{'nelectrons'}-1 );
            &add_PCO_shells($r_par, \%wfngrp_namelist, 2);
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
          else {
            for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
              $found = 0;
              while ($found == 0) {
                for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
                  if ($group_table[$i]->[$j] == $statesym) {
                    $found = 1;
                    next if ($r_par->{'data'}->{'orbitals'}->{'PCO'}->[$j] <=0);
                    $wfngrp_namelist{'GNAME'} = "GS^N-1 PCO^1 positron-PCO^1";
                    $wfngrp_namelist{'QNTAR'} = "1,$i,0";
                    $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
                    $wfngrp_namelist{'LSQUARE'} = 1;
                    &add_GS_shells($r_par, \%wfngrp_namelist, $r_par->{'model'}->{'nelectrons'}-1 );
                    &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
                    &add_PCO_shell_for_given_target_symmetry($r_par, \%wfngrp_namelist, $j, 1, $iposit);
                    $wfngrp = $wfngrp_template;
                    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                    $$r_str .= $wfngrp;
                  }
                }
              }
            }
          }

        }else{
          $wfngrp = $wfngrp_template;
          if ($iposit == 0) {
            $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc PCO^1";
            $wfngrp_namelist{'QNTAR'} = "-1,0,0";
            $wfngrp_namelist{'LSQUARE'} = 1;
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
            &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
          else {
            for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
              $found = 0;
              while ($found == 0) {
                for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
                  if ($group_table[$i]->[$j] == $statesym) {
                    $found = 1;
                    next if ($r_par->{'data'}->{'orbitals'}->{'PCO'}->[$j] <=0);
		    next if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] <=0);
                    $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc positron-TGT-to-PCO^1";
                    $wfngrp_namelist{'QNTAR'} = "1,$i,0";
                    $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
                    $wfngrp_namelist{'LSQUARE'} = 1;
                    &add_frozen_shells($r_par, \%wfngrp_namelist);
                    &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
                    &add_positron_TGTPCO_shell_for_given_target_symmetry($r_par, \%wfngrp_namelist, $j);
                    $wfngrp = $wfngrp_template;
                    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                    $$r_str .= $wfngrp;
                  }
                }
              }
            }
          }

          $wfngrp = $wfngrp_template;
          if ($iposit == 0) {
            $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc-1 PCO^2";
            $wfngrp_namelist{'QNTAR'} = "-1,0,0";
            $wfngrp_namelist{'LSQUARE'} = 1;
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons-1);
            &add_PCO_shells($r_par, \%wfngrp_namelist, 2);
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
          else {
            for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
              $found = 0;
              while ($found == 0) {
                for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
                  if ($group_table[$i]->[$j] == $statesym) {
                    $found = 1;
                    next if ($r_par->{'data'}->{'orbitals'}->{'PCO'}->[$j] <=0);
		    next if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] <=0);
                    $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc-1 PCO^1 positron-PCO^1";
                    $wfngrp_namelist{'QNTAR'} = "1,$i,0";
                    $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
                    $wfngrp_namelist{'LSQUARE'} = 1;
                    &add_frozen_shells($r_par, \%wfngrp_namelist);
                    &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons-1);
                    &add_PCO_shells($r_par, \%wfngrp_namelist, 1);
                    &add_PCO_shell_for_given_target_symmetry($r_par, \%wfngrp_namelist, $j, 1, $iposit);
                    $wfngrp = $wfngrp_template;
                    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                    $$r_str .= $wfngrp;
                  }
                }
              }
            }
          }
        }
      }
      # then add (core+cas)^N+1
      if ($iposit == 0) {
        $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc+1";
        $wfngrp_namelist{'QNTAR'} = "-1,0,0";
        $wfngrp_namelist{'LSQUARE'} = 1;
        &add_frozen_shells($r_par, \%wfngrp_namelist);
        &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons + 1);
        $wfngrp = $wfngrp_template;
        &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
        $$r_str .= $wfngrp;
      }
      else {
        if ($r_par->{'model'}->{'use_PCO'} == 0 and $r_par->{'model'}->{'nvirtual'} == 0) {
        #This if is because: if using VOs, these configurations are included in core^Nc cas^N-Nc positron-TGT-to-virtual^1
        #                    if using PCOs, these configurations are included in core^Nc cas^N-Nc positron-TGT-to-PCO^1
          for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
            $found = 0;
            while ($found == 0) {
              for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
                if ($group_table[$i]->[$j] == $statesym) {
                  $found = 1;
                  next if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$j] <=0);
                  $wfngrp_namelist{'GNAME'} = "core^Nc cas^N-Nc positron^1";
                  $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
                  $wfngrp_namelist{'LSQUARE'} = 1;
                  $wfngrp_namelist{'QNTAR'} = "1,$i,0";
                  &add_frozen_shells($r_par, \%wfngrp_namelist);
                  &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
                  &add_positron_frozact_shells($r_par, \%wfngrp_namelist, $j);
                  $wfngrp = $wfngrp_template;
                  &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                  $$r_str .= $wfngrp;
                }
              }
            }
          }
        }
      }


      # and finally (core+cas)^N-h x virtual^1+h if there are any virtual orbitals
      if ($r_par->{'model'}->{'nvirtual'} > 0) {
        if ($r_par->{'model'}->{'model'} =~ /CAS-A/) {
          if ($iposit == 0 or ($iposit == -1 and $r_par->{'model'}->{'use_PCO'} == 0)) { #if doing positron scattering with PCOs, these configurations are included in core^Nc cas^N-Nc positron-TGT-to-PCO^1
            for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
              foreach my $targetstatespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
              if ($iposit == 0) {$orbs = 'virtual'}
              else {$orbs = 'target_used'}
                if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$spin_multiplicity{lc($targetstatespin)}.".$i"} &&                          # skip if target states of a given spin-symmetry are not used
                    $r_par->{'data'}->{'orbitals'}->{$orbs}->[get_orbitals_symmetry($i, $r_par->{'data'}->{'scattering'}->{'symmetry'})] > 0) { #   or if there are no virtual orbitals of required symmetry
                  my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i];
                  $wfngrp_namelist{'GNAME'} = "$targetstatespin $str_ir x virtual^1";
                  $wfngrp_namelist{'LSQUARE'} = 1;
                  &add_frozen_shells($r_par, \%wfngrp_namelist);
                  &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
                  if ($iposit == 0 ) {
                    $wfngrp_namelist{'QNTAR'} = "$spin_multiplicity{lc($targetstatespin)},$i,0";
                    &add_virtual_shell_for_given_target_symmetry($r_par, \%wfngrp_namelist, $i);
                  }
                  else {
                    $wfngrp_namelist{'QNTAR'} = "$spin_multiplicity{lc($targetstatespin)},$i,0";
                    $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
                    &add_positron_frozact_shells($r_par, \%wfngrp_namelist, $i);
                  }
                  $wfngrp = $wfngrp_template;
                  &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                  $$r_str .= $wfngrp;
                }
              }
            } # foreach spin state (multiplicity)
          } # for each IRs
        }
        else { # CAS or CAS-B or CAS-C
          if ($iposit == 0) {
            $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc virtual^1";
            $wfngrp_namelist{'LSQUARE'} = 1;
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
            &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
            $wfngrp = $wfngrp_template;
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
          else {
            if ($r_par->{'model'}->{'use_PCO'} == 0) { #if using PCOs then these configurations are included in core^Nc cas^N-Nc positron-TGT-to-PCO^1
              for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
                $found = 0;
                while ($found == 0) {
                  for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
                    if ($group_table[$i]->[$j] == $statesym) {
                      $found = 1;
                      next if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$j] <=0);
                      $wfngrp_namelist{'QNTAR'} = "1,$i,0";
                      $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc positron-TGT-to-virtual^1";
                      $wfngrp_namelist{'LSQUARE'} = 1;
                      $wfngrp_namelist{'L2_POSITRON_FLAG'} = '';
                      &add_frozen_shells($r_par, \%wfngrp_namelist);
                      &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons);
                      &add_positron_frozact_shells($r_par, \%wfngrp_namelist, $j);
                      $wfngrp = $wfngrp_template;
                      &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                      $$r_str .= $wfngrp;
                    }
                  }
                }
              }
            }
          }
        }

        # for CAS-C add (core+cas)^N-1 x virtual^2 yet
        if ($r_par->{'model'}->{'model'} =~ /CAS-C/) {
          if ($iposit == 0) {
            $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc-1 virtual^2";
            $wfngrp_namelist{'LSQUARE'} = 1;
            &add_frozen_shells($r_par, \%wfngrp_namelist);
            &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
            &add_virtual_shells($r_par, \%wfngrp_namelist, 2);
            $wfngrp = $wfngrp_template;
            &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
            $$r_str .= $wfngrp;
          }
          else {
            for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
              $found = 0;
              while ($found == 0) {
                for (my $j = 0; $j < $r_par->{'data'}->{'nir'}; $j++) {
                  if ($group_table[$i]->[$j] == $statesym) {
                    $found = 1;
                    next if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$j] <=0);
                    $wfngrp_namelist{'QNTAR'} = "1,$i,0";
                    $wfngrp_namelist{'GNAME'}  = "core^Nc cas^N-Nc-1 virtual^2";
                    $wfngrp_namelist{'LSQUARE'} = 1;
                    &add_frozen_shells($r_par, \%wfngrp_namelist);
                    &add_active_shells($r_par, \%wfngrp_namelist, $nactive_electrons - 1);
                    &add_virtual_shells($r_par, \%wfngrp_namelist, 1);
                    &add_positron_frozact_shells($r_par, \%wfngrp_namelist, $j);
                    $wfngrp = $wfngrp_template;
                    &replace_all_in_template(\$wfngrp, \%wfngrp_namelist);
                    $$r_str .= $wfngrp;
                  }
                }
              }
            }
          }
        }
      }

    }

  } # end of scattering calculation setting

  return 1;
}

sub add_frozen_shells {
  my ($r_par, $r_wfngrp_namelist) = @_;
  my ($nelecp, $nshlp) = (0, 0);

  # clean from previous &wfngrp
  $r_wfngrp_namelist->{'NDPROD'} = 0;
  $r_wfngrp_namelist->{'NELECP'} = "";
  $r_wfngrp_namelist->{'NSHLP'}  = "";
  $r_wfngrp_namelist->{'PQN'}    = "";
  $r_wfngrp_namelist->{'MSHL'}   = "";

  if ($r_par->{'model'}->{'nfrozen'} > 0) {
    $r_wfngrp_namelist->{'NDPROD'}++;
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      if ($r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i] > 0) {
        $nelecp += 2 * $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i];
        $nshlp++;
        $r_wfngrp_namelist->{'PQN'}  .= "0,1,$r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i], ";
        $r_wfngrp_namelist->{'MSHL'} .= "  $i,   ";
      }
    }
    $r_wfngrp_namelist->{'NELECP'} .= "$nelecp,";
    $r_wfngrp_namelist->{'NSHLP'}  .= "$nshlp,";
  }

  return 1;
}

sub add_active_shells {
  my ($r_par, $r_wfngrp_namelist, $nel) = @_;

  if ($r_par->{'model'}->{'nactive'} > 0) {
    $r_wfngrp_namelist->{'NDPROD'}++;
    my $nshlp = 0;
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      if ($r_par->{'data'}->{'orbitals'}->{'active'}->[$i] > 0) {
        $nshlp++;
        my $nfrozen = $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i];
        $r_wfngrp_namelist->{'PQN'}  .= "0,".($nfrozen + 1).",".($nfrozen + $r_par->{'data'}->{'orbitals'}->{'active'}->[$i]).", ";
        $r_wfngrp_namelist->{'MSHL'} .= "  $i,   ";
      }
    }
    $r_wfngrp_namelist->{'NELECP'} .= "$nel,";
    $r_wfngrp_namelist->{'NSHLP'}  .= "$nshlp,";
  }

  return 1;
}

sub add_positron_frozact_shells {
  # THIS IS FOR POSITRON SCATTERING - POSITRON OCCUPYING TARGET TO ACTIVE ORBITALS (AND VIRTUAL IF PRESENT)
  my ($r_par, $r_wfngrp_namelist, $sym) = @_;

  if ($r_par->{'data'}->{'orbitals'}->{'target_used'}->[$sym] > 0  and $r_par->{'model'}->{'positron_flag'} == 1) {
    $r_wfngrp_namelist->{'NDPROD'}++;
    my $nshlp = 0;
    $nshlp++;
    my $nstart = $r_par->{'data'}->{'orbitals'}->{'all_used'}->[$sym];
    my $nfrozen = $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$sym];
    my $nactive = $r_par->{'data'}->{'orbitals'}->{'active'}->[$sym];
    my $nvirtual = $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$sym];
    $r_wfngrp_namelist->{'PQN'}  .= "0,".($nstart + 1).",".($nstart + $nfrozen + $nactive + $nvirtual).", ";
    $r_wfngrp_namelist->{'MSHL'} .= "  $sym,   ";

    $r_wfngrp_namelist->{'NELECP'} .= "1,";
    $r_wfngrp_namelist->{'NSHLP'}  .= "1,";
  }

  return 1;
}

sub add_PCO_shells {
  my ($r_par, $r_wfngrp_namelist,$nelec) = @_;

  my $nshlp = 0;
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    my $nPCO = $r_par->{'data'}->{'orbitals'}->{'PCO'}->[$i];
    if ($nPCO > 0) {
      my $nTGT = $r_par->{'data'}->{'orbitals'}->{'TGT'}->[$i];
      $r_wfngrp_namelist->{'PQN'}  .= "0,".($nTGT + 1).",".($nTGT+$nPCO).", ";
      $r_wfngrp_namelist->{'MSHL'} .= "  $i,   ";
      $nshlp++;
    }
  }
  if ($nshlp > 0) {
    $r_wfngrp_namelist->{'NDPROD'}++;
    $r_wfngrp_namelist->{'NSHLP'}  .= "$nshlp,";
    $r_wfngrp_namelist->{'NELECP'} .= "$nelec,";
  }

  return 1;
}

sub add_PCO_shell_for_given_target_symmetry {
  my ($r_par, $r_wfngrp_namelist, $target_state_sym, $nelec, $iposit) = @_;

  my $PCO_sym = get_orbitals_symmetry($target_state_sym, $r_par->{'data'}->{'scattering'}->{'symmetry'});
  my $nPCO = $r_par->{'data'}->{'orbitals'}->{'PCO'}->[$PCO_sym];
  if ( $nPCO > 0) {
    my $nTGT = $r_par->{'data'}->{'orbitals'}->{'TGT'}->[$PCO_sym];
    my $nelectron = 0;
    if ($iposit != 0){$nelectron = $r_par->{'data'}->{'orbitals'}->{'all_used'}->[$PCO_sym];}
    $r_wfngrp_namelist->{'PQN'}  .= "0,".($nelectron + $nTGT + 1).",".($nelectron + $nTGT+$nPCO).", ";
    $r_wfngrp_namelist->{'MSHL'} .= "  $PCO_sym,   ";
    $r_wfngrp_namelist->{'NELECP'} .= "$nelec,";
    $r_wfngrp_namelist->{'NSHLP'}  .= "1,";
    $r_wfngrp_namelist->{'NDPROD'}++;
  }

  return 1;
}

sub add_positron_TGTPCO_shell_for_given_target_symmetry {
  # POSITRON OCCUPIES ORBITALS FROM TARGET TO PCO
  my ($r_par, $r_wfngrp_namelist, $target_state_sym) = @_;

  my $PCO_sym = get_orbitals_symmetry($target_state_sym, $r_par->{'data'}->{'scattering'}->{'symmetry'});
  my $nPCO = $r_par->{'data'}->{'orbitals'}->{'PCO'}->[$PCO_sym];
  if ( $nPCO > 0) {
    my $nTGT = $r_par->{'data'}->{'orbitals'}->{'TGT'}->[$PCO_sym];
    my $nelectron = $r_par->{'data'}->{'orbitals'}->{'all_used'}->[$PCO_sym];
    $r_wfngrp_namelist->{'PQN'}  .= "0,".($nelectron + 1).",".($nelectron + $nTGT+$nPCO).", ";
    $r_wfngrp_namelist->{'MSHL'} .= "  $PCO_sym,   ";
    $r_wfngrp_namelist->{'NELECP'} .= "1,";
    $r_wfngrp_namelist->{'NSHLP'}  .= "1,";
    $r_wfngrp_namelist->{'NDPROD'}++;
  }

  return 1;
}


sub add_reduced_PCO_shells {
  my ($r_par, $r_wfngrp_namelist, $nelec) = @_;

  my $nshlp = 0;
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    my $nPCO = $r_par->{'data'}->{'orbitals'}->{'PCO'}->[$i];
    if ($nPCO > 0 ) {
      my $nTGT = $r_par->{'data'}->{'orbitals'}->{'TGT'}->[$i];
      $r_wfngrp_namelist->{'PQN'}  .= "0,".($nTGT + 1).",".($nTGT+$nPCO).", ";
      $r_wfngrp_namelist->{'MSHL'} .= "  $i,   ";
      $nshlp++;
    }
  }
  if ($nshlp > 0) {
    $r_wfngrp_namelist->{'NDPROD'}++;
    $r_wfngrp_namelist->{'NSHLP'}  .= "$nshlp,";
    $r_wfngrp_namelist->{'NELECP'} .= "$nelec,";
  }

  return 1;
}

sub add_GS_shells {
  my ($r_par, $r_wfngrp_namelist, $nelec) = @_;

  # clean from previous &wfngrp
  $r_wfngrp_namelist->{'NDPROD'} = 0;
  $r_wfngrp_namelist->{'NELECP'} = "";
  $r_wfngrp_namelist->{'NSHLP'}  = "";
  $r_wfngrp_namelist->{'PQN'}    = "";
  $r_wfngrp_namelist->{'MSHL'}   = "";

  my $nshlp = 0;
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    my $nGS = $r_par->{'data'}->{'orbitals'}->{'GS'}->[$i];
    if ($nGS > 0 ) {
      $r_wfngrp_namelist->{'PQN'}  .= "0,1,$nGS, ";
      $r_wfngrp_namelist->{'MSHL'} .= "  $i,   ";
      $nshlp++;
    }
  }
  if ($nshlp > 0) {
    $r_wfngrp_namelist->{'NDPROD'}++;
    $r_wfngrp_namelist->{'NSHLP'}  .= "$nshlp,";
    $r_wfngrp_namelist->{'NELECP'} .= "$nelec,";
  }

  return 1;
}

sub add_continuum_shells {
  my ($r_par, $r_wfngrp_namelist, $target_state_sym, $iposit) = @_;

  $r_wfngrp_namelist->{'NDPROD'}++;
  my $cont_sym = get_orbitals_symmetry($target_state_sym, $r_par->{'data'}->{'scattering'}->{'symmetry'});
  my $ntarget = $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$cont_sym];
  my $nelectron = 0;
  if ($iposit != 0) {$nelectron = ($r_par->{'data'}->{'orbitals'}->{'all_used'}->[$cont_sym]);}
  $r_wfngrp_namelist->{'PQN'}  .= "0,".($nelectron + $ntarget + 1).",".($nelectron + $ntarget + 2).", ";
  $r_wfngrp_namelist->{'MSHL'} .= "  $cont_sym,   ";
  $r_wfngrp_namelist->{'NELECP'} .= "1,";
  $r_wfngrp_namelist->{'NSHLP'}  .= "1,";

  return 1;
}

sub add_virtual_shells {
  my ($r_par, $r_wfngrp_namelist, $nel) = @_;

  $r_wfngrp_namelist->{'NDPROD'}++;
  my $nshlp = 0;
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    if ($r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i] > 0) {
      $nshlp++;
      my $ntarget = $r_par->{'data'}->{'orbitals'}->{'target'}->[$i];
      $r_wfngrp_namelist->{'PQN'}  .= "0,".($ntarget + 1).",".($ntarget + $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i]).", ";
      $r_wfngrp_namelist->{'MSHL'} .= "  $i,   ";
    }
  }
  $r_wfngrp_namelist->{'NELECP'} .= "$nel,";
  $r_wfngrp_namelist->{'NSHLP'}  .= "$nshlp,";

  return 1;
}

sub add_virtual_shell_for_given_target_symmetry {
  my ($r_par, $r_wfngrp_namelist, $target_state_sym) = @_;

  $r_wfngrp_namelist->{'NDPROD'}++;
  my $virt_sym = get_orbitals_symmetry($target_state_sym, $r_par->{'data'}->{'scattering'}->{'symmetry'});
  my $ntarget = $r_par->{'data'}->{'orbitals'}->{'target'}->[$virt_sym];
  if ($r_par->{'data'}->{'orbitals'}->{'virtual'}->[$virt_sym] > 1) {
    $r_wfngrp_namelist->{'PQN'}  .= "0,".($ntarget + 1).",".($ntarget + 2).", ";
  }
  else {
    $r_wfngrp_namelist->{'PQN'}  .= "0,".($ntarget + 1).",".($ntarget + 1).", ";
  }
  $r_wfngrp_namelist->{'MSHL'} .= "  $virt_sym,   ";
  $r_wfngrp_namelist->{'NELECP'} .= "1,";
  $r_wfngrp_namelist->{'NSHLP'}  .= "1,";

  return 1;
}

sub add_active_shell_for_CHF {
  my ($r_par, $r_wfngrp_namelist, $sym, $j, $nel) = @_;

  my $nfrozen = $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$sym];

  $r_wfngrp_namelist->{'NDPROD'}++;
  $r_wfngrp_namelist->{'PQN'}  .= "0,".($nfrozen + $j).",".($nfrozen + $j).", ";
  $r_wfngrp_namelist->{'MSHL'} .= "  $sym,   ";
  $r_wfngrp_namelist->{'NELECP'} .= "$nel,";
  $r_wfngrp_namelist->{'NSHLP'}  .= "1,";

  return 1;
}

sub add_virtual_shell_for_CHF {
  my ($r_par, $r_wfngrp_namelist, $vsym, $virt_no) = @_;

  $r_wfngrp_namelist->{'NDPROD'}++;
  my $ntarget = $r_par->{'data'}->{'orbitals'}->{'target'}->[$vsym];
  $r_wfngrp_namelist->{'PQN'}  .= "0,".($ntarget + $virt_no).",".($ntarget + $virt_no).", ";
  $r_wfngrp_namelist->{'MSHL'} .= "  $vsym,   ";
  $r_wfngrp_namelist->{'NELECP'} .= "1,";
  $r_wfngrp_namelist->{'NSHLP'}  .= "1,";

  return 1;
}

# Auxiliary subroutine to determine symmetry of orbitals
# from group multiplication table in such a way that
# target state symmetry  x  orbitals symmetry = scattering state symmetry
sub get_orbitals_symmetry {
  my ($target_sym, $scat_sym) = @_;
  my $orb_sym = 0;
  IR: for (my $ir = 0; $ir < scalar @group_table; $ir++) {
    if ($group_table[$target_sym]->[$ir] == $scat_sym) { $orb_sym = $ir; last IR; }
  }
  return $orb_sym;
}

# subroutine which returns reference orbitals of a CSF with a given total spin and symmetry
sub get_reference_orbitals {
  my ($r_par, $nelectrons, $state_sym, $state_spin) = @_;

  my @frozen = @{$r_par->{'data'}->{'orbitals'}->{'frozen'}};
  my @reference = @{$r_par->{'data'}->{'orbitals'}->{'reference'}};

  my @el_pos  = (); # current positions of all electrons
  my @orb_sym = (); # array of symmetries of orbitals

  # first fill up frozen orbitals
  my $nel = $nelectrons;  # number of remaining electrons to put in orbitals
  my $norb = 0;           # number of all orbitals
  for (my $i = 0; $i < scalar @frozen; $i++) { # loop over all IRs
    for (my $j = 0; $j < @frozen[$i]; $j++) { # loop over all orbitals in IR $i
      if    ($nel >  1) { $el_pos[$nelectrons - $nel] = $norb + 1; $el_pos[$nelectrons - $nel + 1] = -$norb - 1; }
      elsif ($nel == 1) { $el_pos[$nelectrons - $nel] = $norb + 1; }
      $orb_sym[$norb] = $i;
      $nel -= 2;
      $norb++;
    }
  }

  # then fill remaining reference orbitals
  for (my $i = 0; $i < scalar @reference; $i++) { # loop over all IRs
    for (my $j = $frozen[$i]; $j < @reference[$i]; $j++) { # loop over all orbitals in IR $i
      if    ($nel >  1) { $el_pos[$nelectrons - $nel] = $norb + 1; $el_pos[$nelectrons - $nel + 1] = -$norb - 1; }
      elsif ($nel == 1) { $el_pos[$nelectrons - $nel] = $norb + 1; }
      $orb_sym[$norb] = $i;
      $nel -= 2;
      $norb++;
    }
  }

  my $found = 0;
  while ($found == 0) {
    # first check the configuration, if symmmetry and spin are correct then stop
#    &print_info("@el_pos\n", $r_par);
    my ($sym, $spin) = &get_csf_symmetry_and_spin(\@orb_sym, \@el_pos);
#    &print_info("Total sym and spin is $sym and $spin\n", $r_par);
    if ($sym == $state_sym && $spin == $state_spin) {
      $found = 1;
    }
    else {
      # find electron to move
      my $move_el = $nelectrons;
      my $n = scalar @orb_sym;
      if ($el_pos[$move_el - 1] == -$n) {                              # if the last electron is in the last orbital with spin down
        while (abs($el_pos[$move_el - 1]) == $n && $move_el > 0) {     # then look for an electron (in reverse order) untill there is one which can be moved
          if    ($nelectrons % 2 == 0 && $move_el % 2 == 1) { $n--; }  # along the way we have to change the highest possible occupied orbital
          elsif ($nelectrons % 2 == 1 && $move_el % 2 == 0) { $n--; }  # bacause there can be only two electrons in each orbital
          $move_el--;
        }
      }
      if ($move_el == 0) { # if no electron can be moved there is no configuration for a given state spin and symmetry
        $found = -1;
      }
      else {               # otherwise move the electron
        $move_el--;        # shift because the first element has an index 0
        if ($el_pos[$move_el] > 0) {                            # if the electron has spin up then
          $el_pos[$move_el] = -$el_pos[$move_el];               # change spin
        }
        else {                                                  # if the electron has spin down then
          $el_pos[$move_el] = -$el_pos[$move_el] + 1;           # move it
        }
        for (my $i = $move_el + 1; $i < $nelectrons; $i++ ) {   # the rest of electrons place to the next available orbitals
          if ($el_pos[$i - 1] < 0) { $el_pos[$i] = -$el_pos[$i - 1] + 1; }
          else                     { $el_pos[$i] = -$el_pos[$i - 1]; }
        }
      } # if ($move_el == 0)
    } # if ($sym == $state_sym && $spin == $state_spin)
  } # while ($found == 0)
  if ($found == -1) {
    &print_info("\n  Warning: no configuration found for $nelectrons electrons\n", $r_par);
    &print_info("             for a state with (sym, spin) = ($state_sym, $state_spin) !!!\n", $r_par);
    &print_info("             Orbitals used: @reference\n\n", $r_par);
    return (0, "");
  }
  else {
    return &get_reforb_quintets_for_csf(\@orb_sym, \@el_pos);
  }
}

# subroutine which returns total spin and symmetry from @orb_sym and @el_pos
sub get_csf_symmetry_and_spin {
  my ($r_orb_sym, $r_el_pos) = @_;
  my ($sym, $spin) = (0, 0);

  for (my $i = 0; $i < scalar @$r_el_pos; $i++) { # loop over all electrons
    if ($r_el_pos->[$i] > 0) { $spin++; }         # add spin up
    else                     { $spin--; }         # add spin down
    # multiply the current total state symmetry by symmetry of the orbital in which the electron is
    $sym = $group_table[$sym]->[$r_orb_sym->[abs($r_el_pos->[$i]) - 1]];
  }
  if ($spin < 0) { $spin = -$spin; }
  $spin++; # to have 1 for singlet etc.

  return ($sym, $spin);
}

# subroutine which returns nrefo and quintets in reforb
# from arrays of orbitals symmetries @orb_sym and positions of electrons @el_pos
sub get_reforb_quintets_for_csf {
  my ($r_orb_sym, $r_el_pos) = @_;
  my $nrefo = 0;
  my @reforb = ();

  my $nelectrons = scalar @$r_el_pos;
  my $nel = 1;                                         # number of electrons for a quintet
  my @start_orb = (1,1,1,1,1,1,1,1);                   # starting orbital to be used in a quintet
  my ($prev_orb, $prev_orb_sym);                       # number and symmetry of the previous orbital
  my $curr_orb = $r_el_pos->[0];                       # number of the current orbital
  my $curr_orb_sym = $r_orb_sym->[abs($curr_orb) - 1]; # symmetry of the current orbital
  for (my $i = 1; $i < $nelectrons; $i++) {            # loop over all electrons but the first

    $prev_orb = $curr_orb;
    $prev_orb_sym = $curr_orb_sym;
    $curr_orb = $r_el_pos->[$i];
    $curr_orb_sym = $r_orb_sym->[abs($curr_orb) - 1];

    if ($curr_orb_sym == $prev_orb_sym                 # if the same symmetry of the current orbital and previous one
       && ($curr_orb == -$prev_orb                     # and the electron is at the same orbital as the previous one
           || (abs($curr_orb) == abs($prev_orb) + 1    # or it is in the next orbital
               && $nel % 2 == 0) )) {                  #    and at the same time there is even number of electrons in the quintet
      $nel++;                                          # then just add one electron to the quintet
    }
    else {                                             # otherwise we have to save the quintet and start a new one
      $nrefo++;
      push @reforb, $prev_orb_sym, @start_orb[$prev_orb_sym], $nel, ($nel % 2 == 1 && $prev_orb < 0) ? 1 : 0, 0;
      @start_orb[$prev_orb_sym] += int(($nel+1)/2);
      $nel = 1;
    }
  }
  $nrefo++;
  push @reforb, $curr_orb_sym, @start_orb[$curr_orb_sym], $nel, ($nel % 2 == 1 && $curr_orb < 0) ? 1 : 0, 0;
  return ($nrefo, @reforb);
}


# ------------------------ scatci ---------------------------

sub make_scatci_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};
  my $multiplicity = $spin_multiplicity{lc($statespin)};

  my %input_namelist = (
    'MOLECULE',         $r_par->{'model'}->{'molecule'},
    'MEGUL',            "${lucsf_base}${multiplicity}${statesym}",
    'NFTE',             "${luham_base}${multiplicity}${statesym}",
    'MPI_SCATCI',       (($r_par->{'run'}->{'mpi_scatci'} eq "" and $r_par->{'run'}->{"parallel_symm"} >= 1) ? "!" : ""),
    'SPIN',             $statespin,
    'SYMMETRY',         $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym],
    'IPOSIT',           -$r_par->{'model'}->{'positron_flag'},
  );

  # For each target csf we calculate given number of eigenvalues
  # from which later in make_denprop_input select states with lowest energies
  if ($task eq "target") {

    # First update auxiliary variables used later for denprop and scattering.scatci inputs
    push @{$r_par->{'data'}->{'target'}->{'spinsym_order'}}, "$multiplicity.$statesym";
    my $nciset = scalar @{$r_par->{'data'}->{'target'}->{'spinsym_order'}};

    if ($r_par->{'run'}->{'parallel_symm'} eq 0) {
      # force dense diagonalization when using mpi-scatci
      $input_namelist{'IGHT'} = '';
      $input_namelist{'IGH'} = 1;
    }
    elsif ($r_par->{'run'}->{'ight'} eq 2) {
      $input_namelist{'IGHT'} = '!';
      $input_namelist{'IGH'} = 1;
    }
    else {
      $input_namelist{'IGHT'} = '';
      $input_namelist{'IGH'} = $r_par->{'run'}->{'ight'};
    }

    $input_namelist{'IPOSIT'} = 0;
    if ($r_par->{'model'}->{'model'} eq "CHF-A") {
      my $megul = sprintf('%03s', $r_par->{'data'}->{'CHF'}->{'current_state'});
      $input_namelist{'MEGUL'} = "${lucsf_base}${megul}";
      $input_namelist{'NFTE'} = "${luham_base}${megul}";
    }
    $input_namelist{'NSTAT'} = $r_par->{'model'}->{'ntarget_states'}->{$statespin}->[$statesym];
    $input_namelist{'LUCITGT'} = $lucitgt;
    $input_namelist{'NCISET'} = $nciset;

    # This option is only used by MPI_SCATCI and is commented out from the input otherwise.
    $input_namelist{'VECSTORE'} = 5;
    if ($r_par->{'run'}->{'parallel_symm'} <= 0) {
      $input_namelist{'MPI_SCATCI'} = '';
      $input_namelist{'VECSTORE'} = 7;
    }

  }
  # For scattering calculation we need all eigenvalues
  else {

    my $nciset = 1;
    push @{$r_par->{'data'}->{'scattering'}->{'spinsym_order'}}, "$multiplicity.$statesym";
    if ($r_par->{'run'}->{'rmt_interface'} or $r_par->{'run'}->{'parallel_symm'} <= 0) {
      $nciset = scalar @{$r_par->{'data'}->{'scattering'}->{'spinsym_order'}};
    }

    my @mcont = map($group_table[$statesym]->[$_], @{$r_par->{'data'}->{'target'}->{'mcont'}});
    my $notgt = "";
    for (my $i = 0; $i < scalar @mcont; $i++) {
      if ($r_par->{'data'}->{'orbitals'}->{'virt_cont_used'}){
        # For contraction of virtuals with continuum
        $notgt .= "$r_par->{'data'}->{'orbitals'}->{'virt_cont_used'}->[$mcont[$i]],";
      } else {
        $notgt .= "$r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$mcont[$i]],";
      }
    }

    if ($r_par->{'run'}->{'parallel_symm'} eq 0) {
      # force dense diagonalization when using mpi-scatci
      $input_namelist{'IGHS'} = '';
      $input_namelist{'IGH'} = 1;
    }
    if ($r_par->{'run'}->{'ighs'} eq 2) {
      $input_namelist{'IGHS'} = '!';
      $input_namelist{'IGH'} = 1;
    }
    else {
      $input_namelist{'IGHS'} = '';
      $input_namelist{'IGH'} = $r_par->{'run'}->{'ighs'};
    }

    $input_namelist{'NTGSYM'} = $r_par->{'data'}->{'target'}->{'ntgt'};
    $input_namelist{'MCONT'}  = join(",", @mcont);
    $input_namelist{'NOTGT'}  = $notgt;
    $input_namelist{'NUMTGT'} = join(",", @{$r_par->{'data'}->{'target'}->{'ntgtl'}});
    $input_namelist{'ICIDG'} = ($r_par->{'run'}->{'parallel_diag'} == 1 ? 0 : 1);
    $input_namelist{'LUCI'} = $luci;
    $input_namelist{'NCISET'} = $nciset;

    # We need the full CI vectors only for photoionization and RMT interface (= 5).
    # For scattering only the continuum part is needed (= 1).
    # This option is only used by MPI_SCATCI and is commented out from the input otherwise.
    $input_namelist{'VECSTORE'} = 1;
    if ($r_par->{'run'}->{'photoionization'} == 1 or $r_par->{'run'}->{'rmt_interface'} == 1) {
      $input_namelist{'VECSTORE'} = 5;
    }
    if ($r_par->{'run'}->{'parallel_symm'} <= 0) {
      $input_namelist{'MPI_SCATCI'} = '';
      $input_namelist{'VECSTORE'} = 7;
    }

    # Ensuring phase consistency for transition dipoles and outer-region potential coefficients
    #   NFTG  ... unit with target eigenstates
    #   NTGTF ... for each used target eigenstate (in congen order): dataset in NFTG where it is stored
    #   NTGTS ... for each used target eigenstate (in congen order): position of the eigenstate in the dataset NTGTF
    my $ntgt  = $r_par->{'data'}->{'target'}->{'ntgt'};      # number of spin-symmetries used
    my @ntgtl = @{$r_par->{'data'}->{'target'}->{'ntgtl'}};  # number of target states per used spin-symmetry
    my @ntgtf = @{$r_par->{'data'}->{'target'}->{'ntgtf'}};  # dataset index where they are stored
    $input_namelist{'NFTG'} = $lucitgt;
    $input_namelist{'NTGTF'} = join(",", map { ($ntgtf[$_ - 1]) x $ntgtl[$_ - 1] } 1..$ntgt);
    $input_namelist{'NTGTS'} = join(",", map { 1..$_ } @ntgtl);

  }
  &replace_all_in_template($r_str, \%input_namelist);
  return 1;
}

# ----------------------- denprop ---------------------------

sub make_denprop_input {
  my ($r_par, $r_str) = @_;
  my $r_states = $r_par->{'data'}->{'target'}->{'states'};

  # According to number of target calculations (different spin-symmetry states) set NTGTF as "1,2,3,..."
  my $ntgtf = "";
  for (my $i = 1; $i <= $r_par->{'data'}->{'target'}->{'ntgt'}; $i++) {
    $ntgtf .= "$i,";
  }

  # Get number of target states used in each spin-symmetry from the number of states
  # as required in settings -> ntarget_states_used (only states with lowest energies are used)

  # First order states by energy and check whether there is enough states as required if not use only calculated ones
  my @ordered_states = sort { $r_states->{$a} <=> $r_states->{$b} } keys %{$r_states};
  $r_par->{'data'}->{'target'}->{'ordered_states'} = \@ordered_states;
  $r_par->{'data'}->{'target'}->{'ground_state'} = $ordered_states[0];
  my $r_ordered = $r_par->{'data'}->{'target'}->{'ordered_states'};
  my $nstates = scalar @$r_ordered;
  if ($nstates < $r_par->{'model'}->{'ntarget_states_used'}) {
    &print_info("  Warning: There is less target states ($nstates) in scatci output(s) than required ($r_par->{'model'}->{'ntarget_states_used'}) !!!\n", $r_par);
    &print_info("           Increase numbers of states in settings->ntarget_states.\n", $r_par);
    &print_info("           Only those states which were found will be used !\n", $r_par);
    $r_par->{'model'}->{'ntarget_states_used'} = $nstates;
  }

  if ($r_par->{'model'}->{'model'} eq "CHF-A") {
    my @states = @{$r_par->{'data'}->{'CHF'}->{'states'}};
    my $no_states = scalar @states;
    $r_par->{'data'}->{'target'}->{'ntgt'} = $no_states;
    # Now get number of states of given spin and symmetry within the lowest-lying states
    $r_par->{'data'}->{'target'}->{'used_tgt_states'} = {};
    for (my $i = 0; $i < scalar @{$r_ordered}; $i++) {
      my ($spin, $sym, $number) = split(/\./, $r_ordered->[$i]);
      $number = $number + 0; #Make sure it is a number and NOT a string
      my $current = $r_par->{'data'}->{'target'}->{'used_tgt_states'}->{"$spin.$sym"};
      if ( $number > $current ){
        $r_par->{'data'}->{'target'}->{'used_tgt_states'}->{"$spin.$sym"} = $number;
      }
    }
  }
  else {
    # Now get number of states of given spin and symmetry within the lowest-lying states
    $r_par->{'data'}->{'target'}->{'used_tgt_states'} = {};
    for (my $i = 0; $i < $r_par->{'model'}->{'ntarget_states_used'}; $i++) {
      my ($spin, $sym, $number) = split(/\./, $r_ordered->[$i]);
      $number = $number + 0; #Make sure it is a number and NOT a string
      my $current = $r_par->{'data'}->{'target'}->{'used_tgt_states'}->{"$spin.$sym"};
      if ($number > $current) {
        $r_par->{'data'}->{'target'}->{'used_tgt_states'}->{"$spin.$sym"} = $number;
      }
    }
    $r_par->{'data'}->{'target'}->{'ntgt'} = scalar (keys %{$r_par->{'data'}->{'target'}->{'used_tgt_states'}});
  }

  # Finally set variables for scatci input
  # and also auxiliary arrays 'mcont' and 'ntgtf' used later in scattering.scatci inputs
  my $nftsor = "";
  $r_par->{'data'}->{'target'}->{'ntgtl'} = [];
  $r_par->{'data'}->{'target'}->{'ntgtf'} = [];
  $r_par->{'data'}->{'target'}->{'mcont'} = [];
  if ($r_par->{'model'}->{'model'} eq "CHF-A") {
    my @states = @{$r_par->{'data'}->{'CHF'}->{'states'}};
    for $state (@states){
      my $id = $state->{'id'};
      my $sym = $state->{'sym'};
      my $ncsfs = $state->{'no_csfs'};
      my $state_id = sprintf('%03s', $id);
      $nftsor .= "$lucsf_base$state_id,";
      push @{$r_par->{'data'}->{'target'}->{'ntgtf'}}, $id;
      push @{$r_par->{'data'}->{'target'}->{'ntgtl'}}, $ncsfs;
      push @{$r_par->{'data'}->{'target'}->{'mcont'}}, $sym;
    }
  }
  else {
    my $r_spinsym_order = $r_par->{'data'}->{'target'}->{'spinsym_order'};
    for (my $i = 0; $i < scalar @{$r_spinsym_order}; $i++) {
      if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$r_spinsym_order->[$i]}) {
        my ($spin, $sym) = split(/\./, $r_spinsym_order->[$i]);
        $nftsor .= "$lucsf_base$spin$sym,";
        push @{$r_par->{'data'}->{'target'}->{'ntgtf'}}, $i + 1;
        push @{$r_par->{'data'}->{'target'}->{'ntgtl'}}, $r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$r_spinsym_order->[$i]};
        push @{$r_par->{'data'}->{'target'}->{'mcont'}}, $sym;
      }
    }
  }

  if ($r_par->{'run'}->{'photoionization'} == 1 or $r_par->{'run'}->{'rmt_interface'} == 1) {
    $electronic_properties_only = 1;  # exclude nuclear properties for electron-laser interaction
  } else {
    $electronic_properties_only = 0;  # include nuclear properties for electron-molecule interaction
  }

  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "IPOL",     $r_par->{'model'}->{'max_multipole'});
  &replace_in_template($r_str, "NTGT",     $r_par->{'data'}->{'target'}->{'ntgt'});
  &replace_in_template($r_str, "NFTSOR",   $nftsor);
  &replace_in_template($r_str, "NFTINT",   $luamp);
  &replace_in_template($r_str, "LUCITGT",  $lucitgt);
  &replace_in_template($r_str, "LUPROP",   $luprop);
  &replace_in_template($r_str, "NTGTF",    join(",", @{$r_par->{'data'}->{'target'}->{'ntgtf'}}));
  &replace_in_template($r_str, "NTGTL",    join(",", @{$r_par->{'data'}->{'target'}->{'ntgtl'}}));
  &replace_in_template($r_str, "ISW",      $electronic_properties_only);
  return 1;
}

# ----------------------- mpi-scatci ---------------------------

sub make_mpi_scatci_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $state_selection = $task eq "target" ? "ntarget_states" : "scattering_states";

  # Read scatci input templates
  my $scatci_template = "";
  if(!read_file("$r_par->{'dirs'}->{'templates'}${bs}$task.scatci.inp", \$scatci_template)) {
    print_info("Warning: no template file for $task.scatci.inp !\n", $r_par);
  }
  my $interface_template = "";
  if(!read_file("$r_par->{'dirs'}->{'templates'}${bs}$task.scatci.interface.inp", \$interface_template)) {
    print_info("Warning: no template file for $task.scatci.interface.inp !\n", $r_par);
  }

  # Generate inputs for the diagonalizations
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{$state_selection}}) {
      $r_par->{'data'}->{$task}->{'spin'} = $statespin;
      # skip spin-symmetries which are not required
      if ($r_par->{'model'}->{$state_selection}->{$statespin}->[$i] > 0) {
        $r_par->{'data'}->{$task}->{'symmetry'} = $i;
        my $str = $scatci_template;
        make_scatci_input($r_par, \$str);
        $$r_str .= $str;
      }
    }
  }

  # Generate inputs for the interface
  replace_in_template(\$interface_template, "WRITE_RMT", $r_par->{'run'}->{'rmt_interface'} ? ".true." : ".false.");
  replace_in_template(\$interface_template, "RMATR",     $r_par->{'model'}->{'rmatrix_radius'});
  replace_in_template(\$interface_template, "LUPROPW",   $luprop);
  replace_in_template(\$interface_template, "ICFORM",    $icform);
  replace_in_template(\$interface_template, "IRFORM",    $irform);
  if (task == "target") {
    if ($r_par->{'run'}->{'photoionization'} or $r_par->{'run'}->{'rmt_interface'}) {
      replace_in_template(\$interface_template, "ISW", 1);
    } else {
      replace_in_template(\$interface_template, "ISW", 0);
    }
  } else {
    replace_in_template(\$interface_template, "ISW", 1);
  }

  $$r_str .= $interface_template;

  return 1;
}

# ----------------------- cdenprop ---------------------------

sub make_cdenprop_input {
  my ($r_par, $r_str) = @_;
  my $r_states = $r_par->{'data'}->{'target'}->{'states'};
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};
  my $initialsym = $r_par->{'model'}->{'initialsym'} - 1;
  my $multiplicity = $spin_multiplicity{lc($statespin)};
  my @spinsym_order = @{$r_par->{'data'}->{'scattering'}->{'spinsym_order'}};
  my $nciseti = 1;
  my $ncisetf = 1;

  &replace_in_template($r_str, "MOLECULE", $r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "LUPROP",   $luamp);
  &replace_in_template($r_str, "LUTDIP",   $luprop);
  &replace_in_template($r_str, "LUCITGT",  $lucitgt);
  &replace_in_template($r_str, "LUPROPW",  $luidip_base.$multiplicity.$statesym);
  &replace_in_template($r_str, "NTGTL",    join(",", @{$r_par->{'data'}->{'target'}->{'ntgtl'}}));
  &replace_in_template($r_str, "GFLAG",    ($r_par->{'model'}->{'model'} eq "CHF-A" ? "0" : "1"));

  if ($r_par->{'run'}->{'parallel_symm'} >= 1) {
    # traditional mode (one symetry per directory)
    &replace_in_template($r_str, "LUCSFI",   "${lucsf_base}00");
    &replace_in_template($r_str, "LUCSFF",   $lucsf_base.$multiplicity.$statesym);
    &replace_in_template($r_str, "LUCII",    "${luci}0");
    &replace_in_template($r_str, "LUCIF",    $luci);
    &replace_in_template($r_str, "NCISETI",  1);
    &replace_in_template($r_str, "NCISETF",  1);
  } else {
    # mpi-scatci mode (everything in the root directory)
    for (my $i = 0; $i < scalar @spinsym_order; $i++) {
      if (@spinsym_order[$i] eq "$multiplicity.$initialsym") { $nciseti = $i + 1; }
      if (@spinsym_order[$i] eq "$multiplicity.$statesym") { $ncisetf = $i + 1; }
    }
    &replace_in_template($r_str, "LUCSFI",   "${lucsf_base}1$initialsym");  # WARNING: only singlets assumed
    &replace_in_template($r_str, "LUCSFF",   $lucsf_base.$multiplicity.$statesym);
    &replace_in_template($r_str, "LUCII",    $luci);
    &replace_in_template($r_str, "LUCIF",    $luci);
    &replace_in_template($r_str, "NCISETI",  $nciseti);
    &replace_in_template($r_str, "NCISETF",  $ncisetf);
  }

  return 1;
}

sub make_cdenprop_all_input {
  my ($r_par, $r_str) = @_;
  my $r_states = $r_par->{'data'}->{'target'}->{'states'};

  my $nftsor = "";
  my $ntgtf = "";
  my $ntgtl = "";

  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] == 1) {
        $nftsor .= "$lucsf_base$spin_multiplicity{$statespin}$i,";
        $ntgtf .= ($i+1).",";
        $ntgtl .= "0,";
      }
    }
  }

  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "NTGT",     $r_par->{'data'}->{'nir'});
  &replace_in_template($r_str, "LUCI",     $luci);
  &replace_in_template($r_str, "LUIDIP",   $luitdip);
  &replace_in_template($r_str, "NFTINT",   $luamp);
  &replace_in_template($r_str, "NFTSOR",   $nftsor);
  &replace_in_template($r_str, "NTGTF",    $ntgtf);
  &replace_in_template($r_str, "NTGTL",    $ntgtl);
  &replace_in_template($r_str, "NUMTGT",   join(",", @{$r_par->{'data'}->{'target'}->{'ntgtl'}}));
  return 1;
}

# ------------------------ outer ----------------------------

sub make_swinterf_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};

  my $ntarg = $r_par->{'model'}->{'ntarget_states_used'};
  my @nvo = (0) x r_par->{'model'}->{'ntarget_states_used'};
  my @idtarg = @{$r_par->{'data'}->{'target'}->{'idtarg'}};

  if ($r_par->{'data'}->{'orbitals'}->{'virt_cont_used'}){
    for (my $i=0; $i<$r_par->{'model'}->{'ntarget_states_used'}; $i++) {
      my $spinsym=$r_par->{'data'}->{'target'}->{'ordered_states'}->[$i];
      my ($spin, $sym) = split(/\./, $spinsym);
      $nvo[$i]= $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$group_table[$sym]->[$statesym]];
    }
  }

  if (exists $r_par->{'model'}->{'ecutoff_target_states_in_outer'} && defined $r_par->{'model'}->{'ecutoff_target_states_in_outer'}) {
    &replace_in_template($r_str, "ECUTOFF", $r_par->{'model'}->{'ecutoff_target_states_in_outer'}/27.2114);
  } else {
    &replace_in_template($r_str, "ECUTOFF", -1.0);
  }

  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "IPOSIT",   -$r_par->{'model'}->{'positron_flag'});
  &replace_in_template($r_str, "SPIN",     $statespin);
  &replace_in_template($r_str, "SYMMETRY", $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym]);
  &replace_in_template($r_str, "MGVN",     $statesym);
  &replace_in_template($r_str, "STOT",     $spin_multiplicity{lc($statespin)});
  &replace_in_template($r_str, "NTARG",    $ntarg);
  &replace_in_template($r_str, "IDTARG",   join(",", @idtarg));
  &replace_in_template($r_str, "NVO",      join(",", @nvo));
  &replace_in_template($r_str, "RMATR",    $r_par->{'model'}->{'rmatrix_radius'});
  &replace_in_template($r_str, "RAF",      $r_par->{'model'}->{'raf'});
  &replace_in_template($r_str, "ISMAX",    $r_par->{'model'}->{'max_multipole'});
  &replace_in_template($r_str, "LCONT",    $r_par->{'data'}->{'scattering'}->{'cont_csf'});
  &replace_in_template($r_str, "USE_LCONT",$r_par->{'data'}->{'scattering'}->{'cont_csf'} ? '' : '!');
  &replace_in_template($r_str, "LUAMP",    $luamp);
  &replace_in_template($r_str, "LUCI",     $luci);
  &replace_in_template($r_str, "LUCHAN",   $luchan);
  &replace_in_template($r_str, "LURMT",    $lurmt);
  &replace_in_template($r_str, "ICFORM",   $icform);
  &replace_in_template($r_str, "IRFORM",   $irform);
  return 1;
}

sub make_rsolve_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};
  my $multiplicity = $spin_multiplicity{$statespin};
  my @spinsym_order = @{$r_par->{'data'}->{'scattering'}->{'spinsym_order'}};
  my $nset = 1;

  if ($r_par->{'run'}->{'parallel_symm'} <= 0) {
    for (my $i = 0; $i < scalar @spinsym_order; $i++) {
      if (@spinsym_order[$i] eq "$multiplicity.$statesym") {
        $nset = $i + 1;
      }
    }
  }

  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "SPIN",     $statespin);
  &replace_in_template($r_str, "SYMMETRY", $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym]);
  &replace_in_template($r_str, "MGVN",     $statesym);
  &replace_in_template($r_str, "STOT",     $spin_multiplicity{lc($statespin)});
  &replace_in_template($r_str, "NERANG",   scalar split(/\s*,\s*/, $r_par->{'model'}->{'nescat'}));
  &replace_in_template($r_str, "NESCAT",   $r_par->{'model'}->{'nescat'});
  &replace_in_template($r_str, "EINC",     $r_par->{'model'}->{'einc'});
  &replace_in_template($r_str, "IEUNIT",   $r_par->{'model'}->{'e_unit'});
  &replace_in_template($r_str, "LUCHAN",   $luchan);
  &replace_in_template($r_str, "ICFORM",   $icform);
  &replace_in_template($r_str, "LURMT",    $lurmt);
  &replace_in_template($r_str, "IRFORM",   $irform);
  &replace_in_template($r_str, "LUKMT",    "$lukmt_base$multiplicity$statesym");
  &replace_in_template($r_str, "IKFORM",   $ikform);
  &replace_in_template($r_str, "MAXF",     $r_par->{'model'}->{'maxf'});
  &replace_in_template($r_str, "NCHSET",   $nset);
  &replace_in_template($r_str, "NRMSET",   $nset);
  &replace_in_template($r_str, "NTARG",    $r_par->{'model'}->{'maxf'});
  &replace_in_template($r_str, "LU_PW_DIPOLES",    "$lupwd_base$multiplicity$statesym");
  &replace_in_template($r_str, "LU_INNER_DIPOLES", "$luidip_base$multiplicity$statesym");

  if ($r_par->{'run'}->{'photoionization'}) {
    &replace_in_template($r_str, "CALCDIP",  1);
    &replace_in_template($r_str, "IGAIL",    2);
    &replace_in_template($r_str, "RAF",      $r_par->{'model'}->{'rmatrix_radius'});
  } else {
    &replace_in_template($r_str, "CALCDIP",  0);
    &replace_in_template($r_str, "IGAIL",    1);
    &replace_in_template($r_str, "RAF",      $r_par->{'model'}->{'raf'});
  }

  return 1;
}

sub make_eigenp_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};
  my $multiplicity = $spin_multiplicity{$statespin};
  my @spinsym_order = @{$r_par->{'data'}->{'scattering'}->{'spinsym_order'}};
  my $nset = 1;

  if ($r_par->{'run'}->{'parallel_symm'} <= 0) {
    for (my $i = 0; $i < scalar @spinsym_order; $i++) {
      if (@spinsym_order[$i] eq "$multiplicity.$statesym") {
        $nset = $i + 1;
      }
    }
  }

  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "SPIN",     $statespin);
  &replace_in_template($r_str, "SYMMETRY", $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym]);
  &replace_in_template($r_str, "IEUNIT",   $r_par->{'model'}->{'e_unit'});
  &replace_in_template($r_str, "LUCHAN",   $luchan);
  &replace_in_template($r_str, "NCHSET",   $nset);
  &replace_in_template($r_str, "ICFORM",   $icform);
  &replace_in_template($r_str, "LUKMT",    "$lukmt_base$spin_multiplicity{$statespin}$statesym");
  &replace_in_template($r_str, "IKFORM",   $ikform);
  &replace_in_template($r_str, "LUPHSO",   "$lueig_base$spin_multiplicity{$statespin}$statesym");

  return 1;
}

sub make_tmatrx_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};
  my $multiplicity = $spin_multiplicity{$statespin};
  my @spinsym_order = @{$r_par->{'data'}->{'scattering'}->{'spinsym_order'}};
  my $nset = 1;

  if ($r_par->{'run'}->{'parallel_symm'} <= 0) {
    for (my $i = 0; $i < scalar @spinsym_order; $i++) {
      if (@spinsym_order[$i] eq "$multiplicity.$statesym") {
        $nset = $i + 1;
      }
    }
  }

  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "SPIN",     $statespin);
  &replace_in_template($r_str, "SYMMETRY", $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym]);
  &replace_in_template($r_str, "MAXI",     $r_par->{'model'}->{'maxi'});
  &replace_in_template($r_str, "MAXF",     $r_par->{'model'}->{'maxf'});
  &replace_in_template($r_str, "LUCHAN",   $luchan);
  &replace_in_template($r_str, "NCHSET",   $nset);
  &replace_in_template($r_str, "ICFORM",   $icform);
  &replace_in_template($r_str, "LUKMT",    "$lukmt_base$spin_multiplicity{$statespin}$statesym");
  &replace_in_template($r_str, "IKFORM",   $ikform);
  &replace_in_template($r_str, "LUTMT",    "$lutmt_base$spin_multiplicity{$statespin}$statesym");
  &replace_in_template($r_str, "ITFORM",   $itform);

  return 1;
}

sub make_ixsecs_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};
  my $multiplicity = $spin_multiplicity{$statespin};
  my @spinsym_order = @{$r_par->{'data'}->{'scattering'}->{'spinsym_order'}};
  my $nset = 1;

  if ($r_par->{'run'}->{'parallel_symm'} <= 0) {
    for (my $i = 0; $i < scalar @spinsym_order; $i++) {
      if (@spinsym_order[$i] eq "$multiplicity.$statesym") {
        $nset = $i + 1;
      }
    }
  }

  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "SPIN",     $statespin);
  &replace_in_template($r_str, "SYMMETRY", $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym]);
  &replace_in_template($r_str, "IXSN",     $r_par->{'model'}->{'x_unit'});
  &replace_in_template($r_str, "MAXI",     $r_par->{'model'}->{'maxi'});
  &replace_in_template($r_str, "MAXF",     $r_par->{'model'}->{'maxf'});
  &replace_in_template($r_str, "IEUNIT",   $r_par->{'model'}->{'e_unit'});
  &replace_in_template($r_str, "LUCHAN",   $luchan);
  &replace_in_template($r_str, "NCHSET",   $nset);
  &replace_in_template($r_str, "ICFORM",   $icform);
  &replace_in_template($r_str, "LUTMT",    "$lutmt_base$spin_multiplicity{$statespin}$statesym");
  &replace_in_template($r_str, "ITFORM",   $itform);
  &replace_in_template($r_str, "LUXSN",    "$luxsn_base$spin_multiplicity{$statespin}$statesym");

  return 1;
}

sub make_reson_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};
  my $multiplicity = $spin_multiplicity{$statespin};
  my @spinsym_order = @{$r_par->{'data'}->{'scattering'}->{'spinsym_order'}};
  my $nset = 1;

  if ($r_par->{'run'}->{'parallel_symm'} <= 0) {
    for (my $i = 0; $i < scalar @spinsym_order; $i++) {
      if (@spinsym_order[$i] eq "$multiplicity.$statesym") {
        $nset = $i + 1;
      }
    }
  }

  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "SPIN",     $statespin);
  &replace_in_template($r_str, "SYMMETRY", $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$statesym]);
  &replace_in_template($r_str, "LUCHAN",   $luchan);
  &replace_in_template($r_str, "ICFORM",   $icform);
  &replace_in_template($r_str, "NCHSET",   $nset);
  &replace_in_template($r_str, "LURMT",    $lurmt);
  &replace_in_template($r_str, "IRFORM",   $irform);
  &replace_in_template($r_str, "NRMSET",   $nset);
  &replace_in_template($r_str, "LUKMT",    "$lukmt_base$spin_multiplicity{$statespin}$statesym");
  &replace_in_template($r_str, "IKFORM",   $ikform);
  &replace_in_template($r_str, "LURES",    "$lures_base$spin_multiplicity{$statespin}$statesym");
  &replace_in_template($r_str, "MGVN",     $statesym);
  &replace_in_template($r_str, "STOT",     $spin_multiplicity{lc($statespin)});

  if ($r_par->{'run'}->{'photoionization'}) {
    &replace_in_template($r_str, "IGAIL",  2);
    &replace_in_template($r_str, "RAF",    $r_par->{'model'}->{'rmatrix_radius'});
  } else {
    &replace_in_template($r_str, "IGAIL",  1);
    &replace_in_template($r_str, "RAF",    $r_par->{'model'}->{'raf'});
  }

  return 1;
}

sub make_time_delay_input {
  my ($r_par, $r_str) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $statespin = $r_par->{'data'}->{$task}->{'spin'};

  &replace_in_template($r_str, "LUKMT", "$lukmt_base$spin_multiplicity{$statespin}$statesym");

  return 1;
}

# ------------------------ dipelm ----------------------------

sub make_dipelm_input {
  my ($r_par, $r_str) = @_;

  my $ntarg = $r_par->{'model'}->{'maxf'} > 0 ?
              $r_par->{'model'}->{'maxf'} :
              $r_par->{'model'}->{'ntarget_states_used'};

  &replace_in_template($r_str, "NTARG", $ntarg);
  if ($r_par->{'run'}->{'dipelm_smooth'} eq 0){
    &replace_in_template($r_str, "SMOOTH", ".false.");
  } else {
    &replace_in_template($r_str, "SMOOTH", ".true.");
  }

  # IP must be in atomic units for dipelm
  my $Ip = $r_par->{'model'}->{'first_Ip'}*1.0;
  if ($r_par->{'model'}->{'e_unit'} == 1) { $Ip = $Ip/2.0; }
  if ($r_par->{'model'}->{'e_unit'} == 2) { $Ip = $Ip/27.211386246; } # value from NIST
  $Ip = sprintf("%.5f", $Ip);

  &replace_in_template($r_str, "IP", $Ip);

  my %irreps_used;
  my @units;
  my @nsets;
  for (my $k = 0; $k < 3; $k++) {
    my $irr = $coor_repr{$r_par->{'model'}->{'symmetry'}}->[$k];
    if(not(exists($irreps_used{$irr}))){     # only use each irrep once
      $irreps_used{$irr} = undef;            # add irr to the list of used irreps
      push(@units, "${lupwd_base}1${irr}");
      push(@nsets, 1);
    }
  }

  &replace_in_template($r_str, "LU_PW_DIPS", join(",", @units));
  &replace_in_template($r_str, "NS_PW_DIPS", join(",", @nsets));

  return 1;
}

# ------------------------ rmt_interface ----------------------------

sub make_rmt_interface_input {
  my ($r_par, $r_str) = @_;

  my $spin = "";
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] == 1) {
        $spin = $statespin;
      }
    }
  }

  &replace_in_template($r_str, "NSYM",     $r_par->{'data'}->{'nir'});
  &replace_in_template($r_str, "MOLECULE", "e + ".$r_par->{'model'}->{'molecule'});
  &replace_in_template($r_str, "SPIN",     $spin);
  &replace_in_template($r_str, "SYMMETRY", "");
  &replace_in_template($r_str, "STOT",     $spin_multiplicity{lc($spin)});
  &replace_in_template($r_str, "NTARG",    $r_par->{'model'}->{'ntarget_states_used'});
  &replace_in_template($r_str, "IDTARG",   join(",", @{$r_par->{'data'}->{'target'}->{'idtarg'}}));
  &replace_in_template($r_str, "RMATR",    $r_par->{'model'}->{'rmatrix_radius'});
  &replace_in_template($r_str, "RAF",      $r_par->{'model'}->{'raf'});
  &replace_in_template($r_str, "ISMAX",    $r_par->{'model'}->{'max_multipole'});
  &replace_in_template($r_str, "LUTARG",   $luprop);
  &replace_in_template($r_str, "LUITDIP",  $luitdip);
  &replace_in_template($r_str, "LUAMP",    $luamp);
  &replace_in_template($r_str, "LUCI",     $luci);

  return 1;
}

# ============= end of subroutines for INPUTS ===============

# ================= subroutines for OUTPUTS =================

# ------------------------- molpro --------------------------

# From molpro the energies of the target molecular orbitals are read in
#   and saved in $r_par->{'data'}->{'orbitals'}->{'energies'}
# the number of occupied orbitals is determined from molpro output
# then the number of frozen, active and virtual orbitals
#   in each IR is determined and the script continues with the next program
# WARNING: at the moment this routine does not handle open shell cases and does not check
#          whether molpro experienced problems with convergence!!!
sub read_molpro_output {
  my ($r_par) = @_;
  my $sym = 0;
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR
  my $r_e = $r_par->{'data'}->{'orbitals'}->{'energies'};
  my $model = $r_par->{'model'}->{'model'};
  my $header = "NATURAL ORBITALS";
  my $nopenshell = 0;
  if ($r_par->{'data'}->{'scf_ok'} == 0 || $model =~ m/^SE/ || $r_par->{'model'}->{'orbitals'} eq "HF") {
    $header = "ELECTRON ORBITALS";
  }
  &print_info("Looking for $header in MOLPRO output...\n", $r_par);

  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    my $look_for_mos = 0;   # auxiliary variable to decide whether to look for orbital energies
    my $in_mos_section = 0;

    $r_par->{'data'}->{'target'}->{'noccupied'} = 0;

    while (my $line = <OUTPUT>) {
      chomp($line);
      if ($header eq "ELECTRON ORBITALS") { # next three items only for HF calculations
        if ($line =~ m/\s*Final (alpha )?occupancy:\s+(\d.*)$/) {
          $r_par->{'data'}->{'orbitals'}->{'occupied'} = [ split(/\s+/, $2) ];
          &print_info("  Occupancy from the Molpro output: ".join(", ", @{$r_par->{'data'}->{'orbitals'}->{'occupied'}})."\n", $r_par);
        }
        elsif ($line =~ m/\s*.RHF STATE\s*(\d+)\.(\d+)\s*Energy\s*([\+\-]?\d+\.\d+)/) {
          $r_par->{'data'}->{'target'}->{'hf_energy'} = $3;
          $r_par->{'data'}->{'target'}->{'hf_symmetry'} = $2 - 1;
          &print_info("  Hartree-Fock energy = $3\n", $r_par);
        }
        elsif ($line =~ m/^\s*NUMBER OF CONTRACTIONS:\s*(\d+)\s*\(\s*([^\)]+)\)\s*$/) { #look for NOB in MOLPRO for INTEGRALS input
          $r_par->{'data'}->{'orbitals'}->{'target_sum'} = $1;
          &check_numbers_of_orbitals($r_par);
          my $str_orbs = $2;
          for (my $i = 0; $i < $nir; $i++) {
            $line =~ s/\s*(\d+)[A-Z][^\+]+//;
            $r_par->{'data'}->{'orbitals'}->{'target_all'}->[$i] = $1;
          }
          &print_info("  Total number of orbitals for each symmetry: ".join(", ", @{$r_par->{'data'}->{'orbitals'}->{'target_all'}})."\n", $r_par);
        }
      }
      if ($line =~ /$header/) { #the data that follow contain the orbital data
        &print_info("  Found $header...\n", $r_par);
        $in_mos_section = 1;
      }
      elsif ($line =~ m/DATASETS/) {
        $in_mos_section = 0; #string 'DATASETS' mark the end of the orbital data
      }
      elsif ($line =~ m/^\s*\? Too many active orbitals/) {
        &print_info("  Error: occurred when running Molpro, probably too many active orbitals !\n", $r_par);
        &print_info("         try to use HF orbitals insted of natural orbitals: set 'orbitals' to 'HF' in model.pl\n", $r_par);
        die "  Error: Molpro failed, probably too many active orbitals !\n  Try to set 'orbitals' to 'HF' in model.pl\n";
      }

      if ($in_mos_section == 1) { #look for the orbital energies if we are in the section containg the orbital coefficients.
        if ($line =~ /^\s*$/) { #we only look for the info on the orbital energies on the line following a blank line.
          $look_for_mos = 1;
        }
        if ($look_for_mos == 1 && $line =~ m/^\s*(\d+)\.(\d)\s+(\d|\+|\-|[\+\-]?\d+\.\d+)\s+([\+\-]?\d+\.\d+)\s+/) {
          $iorb = $1;
          $sym = $2;
          $occ = $3;
          $r_e->{"$sym.$iorb"} = $4;
          $val = $r_e->{"$sym.$iorb"};
          if ($occ eq "+" or $occ eq "-") {
             $nopenshell++;
          }
          if ($occ eq "+" or $occ eq "-" or $occ != 0) { #occupancy of a particular orbital is only YES/NO
             $r_par->{'data'}->{'target'}->{'noccupied'}++;
          }
          $look_for_mos = 0;
        }
      }
    }
    close(OUTPUT);
  }
  else {
    die "  Error: no output of Molpro !\n";
  }

  if ($nopenshell > 0) {
    &print_info("  Number of singly occupied HF orbitals: $nopenshell\n", $r_par);
  }

  # Checking occupied orbitals
  my @sorted_orb;
  if ($r_par->{'model'}->{'use_MAS'} || $r_par->{'model'}->{'use_MASSCF'}){
    # Get the orbitals from the molden files instead of molpro output so that
    # we have the right ordering (subspace-occupation-energy)
    my $prefix = lc($r_par->{'model'}->{molecule});
    my $file = "$prefix.molden";
    my $point_group = $r_par->{'model'}->{'symmetry'};
    my @sorted_orbitals = MultiSpace::read_orbitals_from_molden("molpro", $file, $point_group);

    # If MASSCF is not used allow re-ordering by energy if requested.
    if ($r_par->{'model'}->{'select_orb_by'} eq "energy"){
      if (!$r_par->{'model'}->{'use_MASSCF'}) {
        @sorted_orbitals = sort {$a->{'Ene='} <=> $b->{'Ene='}} @sorted_orbitals;
      } else {
        die "Error: 'select_orb_by' must be set to 'molden' when MASSCF is used!";
      }
    }

    foreach my $orb (@sorted_orbitals){
      push(@sorted_orb, "$orb->{'Sym='}.$orb->{'Idx_in_sym='}");
    }
    $r_par->{'data'}->{'MAS'}->{'sorted_orb'} = \@sorted_orb; # Just the sym and index
    $r_par->{'data'}->{'MAS'}->{'sorted_orbitals'} = \@sorted_orbitals; # Everything
    my $r_occ = [];
    &get_number_of_mos_automatically($r_occ, $nir, \@sorted_orb,
                                   1, $r_par->{'data'}->{'target'}->{'noccupied'});
    if ($r_par->{'model'}->{'norbitals_to_print'} > 0) {
      &print_info("  Occupancy and energy of the first $r_par->{'model'}->{'norbitals_to_print'} molecular orbitals: in X.Y  X stands for IR and Y counts MOs\n", $r_par);
      my $p = 0;
      foreach my $orb (@sorted_orbitals){
        &print_info("    $orb->{'Sym='}.$orb->{'Idx_in_sym='}    $orb->{'Occup='}    $orb->{'Ene='}\n", $r_par);
        if ($r_par->{'model'}->{'norbitals_to_print'} == $p++) { last; }
      }
    }
    # Read the CASSCF (or ORMASSCF) state energies for later comparison to state energies
    # produced in the target scatci run.
    if (!($r_par->{'data'}->{'scf_ok'} == 0 || $model =~ m/^SE/ || $r_par->{'model'}->{'orbitals'} eq "HF")) {
      read_state_energies_from_molpro($r_par);
    }

  } else {
    @sorted_orb = sort { $r_e->{$a} <=> $r_e->{$b} } keys %$r_e;
    my $r_occ = [];
    &get_number_of_mos_automatically($r_occ, $nir, \@sorted_orb,
                                    1, $r_par->{'data'}->{'target'}->{'noccupied'});

    if ($r_par->{'model'}->{'norbitals_to_print'} > 0) {
      &print_info("  Energies of the first $r_par->{'model'}->{'norbitals_to_print'} molecular orbitals: in X.Y  X stands for IR and Y counts MOs\n", $r_par);
      my $p = 0;
      foreach my $key (@sorted_orb) {
        &print_info("    $key    $r_e->{$key}\n", $r_par);
        if ($r_par->{'model'}->{'norbitals_to_print'} == $p++) { last; }
      }
    }
  }
  $r_par->{'data'}->{'scf_ok'} = 1;
  $r_par->{'data'}->{'target'}->{'nfull_occ'} = int($r_par->{'model'}->{'nelectrons'} / 2); # number of fully occupied orbitals
  $r_par->{'data'}->{'target'}->{'nhalf_occ'} = $r_par->{'model'}->{'nelectrons'} % 2;      # half occupied orbital

  &set_orbitals($r_par, \@sorted_orb);
  return 1;
}

sub read_state_energies_from_molpro{
    my ($r_par) = @_;
    my $file = $r_par->{'data'}->{'outputfile'};

    open my $in,  '<',   $file or die "Could not open $file: $!";

    my %states;
    my %nocsfs;
    my @statespins;
    while( my $line = <$in>)  {

      if ($line =~ m/\s*Number of electrons:\s*\d+\s*Spin symmetry=\s*(\w+)\s*Space symmetry=(\d)/){

        my $statesym = $2 - 1;
        my $statespin = $spin_multiplicity{lc($1)};
        $line = <$in>;
        if ($line =~ m/\s*Number of states:\s*(\d+)\s*/){
          for my $id (1..$1){
              push(@statespins, $statespin)
          }
        }
        $line = <$in>;
        if ($line =~ m/\s*Number of CSFs:\s*(\d+)\s*/){
          $nocsfs{"$statespin.$statesym"} = $1;
        }
      }

      if ($line =~ m/\s*.MCSCF STATE\s*(\d+)\.(\d+)\s*Energy\s*([\+\-]?\d+\.\d+)/){
        my $id = $1;
        my $statesym = $2 - 1;
        my $statespin = shift @statespins;
        $states{"$statespin.$statesym.$id"} = $3;
      }
    }

    close $in;

    $r_par->{'data'}->{'qchem'}->{'states'} = \%states;
    $r_par->{'data'}->{'qchem'}->{'nocsfs'} = \%nocsfs;

    return 1;
}
# ------------------------- psi4 --------------------------

sub read_psi4_output {
  my ($r_par) = @_;
  my $nir = $r_par->{'data'}->{'nir'};
  my $r_e = $r_par->{'data'}->{'orbitals'}->{'energies'};
  my $grp = $r_par->{'model'}->{'symmetry'};
  my $occupation = -1;

  $r_par->{'data'}->{'target'}->{'noccupied'} = 0;
  $r_par->{'data'}->{'orbitals'}->{'target_sum'} = 0;

  for (my $j = 0; $j < $nir; $j++) {
    $r_par->{'data'}->{'orbitals'}->{'target_all'}->[$j] = 0;
  }

  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    while (my $line = <OUTPUT>) {
      chomp($line);
      if    ($line =~ m/\s*Doubly Occupied:\s*$/)           { $occupation =  2; }
      elsif ($line =~ m/\s*Singly Occupied:\s*$/)           { $occupation =  1; }
      elsif ($line =~ m/\s*Virtual:\s*$/)                   { $occupation =  0; }
      elsif ($line =~ m/\s*Final Occupation by Irrep:\s*$/) {             last; }
      elsif ($occupation >= 0) {
        my @tokens = split ' ', $line;
        my $ntokens = @tokens;
        for (my $i = 0; $i < $ntokens/2; $i++) {
          my $orb = $tokens[2*$i+0];  $orb =~ s/(\d+)(.*)/\1/;
          my $irr = $tokens[2*$i+0];  $irr =~ s/(\d+)(.*)/\2/;
          my $ene = $tokens[2*$i+1];
          for (my $j = 0; $j < $nir; $j++) {
            my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$j];
            if (lc($str_ir) eq lc($irr)) {
              my $sym = 1 + $j;
              $r_e->{"$sym.$orb"} = $ene;
              if ($occupation > 0) { $r_par->{'data'}->{'target'}->{'noccupied'}++; }
              $r_par->{'data'}->{'orbitals'}->{'target_sum'}++;
              $r_par->{'data'}->{'orbitals'}->{'target_all'}->[$j]++;
              last;
            }
          }
        }
      }
    }
    close(OUTPUT);
  }

  &print_info("  Total number of orbitals for each symmetry: ".join(", ", @{$r_par->{'data'}->{'orbitals'}->{'target_all'}})."\n", $r_par);

  # Checking occupied orbitals
  my @sorted_orb = sort { $r_e->{$a} <=> $r_e->{$b} } keys %$r_e;
  my $r_occ = [];
  &check_numbers_of_orbitals($r_par);
  &get_number_of_mos_automatically($r_occ, $nir, \@sorted_orb,
                                   1, $r_par->{'data'}->{'target'}->{'noccupied'});

  if ($r_par->{'model'}->{'norbitals_to_print'} > 0) {
    &print_info("  Energies of the first $r_par->{'model'}->{'norbitals_to_print'} molecular orbitals: in X.Y  X stands for IR and Y counts MOs\n", $r_par);
    my $p = 0;
    foreach my $key (@sorted_orb) {
      &print_info("    $key    $r_e->{$key}\n", $r_par);
      if ($r_par->{'model'}->{'norbitals_to_print'} == $p++) { last; }
    }
  }

  $r_par->{'data'}->{'scf_ok'} = 1;
  $r_par->{'data'}->{'target'}->{'nfull_occ'} = int($r_par->{'model'}->{'nelectrons'} / 2); # number of fully occupied orbitals
  $r_par->{'data'}->{'target'}->{'nhalf_occ'} = $r_par->{'model'}->{'nelectrons'} % 2;      # half occupied orbital

  &set_orbitals($r_par, \@sorted_orb);
  return 1;
}

# ------------------------- molcas --------------------------
# TODO: If frozen orbitals are used then energies of those orbital
#       must be taken from the scf.molden file as they are zero in the
#       rasscf.molden file. This is not implemented yet.
#       At the moment, despite the confusing 'frozen_orbs'
#       key name, frozen orbs as they are defined in molpro/molcas are not used
#       in the scripts, 'frozen_orbs' infact refers to closed/inactive orbs.
#       However for the large calcs that will be attempted with molcas it may
#       be desirable to use frozen orbitals.
sub read_molcas_output {
  my ($r_par) = @_;
  my $nir = $r_par->{'data'}->{'nir'};
  my $r_e = $r_par->{'data'}->{'orbitals'}->{'energies'};
  my $point_group = $r_par->{'model'}->{'symmetry'};

  $r_par->{'data'}->{'target'}->{'noccupied'} = 0;
  $r_par->{'data'}->{'orbitals'}->{'target_sum'} = 0;

  for (my $j = 0; $j < $nir; $j++) {
    $r_par->{'data'}->{'orbitals'}->{'target_all'}->[$j] = 0;
  }

  my $prefix = lc($r_par->{'model'}->{molecule});
  my $file = "$prefix.molden";

  my @orbitals = MultiSpace::read_orbitals_from_molden("molcas", $file, $point_group);
  foreach my $orb (@orbitals){
    $r_par->{'data'}->{'orbitals'}->{'target_all'}->[$orb->{'Sym='}-1]++;
    $r_par->{'data'}->{'orbitals'}->{'target_sum'}++;
    if ($orb->{'Occup='} > 0) {
      $r_par->{'data'}->{'target'}->{'noccupied'}++;
    }
  }

  my $subspaces;
  if ($r_par->{'model'}->{'model'} =~ /^SE/ || $r_par->{'data'}->{'scf_ok'} == 0 || $r_par->{'model'}->{'orbitals'} eq "HF") {
    $subspaces = [$r_par->{'data'}->{'orbitals'}->{'target_all'}];
  } else {
    $subspaces = $r_par->{'data'}->{'MAS'}->{'qchem'}->{'subspaces'};
  }

  my @sorted_orbitals = MultiSpace::sort_molden_orbitals(orbitals => \@orbitals, subspaces => $subspaces);

  my @sorted_orb;
  foreach my $orb (@sorted_orbitals){
    push(@sorted_orb, "$orb->{'Sym='}.$orb->{'Idx_in_sym='}");
  }
  $r_par->{'data'}->{'MAS'}->{'sorted_orb'} = \@sorted_orb;

  &print_info("  Total number of orbitals for each symmetry: ".join(", ", @{$r_par->{'data'}->{'orbitals'}->{'target_all'}})."\n", $r_par);

  # Checking occupied orbitals
  my $r_occ = [];
  &check_numbers_of_orbitals($r_par);
  &get_number_of_mos_automatically($r_occ, $nir, \@sorted_orb,
                                   1, $r_par->{'data'}->{'target'}->{'noccupied'});

  if ($r_par->{'model'}->{'norbitals_to_print'} > 0) {
    &print_info("  Occupancy and energy of the first $r_par->{'model'}->{'norbitals_to_print'} molecular orbitals: in X.Y  X stands for IR and Y counts MOs\n", $r_par);
    my $p = 0;
    foreach my $orb (@sorted_orbitals){
      &print_info("    $orb->{'Sym='}.$orb->{'Idx_in_sym='}    $orb->{'Occup='}    $orb->{'Ene='}\n", $r_par);
      if ($r_par->{'model'}->{'norbitals_to_print'} == $p++) { last; }
    }
  }

  $r_par->{'data'}->{'scf_ok'} = 1;
  $r_par->{'data'}->{'target'}->{'nfull_occ'} = int($r_par->{'model'}->{'nelectrons'} / 2); # number of fully occupied orbitals
  $r_par->{'data'}->{'target'}->{'nhalf_occ'} = $r_par->{'model'}->{'nelectrons'} % 2;      # half occupied orbital

  &set_orbitals($r_par, \@sorted_orb);
  return 1;
}
# ------------------------- integrals --------------------------

sub read_scatci_integrals_output {
  my ($r_par) = @_;
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR

  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    my $found_orbs = 0;
    while ($line = <OUTPUT>) {
      chomp($line);
      if ($line =~ m/\s*Final number of orbitals for target/) {
        $found_orbsi = 1;
        &print_info("Final number of orbitals for target: ", $r_par);
        $line = <OUTPUT>;
        for (my $i = 0; $i < $nir; $i++) {
          $line =~ s/\s*(\d+)//;
          $r_par->{'data'}->{'orbitals'}->{'TGT'}->[$i] = $1;
          &print_info("$1, ", $r_par);
        }
        &print_info("\n", $r_par);
      }
      if ($line =~ m/\s*Final number of orbitals for PCOs/) {
        $found_orbsi = 1;
        &print_info("Final number of orbitals for PCOs: ", $r_par);
        $line = <OUTPUT>;
        for (my $i = 0; $i < $nir; $i++) {
          $line =~ s/\s*(\d+)//;
          $r_par->{'data'}->{'orbitals'}->{'PCO'}->[$i] = $1;
          &print_info("$1, ", $r_par);
        }
        &print_info("\n", $r_par);
      }
      if ($line =~ m/\s*Final number of orbitals for TGT\+PCOs/) {
        $found_orbsi = 1;
        &print_info("Final number of orbitals for TGT\+PCOs: ", $r_par);
        $line = <OUTPUT>;
        for (my $i = 0; $i < $nir; $i++) {
          $line =~ s/\s*(\d+)//;
          $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] = $1;
          &print_info("$1, ", $r_par);
        }
        &print_info("\n", $r_par);
      }
      if ($line =~ m/\s*Final number of orbitals for continuum/) {
        $found_orbsi = 1;
        &print_info("Final number of orbitals for continuum: ", $r_par);
        $line = <OUTPUT>;
        for (my $i = 0; $i < $nir; $i++) {
          $line =~ s/\s*(\d+)//;
          $r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$i] = $1;
          &print_info("$1, ", $r_par);
        }
        &print_info("\n", $r_par);
      }
      if ($line =~ m/\s*Final number of orbitals for total/) {
        $found_orbs = 1;
        &print_info("Number of molecular orbitals in each irreducible representation: ", $r_par);
        $line = <OUTPUT>;
        for (my $i = 0; $i < $nir; $i++) {
          $line =~ s/\s*(\d+)//;
          $r_par->{'data'}->{'orbitals'}->{'all'}->[$i] = $1;
          &print_info("$1, ", $r_par);
        }
        &print_info("\n", $r_par);
      }
    }
    close(OUTPUT);
    if ($found_orbs == 0) {
      &print_info("\nError: There were no information about molecular orbitals found in scatci_integrals output!\n");
      &print_info("Exiting! Please check inputs and outputs to find out what is wrong.\n");
      die "  Error: no information about molecular orbitals !\n";
    }
  }
  else {
    die "  Error: no output of scatci_integrals !\n";
  }

  # We also determine number of used continuum orbitals = all - target_all (= all_used - target)
  for (my $i = 0; $i < $nir; $i++) {
    $r_par->{'data'}->{'orbitals'}->{'target_all'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i];
    $r_par->{'data'}->{'orbitals'}->{'cont_all'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$i];
  }
  &print_info("  Cont used orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'cont_used'}})."\n", $r_par);

  # All used orbitals
  for (my $i = 0; $i < $nir; $i++) {
    $r_par->{'data'}->{'orbitals'}->{'all_used'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'all'}->[$i];  #$r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] + $r_par->{'data'}->{'orbitals'}->{'cont_used'}->[$i];
  }

  for (my $multi2 = 0; $multi2 < $nir; $multi2++) {
    my $a = $r_par->{'data'}->{'orbitals'}->{'all_used'}->[$multi2];
    $r_par->{'data'}->{'orbitals'}->{'all_used_positron'}->[$multi2] = $a * 2;
    $r_par->{'data'}->{'orbitals'}->{'empty'}->[$multi2] = 0;
  }

  &print_info("  All  used orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'all_used'}})."\n", $r_par);

  return 1;
}

# This auxiliary subroutine sets and prints number of various
# orbitals in each irreducible representation
sub set_orbitals {
  my ($r_par, $r_sorted_orb) = @_;
  my $nir = $r_par->{'data'}->{'nir'};

  $r_par->{'data'}->{'MAS'}->{'sorted_orb'} = $r_sorted_orb;

  # Frozen orbitals
  my $starting_orb = 1;
  &get_number_of_mos('frozen', $r_par, $r_sorted_orb, $starting_orb);
  &print_info("  Frozen orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'frozen'}})."\n", $r_par);

  # Active orbitals
  $starting_orb += $r_par->{'model'}->{'nfrozen'};
  &get_number_of_mos('active', $r_par, $r_sorted_orb, $starting_orb);
  &print_info("  Active orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'active'}})."\n", $r_par);

  # Target orbitals = frozen + active
  for (my $i = 0; $i < $nir; $i++) {
    $r_par->{'data'}->{'orbitals'}->{'target'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'frozen'}->[$i]
                                                     + $r_par->{'data'}->{'orbitals'}->{'active'}->[$i];
  }
  &print_info("  Target  orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'target'}})."\n", $r_par);

  # Virtual orbitals
  $starting_orb += $r_par->{'model'}->{'nactive'};
  &get_number_of_mos('virtual', $r_par, $r_sorted_orb, $starting_orb);
  &print_info("  Virtual orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'virtual'}})."\n", $r_par);

  # Target used orbitals = target + virtual
  for (my $i = 0; $i < $nir; $i++) {
    $r_par->{'data'}->{'orbitals'}->{'target_used'}->[$i] = $r_par->{'data'}->{'orbitals'}->{'target'}->[$i]
                                                          + $r_par->{'data'}->{'orbitals'}->{'virtual'}->[$i];
  }
  &print_info("  Used    orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'target_used'}})."\n", $r_par);

  # Reference orbitals
  $starting_orb = 1;
  &get_number_of_mos('reference', $r_par, $r_sorted_orb, $starting_orb);
  &print_info("  Reference orbitals: ".join(",", @{$r_par->{'data'}->{'orbitals'}->{'reference'}})."\n", $r_par);
  my $norbGS = int( $r_par->{'model'}->{'nelectrons'} / 2 ) + $r_par->{'model'}->{'nelectrons'}%2;
  for (my $i = 0; $i < $norbGS ; $i++) {
    $r_sorted_orb->[$i] =~ m/(\d+)\.(\d+)/;
    if( $1 > 0 ){ ${$r_par->{'data'}->{'orbitals'}->{'GS'}}[$1 - 1]++; }
  }

  return 1;
}

# This auxiliary subroutine determines number of orbitals
# in each irreducible representation
# $which_orbitals = frozen, active, virtual, or reference
# it is assumed that $r_sorted_orb contains sorted molecular orbitals from HF or CASSCF calculations
# $starting_orb is the number of the first orbital to use in $r_sorted_orb
sub get_number_of_mos {
  my ($which_orbitals, $r_par, $r_sorted_orb, $starting_orb) = @_;
  my $nir = $r_par->{'data'}->{'nir'};

  my @given_orbs = @{$r_par->{'model'}->{$which_orbitals.'_orbs'}}[0..$nir-1];
  my $ngiven_orbs = 0;
  $ngiven_orbs += $_ for @given_orbs;

  my $get_orbitals_automatically = 1;
  # if user specified orbitals are consistent with the total number of orbitals then ...
  if ($ngiven_orbs > 0) {
    if ($ngiven_orbs == $r_par->{'model'}->{'n'.$which_orbitals}) {
      $get_orbitals_automatically = 0;
      @{$r_par->{'data'}->{'orbitals'}->{$which_orbitals}} = @given_orbs;
    }
  }
  # otherwise orbitals are determined automatically ...
  if ($get_orbitals_automatically == 1) {
    &get_number_of_mos_automatically($r_par->{'data'}->{'orbitals'}->{$which_orbitals}, $nir, $r_sorted_orb,
                                     $starting_orb, $starting_orb - 1 + $r_par->{'model'}->{'n'.$which_orbitals});
  }
  return 1;
}

# This auxiliary subroutine determines automatically number of orbitals
# in each irreducible representation to an array specified using reference $r_orb
# starting from the $first_orb to the $last_orb orbital ordered by energy
sub get_number_of_mos_automatically {
  my ($r_orb, $nir, $r_sorted_orb, $first_orb, $last_orb) = @_;
  for (my $i = 0; $i < $nir; $i++) { $r_orb->[$i] = 0; }
  for (my $i = $first_orb - 1; $i < $last_orb; $i++) {
    $r_sorted_orb->[$i] =~ m/(\d+)\.(\d+)/;
    $r_orb->[$1 - 1]++;
  }
  return 1;
}


# ------------------------ congen --------------------------

# In congen output we only check that some configuartions were generated
# if there is no configuration
#   $r_par->{'data'}->{'no_scatci'} is set to 1
# to skip scatci calculation
sub read_congen_output {
  my ($r_par) = @_;

  my $str = "";
  if (read_file($r_par->{'data'}->{'outputfile'}, \$str)) {
    if ($str =~ m/TOTAL NUMBER OF CSF\'S GENERATED IS\s+(\d+)/is) {
      &print_info("  Total number of CSF's generated is $1\n", $r_par);
      my $ncsfs = $1;
      my $task = $r_par->{'data'}->{'task'};
      if ($r_par->{'model'}->{'use_MAS'} > 0 && $task eq "target") {
        my $sym = $r_par->{'data'}->{$task}->{'symmetry'};
        my $spin = $spin_multiplicity{$r_par->{'data'}->{$task}->{'spin'}};
        $r_par->{'data'}->{'MAS'}->{'nocsfs'}->{"$spin.$sym"} = $ncsfs;
      }
      if ($r_par->{'model'}->{'model'} =~ /(^CHF)/ && $task eq "target") {
        my $sym = $r_par->{'data'}->{$task}->{'symmetry'};
        my $spin = $r_par->{'data'}->{$task}->{'spin'};
        $state_id = $r_par->{'data'}->{'CHF'}->{'current_state'}-1;
        $r_par->{'data'}->{'CHF'}->{'states'}->[$state_id]->{'no_csfs'} = $ncsfs;
        $r_par->{'model'}->{'ntarget_states'}->{$spin}->[$sym] = $ncsfs;
      }
      if ($ncsfs == 0) {
        $r_par->{'data'}->{'no_scatci'} = 1;
        &print_info("\n  Warning:  Other runs for this spin-symmetry will be skipped!\n", $r_par);
        &print_info("\n            The number of target states given by the user can be different\n", $r_par);
        &print_info("\n            from the actual number of target states used in scattering calculations!\n", $r_par);
      }
      else {
        &print_info(".\n", $r_par);
      }
    }
    elsif ($str =~ m/ERROR IN REFORB DATA/is) {
      $r_par->{'data'}->{'no_scatci'} = 1;
      &print_info("  Error in REFORB in congen input !!!\n  Probably there is no available orbital of required symmetry !\n  Other runs for this spin-symmetry will be skipped!\n", $r_par);
    }
    else {
      die "  Error: incorrect output of congen: file $r_par->{'data'}->{'outputfile'} !\n";
    }
  }
  else {
    die "  Error: no output of congen: file $r_par->{'data'}->{'outputfile'} !\n";
  }

  return 1;
}

# ------------------------ scatci --------------------------

# From scatci output we can determine energies of states
#   saved in $r_par->{'data'}->{'target'}->{'states'}
#   and   in $r_par->{'data'}->{$task}->{'states_all_geom'}
sub read_scatci_output {
  my ($r_par) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $statespin = $spin_multiplicity{lc($r_par->{'data'}->{$task}->{'spin'})};
  my $statesym  = $r_par->{'data'}->{$task}->{'symmetry'};
  my $igeom = $r_par->{'data'}->{'igeom'} - 1; # index of geometry for storage
  my $offset = 0;

  if ($r_par->{'model'}->{'model'} eq "CHF-A" && $task eq "target") {
    foreach my $key (keys %{$r_par->{'data'}->{'target'}->{'states'}}){
      if ($key =~ /$statespin.$statesym.\d+/){
        $offset++;
      }
    }
  }

  my $str = "";
  if (read_file($r_par->{'data'}->{'outputfile'}, \$str)) {
    if ($str =~ m/^.*EIGEN-ENERGIES\s+(([\+\-]?\d+\.\d+\s*)+)/is) {
      $energies_list = $1;
      $energies_list =~ s/\s*$//;
      my @energies = split(/\s+/, $energies_list);
      my $ne = scalar @energies;
      &print_info("  Found $ne energies");
      if ($ne > $r_par->{'model'}->{'nstates_to_print'}) {
        &print_info(", printing first $r_par->{'model'}->{'nstates_to_print'} energies:\n", $r_par);
      }
      else {
        &print_info(":\n", $r_par);
      }
      for (my $i = 1; $i <= $ne; $i++) {
        if ($task eq "target") {
          my $id = $i + $offset;
          $r_par->{'data'}->{'target'}->{'states'}->{"$statespin.$statesym.$id"} = $energies[$i-1];
        }
        if ($i <= $r_par->{'model'}->{'nstates_to_print'}) {
          &print_info("  $statespin.$statesym.$i -> ".$energies[$i-1]."\n", $r_par);
        }
        if ($igeom > 0) {
          $r_par->{'data'}->{$task}->{'states_all_geom'}->{"$statespin.$statesym.$i"}->[$igeom] = $energies[$i-1];
        }
        else {
          $r_par->{'data'}->{$task}->{'states_all_geom'}->{"$statespin.$statesym.$i"} = [];
          $r_par->{'data'}->{$task}->{'states_all_geom'}->{"$statespin.$statesym.$i"}->[$igeom] = $energies[$i-1];
        }
      }
    }
    else {
      die "  Error: incorrect output of scatci: file $r_par->{'data'}->{'outputfile'} !\n";
    }
    if ($r_par->{'run'}->{'parallel_diag'} == 0 && $task eq 'scattering') {
      my $nocsf = -1, $mocsf = -1, $ncont = -1;
      if ($str =~ m/.*NOCSF =\s*(\d+)\s*prototype CSFs.*/is) { $nocsf = $1; }
      if ($str =~ m/.*MOCSF =\s*(\d+)\s*dimension final Hamiltonian.*/is) { $mocsf = $1; }
      if ($str =~ m/.*Number of last continuum CSF, NCONT =\s*(\d+).*/is) { $ncont = $1; }
      if (mocsf == -1 or nocsf == -1 or ncont == -1) {
        die "  Error: incorrect output of scatci: file $r_par->{'data'}->{'outputfile'} (NOCSF = $nocsf, MOCSF = $mocsf, NCONT = $ncont)!\n";
      } else {
        # LCONT = "Hamiltonian size"
        my $lcont = 0;
        #A SUBSET OF THE CI VECTORS COULD ONLY BE STORED BY MPI-SCATCI
        if (not $r_par->{'run'}->{'mpi_scatci'} eq "") {
           #A SUBSET OF THE CI VECTORS IS SAVED ONLY IN CASE OF SCATTERING CALCULATIONS AND WHEN RMT_INTERFACE IS NOT REQUIRED
           if (not $r_par->{'run'}->{'photoionization'} and not $r_par->{'run'}->{'rmt_interface'}) {
              # LCONT = "Hamiltonian size" minus "number of L2 configurations"
              $lcont = $mocsf - ($nocsf - $ncont);
           }
        }
        $r_par->{'data'}->{'scattering'}->{'cont_csf'} = $lcont;
        &print_info("LCONT = $lcont\n", $r_par);
      }
    } else { # Target and photoionization runs need the full CI vectors in the outer region.
        $r_par->{'data'}->{'scattering'}->{'cont_csf'} = 0;
    }
  }
  else {
    die "  Error: no output of scatci: file $r_par->{'data'}->{'outputfile'} !\n";
  }

  return 1;
}

# ------------------------ mpi-scatci --------------------------

sub read_mpi_scatci_output {
  my ($r_par) = @_;
  my $istates = 0;
  my $task = $r_par->{'data'}->{'task'};
  my $igeom = $r_par->{'data'}->{'igeom'} - 1; # index of geometry for storage

  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    if ($task eq "target") {
      # Read and print information about target states
      &print_info("  Target states:\n", $r_par);
      if (open(PROPS, "$r_par->{'dirs'}->{'geom'}${bs}fort.$luprop")) {
        while (my $line = <PROPS>) {
          chomp($line);
          if ($line =~ m/^5\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+(\S+)\s+State No.\s+(\d+)\s+(\S+)\s+(\S+)/i) {
            my $spin_sym = "$3.$2";
            &print_info("    No. $1 - $6($3)  $7($2)  $4\n", $r_par);
            if ($nstates{$spin_sym}) { $nstates{$spin_sym}++; }
            else                     { $nstates{$spin_sym} = 1; }
            $r_par->{'data'}->{'target'}->{'states'}->{"$3.$2.$5"} = $4;
            $r_par->{'data'}->{'target'}->{'ordered_states'}->[$istates] = "$spin_sym.$nstates{$spin_sym}";
            if ($igeom > 0) {
              $r_par->{'data'}->{$task}->{'states_all_geom'}->{"$spin_sym.$nstates{$spin_sym}"}->[$igeom] = $energies[$nstates{$spin_sym}-1];
            } else {
              $r_par->{'data'}->{$task}->{'states_all_geom'}->{"$spin_sym.$nstates{$spin_sym}"} = [];
              $r_par->{'data'}->{$task}->{'states_all_geom'}->{"$spin_sym.$nstates{$spin_sym}"}->[$igeom] = $energies[$nstates{$spin_sym}-1];
            }
            $istates++;
          }
        }
        close(PROPS);
      } else {
        die "Error: cannot open output property file \"$r_par->{'dirs'}->{'geom'}${bs}fort.$luprop\"!\n";
      }
      # Perform a dummy DENPROP input configuration with an empty template to set up further target information
      my $dummy_str = "";
      make_denprop_input($r_par, \$dummy_str);
    } else {
      # TODO: scattering
    }
  } else {
    die "  Error: no output of mpi-scatci: file \"$r_par->{'data'}->{'outputfile'}\"!\n";
  }

  close OUTPUT;

  return 1;
}

# ------------------------ hamdiag --------------------------

sub read_hamdiag_output {
  my ($r_par) = @_;
  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
     while (my $line = <OUTPUT>) {
        if ($line =~ m/Only the first\s*(\d+)\s*CI coefficients from each eigenvector will be saved./) {
           $r_par->{'data'}->{'scattering'}->{'cont_csf'} = $1;
           &print_info("LCONT = $1\n", $r_par);
        }
     }
  }
  else {
    die "  Error: no output of hamdiag: file $r_par->{'data'}->{'outputfile'} !\n";
  }

  close OUTPUT;

  return 1;
}

# ------------------------ denprop -------------------------

# From denprop output we have to determine order of target states again
# because it can be different from ordering obtained from scatci outputs if there are degenerate states
#   $r_par->{'data'}->{'target'}->{'states'}
# values of dipole moments are also read from the file borndat, x,y,z components are stored in
#   $r_par->{'data'}->{'target'}->{'dipole_all_geom'}
sub read_denprop_output {
  my ($r_par) = @_;
  my $istates = 0;
  my %nstates = ();
  my %str_idtarg = (); # used below for setting idtarg input for outer
  my $igeom = $r_par->{'data'}->{'igeom'} - 1; # index of geometry for storage

  if ($r_par->{'model'}->{'model'} eq "CHF-A") {
    my @prop_eigs = ();
    my @idtarg = ();
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      foreach my $multiplet (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
        my $spin = $spin_multiplicity{lc($multiplet)};
        for (my $j = 1; $j <= $r_par->{'model'}->{'ntarget_states'}->{$multiplet}[$i]; $j++) {
          #          print "$spin.$i.$j  ",$r_par->{'data'}->{'target'}->{'states'}->{"$spin.$i.$j"}, " \n";
          push @prop_eigs, $r_par->{'data'}->{'target'}->{'states'}->{"$spin.$i.$j"}+0.;
        }
      }
    }

    my @temp1 = sort { $prop_eigs[$a] <=> $prop_eigs[$b] } (0 .. $#prop_eigs);
    my @temp2 = sort { $temp1[$a] <=> $temp1[$b] } (0 .. $#temp1);
    for my $id (@temp2) {
      push @idtarg, $id+1;
    }
    $r_par->{'data'}->{'CHF'}->{'idtarg'} = \@idtarg;
  }

  # Read denprop output and print information about target states
  &print_info("  Target states:\n", $r_par);
  if (open(OUTPUT, "$r_par->{'data'}->{'outputfile'}")) {
    while (my $line = <OUTPUT>) {
      chomp($line);
      if ($line =~ m/^\s*\d+\s+(\d+)\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+([^\s]+)\s+state no.\s+\d+\s+([^\s]+)\s+([^\s]+)/i) {
        my $spin_sym = "$3.$2";
        &print_info("    No. $1 - $5($3)  $6($2)  $4\n", $r_par);
        if ($nstates{$spin_sym}) { $nstates{$spin_sym}++; }
        else                     { $nstates{$spin_sym} = 1; }
        $r_par->{'data'}->{'target'}->{'ordered_states'}->[$istates] = "$spin_sym.$nstates{$spin_sym}";
        $istates++;
        $str_idtarg{"$spin_sym.$nstates{$spin_sym}"} = $istates;
      }
    }
    close(OUTPUT);
  }
  else {
    die "Error: cannot open output of denprop !\n";
  }

  # Read the file borndat, store and print information about the dipole moment
  &print_info("  Dipole moment (in a.u.) of the ground state:\n", $r_par);
  my @order_of_components = ('z', 'x', 'y'); # as in the file borndat
  if (open(BORNDAT, "borndat")) {
    my $line = <BORNDAT>; # the first number not needed
    for (my $i = 0; $i < 3; $i++) {
      $line = <BORNDAT>;
      chomp($line);
      $line =~ s/\s//g; # get rid of spaces
      if ($igeom == 0) {
        $r_par->{'data'}->{'target'}->{'dipole_all_geom'}->{$order_of_components[$i]} = [];
      }
      $r_par->{'data'}->{'target'}->{'dipole_all_geom'}->{$order_of_components[$i]}->[$igeom] = $line;
      &print_info("   $order_of_components[$i] component: $line\n", $r_par);
    }
    close(BORNDAT);
  }
  else {
    &print_info("    Warning: cannot open file borndat !\n", $r_par);
  }

  # Finally set auxiliary array 'id_targ' used later in scattering.outer inputs
  my $r_spinsym_order = $r_par->{'data'}->{'target'}->{'spinsym_order'};
  my $it = 0;
  for (my $i = 0; $i < scalar @{$r_spinsym_order}; $i++) {
    if ($r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$r_spinsym_order->[$i]}) {
      my ($spin, $sym) = split(/\./, $r_spinsym_order->[$i]);
      for ($is = 1; $is <= $r_par->{'data'}->{'target'}->{'used_tgt_states'}->{$r_spinsym_order->[$i]}; $is++) {
        $r_par->{'data'}->{'target'}->{'idtarg'}->[$it] = $str_idtarg{"$spin.$sym.$is"};
        $it++;
      }
    }
  }

  if ($r_par->{'model'}->{'model'} eq "CHF-A") {
    $r_par->{'data'}->{'target'}->{'idtarg'} = $r_par->{'data'}->{'CHF'}->{'idtarg'};
  }
  return 1;
}

# ------------------------- outer --------------------------

# From outer output we can determine position and width of the resonance
#   saved in $r_par->{'data'}->{'scattering'}->{'resonance_position'}
#   and   in $r_par->{'data'}->{'scattering'}->{'resonance_width'}
sub read_reson_output {
  my ($r_par) = @_;
  my $state = $spin_multiplicity{lc($r_par->{'data'}->{'scattering'}->{'spin'})}.'.'.$r_par->{'data'}->{'scattering'}->{'symmetry'};
  my $igeom = $r_par->{'data'}->{'igeom'}; # index of geometry
  my ($position, $width) = (1.0e10, 0.0);

  my $str = "";
  my @positions = ();
  my @widths    = ();
  if (read_file($r_par->{'data'}->{'outputfile'}, \$str)) {
    while ($str =~ s/Fitted resonance parameters\s+Positions\s*\/\s*Ryd\s*(.*?)\s*Widths\s*\/\s*Ryd\s*(.*?)\s*Background//is) {
      $str_positions = $1;          $str_widths = $2;
      $str_positions =~ s/D/E/g;    $str_widths =~ s/D/E/g;
      push(@positions, split(/\s+/, $str_positions));
      push(@widths, split(/\s+/, $str_widths));
    }
    for (my $i = 0; $i < scalar @positions; $i++) {
      if ($positions[$i] < $position && $positions[$i] > 0.0) {
        $position = $positions[$i]; $width = $widths[$i];
      }
    }
    if ($position < 1.0e10) {
      &print_info("  Resonance position: $position Ryd\n", $r_par);
      &print_info("  Resonance width:    $width Ryd\n", $r_par);
    }
    else {
      &print_info("  Warning: there is no resonance in output of reson: file $r_par->{'data'}->{'outputfile'} !\n", $r_par);
      $position = 0.0;
    }
    if ($igeom == 1) {
      $r_par->{'data'}->{'scattering'}->{'resonance_position'}->{"$state"} = [];
      $r_par->{'data'}->{'scattering'}->{'resonance_width'}->{"$state"} = [];
    }
    $r_par->{'data'}->{'scattering'}->{'resonance_position'}->{"$state"}->[$igeom - 1] = 0.5 * $position; # 0.5 is for conversion to a.u.
    $r_par->{'data'}->{'scattering'}->{'resonance_width'}->{"$state"}->[$igeom - 1] = 0.5 * $width;       # 0.5 is for conversion to a.u.
  }
  else {
    if ($r_par->{'run'}->{'only'} !~ /reson/) {
       print "  Error: no output of reson: file $r_par->{'data'}->{'outputfile'} !\n";
       print "  Outer run was not selected so this is not a critical error and I am continuing with the run !\n";
    } else {
       die "  Error: no output of reson: file $r_par->{'data'}->{'outputfile'} !\n";
    }
  }

  return 1;
}

# ============= end of subroutines for OUTPUTS ==============

# ----------------- auxiliary subroutines -------------------

# Check whether settings are consistent
# called from main.pl before loop over geometries
sub check_settings {
  my ($r_par) = @_;
  my $str_error = "Script stopped in check_settings()! Modify settings in model.pl accordingly.\n";
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR

  print_info("Checking user settings given in model.pl ...\n", $r_par);

  # check whether a sum of orbitals in each IR equals a given total number of orbitals
  # if any sum differs from the given total number of ordbitals
  # orbitals in each IRs will be determined automatically according to energy
  for my $which_orbitals ('frozen', 'active', 'virtual') {
    my $norbitals = 0;
    $norbitals += $_ for @{$r_par->{'model'}->{$which_orbitals.'_orbs'}}[0..$nir-1];
    if ($norbitals != $r_par->{'model'}->{'n'.$which_orbitals}) {
      if ($norbitals > 0) { # something wrong
        &print_info("\nWarning: Number of $which_orbitals orbitals ".$r_par->{'model'}->{'n'.$which_orbitals}."\n", $r_par);
        &print_info("         inconsistent with settings in '${which_orbitals}_orbs' !\n", $r_par);
        &print_info("         Number of orbitals in each IRs will be determined automatically !\n\n", $r_par);
      }
      if ($which_orbitals =~ /(frozen|active)/) {
        $r_par->{'data'}->{'orbs_ok'} = 0;
      }
    }
  }

  # in SE and SEP models the target is described by one state obtained at the HF level
  if ($r_par->{'model'}->{'model'} =~ /(SE|SEP)/) {
    if ($r_par->{'model'}->{'ntarget_states_used'} > 1) {
      print "Error in settings: ntarget_states_used is $r_par->{'model'}->{'ntarget_states_used'}\n";
      print "                   which is inconsistent with $r_par->{'model'}->{'model'} model !\n";
      print "                   Only one state (the lowest one) must be used !\n";
      die $str_error;
    }
  }

  #Check if triplet states used with positron scattering - not allowed!
  if ($r_par->{'model'}->{'positron_flag'} == 1 ) {
    my $triplet_casscf_states = 0;
    my $triplet_target_states = 0;
    for (my $i = 1; $i <= $r_par->{'data'}->{'nir'}; $i++) {
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'ncasscf_states'}}) {
        my $spin = $spin_multiplicity{$statespin};
        my $nstates = $r_par->{'model'}->{'ncasscf_states'}->{$statespin}->[$i - 1];
        if ($nstates > 0 and $spin == 3) {
          $triplet_casscf_states += $nstates;
        }
      }
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
        my $spin = $spin_multiplicity{$statespin};
        my $nstates = $r_par->{'model'}->{'ntarget_states'}->{$statespin}->[$i - 1];
        if ($nstates > 0 and $spin == 3) {
          $triplet_target_states += $nstates;
        }
      }
    }
    my $ntriplet_state_tot = $triplet_casscf_states + $triplet_target_states;
    if ($ntriplet_state_tot > 0) {
        print "Error in settings: Number of triplet CASSCF states is $triplet_casscf_states\n";
        print "                   Number of triplet target states is $triplet_target_states\n";
        print "                   Positron_flag is $r_par->{'model'}->{'positron_flag'}\n";
        print "                   triplet states cannot be used for positron scattering!\n";
        die $str_error;
    }
  }
  # CHF calculations only support doublet (N-1)-electron states and singlet N-electron states
  if ($r_par->{'model'}->{'model'} =~ /(^CHF)/) {
    use List::Util qw( sum );

    # check photoionization target and scattering states are supported by CHF model
    if ($r_par->{'run'}->{'photoionization'} == 1){
      for my $multiplet (keys %{$r_par->{'model'}->{'ntarget_states'}}) {
        my @targs = @{$r_par->{'model'}->{'ntarget_states'}->{$multiplet}};
        my $ntargs = sum(@targs);
        if ($ntargs > 0 and $multiplet !~ /doublet/) {
          print "Error in settings: CHF model only supports doublet target states\n";
          print "                   You have requested $ntargs $multiplet state(s)\n";
          print "                   $multiplet = [".join(",",@targs)."]\n";
          die $str_error;
        }
      }
      for my $multiplet (keys %{$r_par->{'model'}->{'scattering_states'}}) {
        my @targs = @{$r_par->{'model'}->{'scattering_states'}->{$multiplet}};
        my $ntargs = sum(@targs);
        if ($ntargs > 0 and $multiplet !~ /singlet/) {
          print "Error in settings: CHF model only supports singlet scattering states\n";
          print "                   You have requested $ntargs $multiplet state(s)\n";
          print "                   $multiplet = [".join(",",@targs)."]\n";
          die $str_error;
        }
      }
    }
    # check electron scattering target and scattering states are supported by CHF model
    else {
      for my $multiplet (keys %{$r_par->{'model'}->{'ntarget_states'}}) {
        my @targs = @{$r_par->{'model'}->{'ntarget_states'}->{$multiplet}};
        my $ntargs = sum(@targs);
        my $error = 0;
        if ($multiplet =~ /singlet/) {
          if (sum(@targs[1 .. $#targs]) > 0) {$error = 1;}
          if (@targs[0] != 1) {$error = 1;}
        }
        else{
          if ($ntargs > 0) {$error = 1;}
        }
        if ($error != 0) {
          print "Error in settings: CHF model only supports a singlet ground state\n";
          print "                   You have requested $ntargs $multiplet state(s)\n";
          print "                   $multiplet = [".join(",",@targs)."]\n";
          die $str_error;
        }
      }
      for my $multiplet (keys %{$r_par->{'model'}->{'scattering_states'}}) {
        my @targs = @{$r_par->{'model'}->{'scattering_states'}->{$multiplet}};
        my $ntargs = sum(@targs);
        if ($ntargs > 0 and $multiplet !~ /doublet/) {
          print "Error in settings: CHF model only supports doublet scattering states\n";
          print "                   You have requested $ntargs $multiplet state(s)\n";
          print "                   $multiplet = [".join(",",@targs)."]\n";
          die $str_error;
        }
      }
    }
  }

  # check whether spins of target states are consistent with number of electrons and required total spin of scattering wf
  # e.g. if we use doublets and quartets for target, we can only do scattering calculations for triplets,
  #      for singlets we can use only doublet states of target,
  #      because one additional electron cannot couple singlets with quartets !
  my $states_ok = 1;
  foreach my $scattering_state_spin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
    my $nstates = 0;
    $nstates += $_ for @{$r_par->{'model'}->{'scattering_states'}->{$scattering_state_spin}}[0..$nir-1];
    if ($nstates > 0) {
      if ($spin_multiplicity{lc($scattering_state_spin)} % 2 != $r_par->{'model'}->{'nelectrons'} % 2) {
        print "Error in settings: required scattering '$scattering_state_spin' states\n";
        print "                   are not valid states of ".($r_par->{'model'}->{'nelectrons'} + 1)." electrons !\n";
        $states_ok = 0;
      }
      foreach my $target_state_spin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
        my $ntarget_states = 0;
        $ntarget_states += $_ for @{$r_par->{'model'}->{'ntarget_states'}->{$target_state_spin}}[0..$nir-1];
        if ($ntarget_states > 0) {
          if (($spin_multiplicity{lc($target_state_spin)} + 1) % 2 != $r_par->{'model'}->{'nelectrons'} % 2) {
            print "Error in settings: required target '$target_state_spin' states\n";
            print "                   are not valid states of ".$r_par->{'model'}->{'nelectrons'}." electrons !\n";
            $states_ok = 0;
          }
          if (abs($spin_multiplicity{lc($scattering_state_spin)} - $spin_multiplicity{lc($target_state_spin)}) > 1) {
            print "Error in settings: required scattering '$scattering_state_spin' states\n";
            print "                   are not coupled to target '$target_state_spin' states !\n";
            $states_ok = 0;
          }
        }
      }
    }
  }
  if ($states_ok == 0) {
    die $str_error;
  }

  # check that the GTO continuum basis radius is compatible with the R-matrix radius
  if ($r_par->{'run'}->{'scattering'} == 1) {
    if ($r_par->{'model'}->{'use_BTO'} == 0 and $r_par->{'model'}->{'use_GTO'} == 1 and
        $r_par->{'model'}->{'radius_GTO'} != $r_par->{'model'}->{'rmatrix_radius'}) {
      print "Warning: Your GTO radius is not equal to the R-matrix radius!\n\n"
    }
  }

  # check whether the number of required target states ('ntarget_states_used') is
  # equal or less than the total number of all states specified in 'ntarget_states' for all spin symmetries
  my $ntarget_states = 0;
  foreach my $target_state_spin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
    $ntarget_states += $_ for @{$r_par->{'model'}->{'ntarget_states'}->{$target_state_spin}}[0..$nir-1];
  }
  if ($r_par->{'model'}->{'ntarget_states_used'} > $ntarget_states) {
    print "Error in settings: ntarget_states_used ($r_par->{'model'}->{'ntarget_states_used'}) is greater than number of states in ntarget_states ($ntarget_states) !\n\n";
    die $str_error;
  }
  elsif ($r_par->{'model'}->{'ntarget_states_used'} < $ntarget_states) {
    print "Warning: ntarget_states_used ($r_par->{'model'}->{'ntarget_states_used'}) is smaller than number of states in ntarget_states ($ntarget_states) !\n";
    print "         If you run codes for more geometries then different target states can be used for different geometries !\n";
    print "         To avoid this select only $r_par->{'model'}->{'ntarget_states_used'} states in ntarget_states.\n\n";
  }

  #Check if MAS/MASSCF is being used with positrons
  if ($r_par->{'model'}->{'use_MAS'} != 0 or $r_par->{'model'}->{'use_MASSCF'} != 0) {
    if ($r_par->{'model'}->{'positron_flag'} == 1) {print "MAS/MASSCF cannot be used with positron scattering at the moment.\n"; die $str_error;}
  }

  print_info("  ... done.\n", $r_par);
  return 1;
}

# Check whether numbers of orbitals are consistent.
# This subroutine is called from read_molpro_output
# after total number of orbitals (basis functions) is determined.
sub check_numbers_of_orbitals {
  my ($r_par) = @_;
  my $str_error = "Script stopped in check_numbers_of_orbitals()! Modify settings in model.pl accordingly.\n";
  my $nir = $r_par->{'data'}->{'nir'};                    # number of IR

  if ($r_par->{'model'}->{'nfrozen'} < 0) {
    $r_par->{'model'}->{'nfrozen'}  = 0;
    &print_info("  Number of frozen orbitals set to 0. Automatic determination of frozen orbitals is not implemented yet!\n", $r_par);
  }
  if ($r_par->{'model'}->{'nvirtual'} < 0) {
    $r_par->{'model'}->{'nvirtual'} = 0;
    &print_info("  Number of virtual orbitals set to 0. Automatic determination of virtual orbitals is not implemented yet!\n", $r_par);
  }
  my $nrequired_orbs = $r_par->{'model'}->{'nfrozen'} + $r_par->{'model'}->{'nactive'} + $r_par->{'model'}->{'nvirtual'};
  if ($r_par->{'model'}->{'nactive'}  < 0) {
    $nrequired_orbs = $r_par->{'model'}->{'nfrozen'} + $r_par->{'model'}->{'nvirtual'};
    $r_par->{'model'}->{'nactive'} = $r_par->{'data'}->{'orbitals'}->{'target_sum'}
                                   - $r_par->{'model'}->{'nfrozen'}
                                   - $r_par->{'model'}->{'nvirtual'};
    &print_info("  Number of active orbitals determined from the basis:\n", $r_par);
    &print_info("    nfrozen, nactive, nvirtual = $r_par->{'model'}->{'nfrozen'}, $r_par->{'model'}->{'nactive'}, $r_par->{'model'}->{'nvirtual'}\n", $r_par);
  }
  if ($nrequired_orbs > $r_par->{'data'}->{'orbitals'}->{'target_sum'}) {
    print "Error in settings: number of frozen, active and virtual orbitals required in model.pl is larger\n";
    print "                   than number of available target molecular orbitals: $nrequired_orbs > $r_par->{'data'}->{'orbitals'}->{'target_sum'} !\n";
    die $str_error;
  }

  return 1;
}

# Run a specific program including making an input before
#   and reading an output afterwards
sub run_code {
  my ($program, $r_par) = @_;
  my $r_dirs = $r_par->{'dirs'};
  my $task = $r_par->{'data'}->{'task'};

  # A replace hyphens with underscores for use in subroutine names
  my $safeprogram = $program;
  $safeprogram =~ s/-/_/;

  # Setting auxiliary variables for input and output files
  $r_par->{'data'}->{'program'} = $program;
  $r_par->{'data'}->{'inputfile'} = "$r_dirs->{'inputs'}$bs$task.$program.inp";
  $r_par->{'data'}->{'outputfile'} = "$r_dirs->{'outputs'}$bs$task.$program.out";
  $r_par->{'data'}->{'errorfile'} = "$r_dirs->{'outputs'}$bs$task.$program.err";

  if ($r_par->{'run'}->{'use_cdenprop'} == 1 and $program eq 'denprop') {
    $r_par->{'data'}->{'inputfile'} = "$r_dirs->{'inputs'}$bs$task.denprop.inp";
    $r_par->{'data'}->{'outputfile'} = "$r_dirs->{'outputs'}$bs$task.denprop.out";
    $r_par->{'data'}->{'errorfile'} = "$r_dirs->{'outputs'}$bs$task.denprop.err";
  }

  # If CONGEN, SCATCI, OUTER or TIME-DELAY are to be run then input and output files contain information about spin and symmetry
  # e.g. target.scatci.triplet.A2, or scattering.congen.doublet.B1.inp
  if ($program =~ /^(congen|scatci|cdenprop|swinterf|rsolve|eigenp|tmatrx|ixsecs|reson|time-delay)$/) {
    $statespin = $r_par->{'data'}->{$task}->{'spin'};
    $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$r_par->{'data'}->{$task}->{'symmetry'}];
    if ($r_par->{'model'}->{'model'} eq "CHF-A" && $task eq "target") {
      my $state_id = sprintf('%03s', $r_par->{'data'}->{'CHF'}->{'current_state'});
      $r_par->{'data'}->{'inputfile'} = "$r_dirs->{'geom'}$bs$r_dirs->{'inputs'}$bs$task.$program.$statespin.$str_ir.$state_id.inp";
      $r_par->{'data'}->{'outputfile'} = "$r_dirs->{'geom'}$bs$r_dirs->{'outputs'}$bs$task.$program.$statespin.$str_ir.$state_id.out";
      $r_par->{'data'}->{'errorfile'} = "$r_dirs->{'geom'}$bs$r_dirs->{'outputs'}$bs$task.$program.$statespin.$str_ir.$state_id.err";
    }
    else{
      $r_par->{'data'}->{'inputfile'} = "$r_dirs->{'geom'}$bs$r_dirs->{'inputs'}$bs$task.$program.$statespin.$str_ir.inp";
      $r_par->{'data'}->{'outputfile'} = "$r_dirs->{'geom'}$bs$r_dirs->{'outputs'}$bs$task.$program.$statespin.$str_ir.out";
      $r_par->{'data'}->{'errorfile'} = "$r_dirs->{'geom'}$bs$r_dirs->{'outputs'}$bs$task.$program.$statespin.$str_ir.err";
    }
  }

  # Run the program if it is to be run and print info ...
  if (&run_program_this_time($r_par) == 1) {

    # If 'use_templates' option is set to 0, skip creating inputs
    # otherwise use template to create input file
    if ($r_par->{'run'}->{'use_templates'} == 1) {
      # Read template for input
      my $str = "";
      if (!&read_file("$r_dirs->{'templates'}$bs$task.$program.inp", \$str)) {
        if(!&read_file("$r_dirs->{'templates'}$bs$program.inp", \$str)) {
          &print_info("Warning: no template file for $program.inp !\n", $r_par);
        }
      }
      # Call make_$safeprogram_input to modify input if exists and then save it
      my $subroutine = "";
      if ($program eq "congen" && exists($r_par->{'model'}{'use_MAS'}) && $r_par->{'model'}{'use_MAS'}){
        $subroutine = "make_${safeprogram}_mas_input";
      } else {
        $subroutine = "make_${safeprogram}_input";
      }
      if ($r_par->{'run'}->{'use_cdenprop'} == 1 and $program eq 'denprop') {&make_denprop_input($r_par, \$str);}
      else{
        if (exists &$subroutine) {
          &$subroutine($r_par, \$str);
        }
      }
      &save_file($r_par->{'data'}->{'inputfile'}, \$str);
    }

    if ($program =~ /^(congen|scatci|cdenprop|swinterf|rsolve|eigenp|tmatrx|ixsecs|reson|time-delay)$/) {
      &print_info("Running $program for $r_par->{'data'}->{$task}->{'spin'} $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$r_par->{'data'}->{$task}->{'symmetry'}] ...\n", $r_par);
    }
    else {
      &print_info("Running $program ...\n", $r_par);
    }
    my $dir = "";
    if    ($program eq "molpro")             { $dir = $r_dirs->{'molpro'}; }
    elsif ($program eq "psi4")               { $dir = $r_dirs->{'psi4'}; }
    elsif ($program eq "molcas")             { $dir = $r_dirs->{'molcas'}; }
    elsif ($program =~ /(swinterf|rsolve|eigenp|tmatrx|ixsecs|reson|time-delay|dipelm|rmt_interface)/) { $dir = $r_dirs->{'bin_out'}; }
    else                                     { $dir = $r_dirs->{'bin_in'}; }

    my $command = "";
    my $input   = $r_par->{'data'}->{'inputfile'};
    my $output  = $r_par->{'data'}->{'outputfile'};
    my $error   = $r_par->{'data'}->{'errorfile'};
    my $redirected_master_output = 0;

    if ($program eq "molpro") {
      if ($sys eq win) {
        $input  =~ s/$rebs/\//g; # On Windows, MOLPRO is running under cygwin and there should be used '/' in a file path
        $output =~ s/$rebs/\//g;
        $error  =~ s/$rebs/\//g;
      };
      $command = "$dir$bs$program$ext_exe --no-xml-output $input --output $output 2> $error";
    } elsif ($program eq "psi4") {
      $command = "$dir$bs$program$ext_exe --input $input --output $output 2> $error";
    } elsif ($program eq "molcas") {
      $command = "$dir$bs$program$ext_exe $input -o $output 2> $error";
    } elsif ($program eq "mpi-scatci" or ($program eq "scatci" and not $r_par->{'run'}->{'mpi_scatci'} eq "")) {
      $command = $r_par->{'run'}->{'mpi_scatci'}  . " ${dir}${bs}mpi-scatci$ext_exe $input 1> $output 2> $error";
    } elsif ($program eq "scatci_integrals") {
      $redirected_master_output = 1;
      &run_sub("make_symlink", $r_par, $input, "inp");
      $command = ($r_par->{'run'}->{'mpi_integrals'} eq "" ? "" : $r_par->{'run'}->{'mpi_integrals'}." ")."$dir$bs$program$ext_exe 2> $error";
    } elsif ($program eq "rsolve" ) {
      if ($r_par->{'run'}->{'use_saved_ramps'} == 1) {
        #Link the saved channel and amplitude data and skip running swintf
        my $channels_file = "..${bs}..${bs}collected_scattering_data${bs}channels${bs}channels.geom$r_par->{'data'}->{'igeom'}.$statespin.$str_ir";
        my $amplitudes_file = "..${bs}..${bs}collected_scattering_data${bs}rmat_amplitudes${bs}ramps.geom$r_par->{'data'}->{'igeom'}.$statespin.$str_ir";
        &print_info("WARNING: using R-matrix amplitudes and channel data from the previous calculation!\n", $r_par);
        &print_info("Channels: $channels_file\n", $r_par);
        &print_info("R-matrix amplitudes: $amplitudes_file\n", $r_par);
        &run_system("$rm_cmd fort.$luchan", $r_par);
        &run_system("$rm_cmd fort.$lurmt", $r_par);
        &run_sub("make_symlink", $r_par, $channels_file, "fort.$luchan");
        &run_sub("make_symlink", $r_par, $amplitudes_file, "fort.$lurmt");
      }
      if (not $r_par->{'run'}->{'mpi_rsolve'} eq "") { #use parallel MPI_RSOLVE
        $redirected_master_output = 1;
        &run_sub("make_symlink", $r_par, $input, "inp");
        $command = $command . " $r_par->{'run'}->{'mpi_rsolve'} ${dir}${bs}mpi_rsolve$ext_exe 2> $error";
      } else { #use serial RSOLVE
        $command = $command . "${dir}${bs}rsolve$ext_exe < $input 1> $output 2> $error";
      }
    } else {
      if ($r_par->{'run'}->{'use_cdenprop'} == 1 and $program eq 'denprop') {
        $command = "${dir}${bs}cdenprop_all$ext_exe < $input 1> $output 2> $error";
      } else {
        $command = "$dir$bs$program$ext_exe < $input 1> $output 2> $error";
      }
    }

    # attempt execution and print the error log on failure
    if (system("$command") != 0) {
      &print_info(("!" x 80)."\n", $r_par);
      &print_info("Execution of \"$command\" failed\nError output follows:\n", $r_par);
      &print_info(("-" x 80)."\n", $r_par);
      if (open(my $ferr, $error)) {
        while ($line = <$ferr>) {
          chomp($line);
          &print_info("$line\n", $r_par);
        }
      }
      &print_info(("!" x 80)."\n", $r_par);
    }

    if ($redirected_master_output) { &run_system("$mv_cmd log_file.0 $output", $r_par); }
    if (-e "inp") { &run_system("$rm_cmd inp", $r_par); }
    if ($program eq "molcas") {
      my $prefix =  lc($r_par->{'model'}->{'molecule'});
      if ($r_par->{'model'}->{'model'} =~ /^SE/ || $r_par->{'data'}->{'scf_ok'} == 0 || $r_par->{'model'}->{'orbitals'} eq "HF") {
        &run_system("$mv_cmd target.molcas.scf.molden ".$prefix.".molden", $r_par);
        &run_system("$rm_cmd *molcas.*", $r_par);
        &run_system("$rm_cmd coord.inp pid rc.global rc.local stdin xmldump", $r_par);
      } else {
        &run_system("$mv_cmd target.molcas.rasscf.molden ".$prefix.".molden", $r_par);
        &run_system("$rm_cmd *molcas.*", $r_par);
        &run_system("$rm_cmd CI_Iterations.txt CleanInput coord.inp pid rc.global rc.local stdin xmldump", $r_par);
      }
    }

    # save error log to output log and remove error log
    if (open(my $fout, ">>", $output) and open(my $ferr, "<", $error)) {
      while ($line = <$ferr>) { print($fout $line); }
      close($fout); close($ferr);
    }
    unlink($error);
  }
  else {
    &print_info("Skipping running $program.\n", $r_par);
  }

  # Call read_$safeprogram_output to gather information necessary for later inputs
  my $subroutine = "read_${safeprogram}_output";
  if ($r_par->{'run'}->{'use_cdenprop'} == 1 and $program eq 'denprop') {my $subroutine = "read_denprop_output";}
  if (exists &$subroutine) {
    &print_info("Reading $program output ...\n", $r_par);
    &$subroutine($r_par);
    &print_info(" ... done.\n", $r_par);
  }

  return 1;
}

sub run_program_this_time {
  my ($r_par) = @_;
  my $task = $r_par->{'data'}->{'task'};
  my $program = $r_par->{'data'}->{'program'};

  # If 'only' in the hash array 'run' is not empty, then only those programs will run which are specified
  my $run_this_time = 0;
  if ($r_par->{'run'}->{'only'} eq "") {
    $run_this_time = 1;
  }
  else {
    if ("$task-$program" =~ /^($r_par->{'run'}->{'only'})$/) {
      $run_this_time = 1;
    }
  }
  return $run_this_time;
}

sub run_system {
  my ($command, $r_par) = @_;
  if (&run_program_this_time($r_par) == 1) {
    &print_info("Running command \"$command\" ...\n", $r_par);
    system($command);
  }
  return 1;
}

sub run_sub {
  my ($sub_name, $r_par, @args) = @_;

  #Enforce running of the routines with suffix _ALWAYS, typically make_symlink
  $forced_run = 0;
  if ($sub_name =~ /_ALWAYS/) {
     $forced_run = 1;
     $sub_name =~ s/_ALWAYS//;
  }
  if (&run_program_this_time($r_par) == 1 or $forced_run == 1) {
    &$sub_name(@args);;
  }
  return 1;
}

sub print_info {
  my ($message, $r_par) = @_;
  if ($r_par->{'run'}->{'print_info'} eq "file" or $r_par->{'run'}->{'print_info'} eq "both") {
    open(LOG, ">>$r_par->{'data'}->{'logfile'}");
    print LOG $message;
    close(LOG);
  }
  if ($r_par->{'run'}->{'print_info'} eq "screen" or $r_par->{'run'}->{'print_info'} eq "both") {
    print $message;
  }
  return 1;
}

# Auxiliary subroutine which adds the current working directory to all user specified directories
sub add_cwd_to_dirs {
  my ($r_dirs) = @_;
  foreach my $dir (keys %$r_dirs) {
    if ($dir =~ /(bin_in|bin_out|basis|templates|output|molpro)/) {
      if ($r_dirs->{$dir} eq ".") {
        $r_dirs->{$dir} = $r_dirs->{'cwd'};
      }
      else {
        if ($r_dirs->{$dir} !~ /^(\/|[A-Z]:|[a-z]:)/) { # if it is not a full path ( /... on linux, X:... on win)
          $r_dirs->{$dir} = "$r_dirs->{'cwd'}$bs$r_dirs->{$dir}";
        }
      }
    }
  }
  return 1;
}

sub replace_in_template {
  my ($r_str, $what, $with) = @_;
  $with =~ s/,\s*$//; # Get rid of ',' at the end
  $$r_str =~ s/>>>$what<<</$with/gi;
  return 1;
}

sub replace_all_in_template {
  my ($r_str, $r_namelist) = @_;

  foreach my $what (keys %{$r_namelist}) {
    &replace_in_template($r_str, $what, $r_namelist->{$what});
  }
  return 1;
}

sub add_range {
  my ($r_array, $from, $to, $step) = @_;
  my $epsilon = 1e-10;       # used to deal with rounding errors
  my $n = scalar @$r_array;  # starting number of array elements
  my $value = $from;         # initial value
  my %present = ();          # auxiliary
  while ($value <= $to + $epsilon) {
    $r_array->[$n] = $value;
    $value += $step;
    $n++;
  }
  @$r_array = &delete_duplicates(@$r_array);
  @$r_array = sort {$a <=> $b} (@$r_array);
  return 1;
}

# -- BEGIN: Added for make_conge_mas_input ----
# Routines to insert arrays directly into the input templates with easier to
# read spacing for pqn and mshl.

sub replace_in_template_with_arrays {
  my ($r_str, $what, $with, $sep) = @_;

  my $with_str = "";
  if (ref($with) eq 'ARRAY') {
    if ($what eq 'MSHL'){
      foreach (@{$with}){
        $with_str .= "  $_$sep   ";
      }
    } elsif ($what eq 'PQN'){
      my @pqns = @{$with};
      while (@pqns) {
        my @pqn = splice @pqns, 0, 3;
        $with_str .=  "$pqn[0]$sep$pqn[1]$sep$pqn[2]$sep "
      }
    } elsif ($what eq 'GASSCF'){
      foreach  (@{$with}){
        $with_str .= join($sep, @{$_});
        $with_str .= "\n";
      }
    } else {
      $with_str = join($sep, @{$with});
    }
  } else {
    $with_str = $with;
  };
  $with_str =~ s/,\s*$//; # Get rid of ',' at the end
  $$r_str =~ s/>>>$what<<</$with_str/gi;
  return 1;
}

sub replace_all_in_template_with_arrays {
  my ($r_str, $r_namelist, $sep) = @_;

  foreach my $what (keys %{$r_namelist}) {
    &replace_in_template_with_arrays($r_str, $what, $r_namelist->{$what}, $sep);
  }
  return 1;
}

# -- END: Added for make_conge_mas_input ----

sub delete_duplicates { my %seen; grep !$seen{$_}++, @_ }

# Gather information about eigenphase sums from files eigenph.doublet.A1, ...
# Only those spin symmetries are gathered, which were required in settings->nscat_states
sub gather_eigenphases {
  my ($prefix, $r_par) = @_;
  my @fhs = ();
  my @lines = ();
  # Open all files to filehandles fh0, ...
  my $fhcounter = 0;
  my $data_description = ();
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] == 1) {
        $fhs[$fhcounter] = "fh$fhcounter";
        my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i];
        open($fhs[$fhcounter], "<", "$prefix.$statespin.$str_ir");
        push @data_description, "$statespin $str_ir";
        $fhcounter++;
      }
    }
  }
  open(OUT, ">", "$r_par->{'dirs'}->{'data'}${bs}$prefix${bs}$prefix.all.geom$r_par->{'data'}->{'igeom'}");
  # Print the first line
  print OUT "#Energy    ";
  for (my $i = 0; $i < $fhcounter; $i++) {
    print OUT "  $data_description[$i]";
  }
  print OUT "\n";
  # Read lines in all files at once and save
  my $in = $fhs[0];
  while ($lines[0] = <$in>) {
    for (my $i = 1; $i < $fhcounter; $i++) {
      my $otherin = $fhs[$i];
      $lines[$i] = <$otherin>;
    }
    if ($lines[0] =~ /^\s*\d\.\d+[ED][-\+]\d+/) {
      for (my $i = 0; $i < $fhcounter; $i++) {
        my ($nothing, $energy, $eigenph) = split(/\s+/, $lines[$i]);
        $energy =~ s/D/E/ig;
        $eigenph =~ s/D/E/ig;
        if ($i == 0) { print OUT "$energy  $eigenph"; }
        else         { print OUT "  $eigenph"; }
      }
      print OUT "\n";
    }
  }
  for (my $i = 0; $i < $fhcounter; $i++) {
    close($fhs[$i]);
  }
  return 1;
}

# Print the energy difference between the states produced by the CASSCF quantum
# chemistry calculation used to produce the natural orbitals and the MAS
# target states.
sub print_mas_analysis{
  my ($r_par) = @_;
  my %target_states = %{$r_par->{'data'}->{'target'}->{'states'}};
  my %qchem_states = %{$r_par->{'data'}->{'qchem'}->{'states'}};

  &print_info("\nMAS Analysis: Comparison of qchem and target states:\n", $r_par);
  &print_info("  Energy difference = target energy - qchem energy\n", $r_par);
  &print_info("   ---------------------------------------------------\n", $r_par);
  &print_info("    State  Spin.Sym  Energy diff.  No. CSFs  No. CSFs\n", $r_par);
  &print_info("     No.                 (Ha)       (qchem)  (target)\n", $r_par);
  &print_info("   ---------------------------------------------------\n", $r_par);

  my @sorted_target = sort {$target_states{$a} <=> $target_states{$b}} keys %target_states;

  my %max_diff = ('state', 0, 'state_idx', 0);
  my $i = 0;
  foreach my $state (@sorted_target){
    my $curr_state_diff = 0;
    if (exists($qchem_states{$state})){
      $curr_state_diff = $target_states{$state} - $qchem_states{$state};
      if (abs($curr_state_diff) > $max_diff{'state'} ) {
        $max_diff{'state'} = $curr_state_diff;
        $max_diff{'state_idx'} = $i + 1;
      }
    }
    $i++;
    $idx = sprintf("%5d", $i);
    $curr_state_diff = sprintf("%.6e", $curr_state_diff);
    my ($spin, $sym) = (split /\./, $state)[0, 1];
    my $nocsf_qchem;
    if ($r_par->{'data'}->{'qchem'}->{'nocsfs'}->{"$spin.$sym"} >= 10**8){
      $nocsf_qchem = sprintf("%.2e", $r_par->{'data'}->{'qchem'}->{'nocsfs'}->{"$spin.$sym"});
    } else {
      $nocsf_qchem = sprintf("%8d", $r_par->{'data'}->{'qchem'}->{'nocsfs'}->{"$spin.$sym"});
    }
    my $nocsf_target;
    if ($r_par->{'data'}->{'MAS'}->{'nocsfs'}->{"$spin.$sym"} >= 10**8){
      $nocsf_target = sprintf("%.2e", $r_par->{'data'}->{'MAS'}->{'nocsfs'}->{"$spin.$sym"});
    } else {
      $nocsf_target = sprintf("%8d", $r_par->{'data'}->{'MAS'}->{'nocsfs'}->{"$spin.$sym"});
    }

    &print_info("   $idx      $spin.$sym    $curr_state_diff  $nocsf_qchem  $nocsf_target\n", $r_par);
  }
  &print_info("   ---------------------------------------------------\n", $r_par);
  # $max_diff{'state'} = sprintf("%.6e", $max_diff{'state'});
  my $md = sprintf("%.6e", $max_diff{'state'});
  &print_info("  Max. energy difference: $md for state $max_diff{'state_idx'}\n\n", $r_par);

  return %max_diff;
}

# Convert files with the cross sections in such a way that gnuplot can plot them
# Also gather information about total cross sections from these files
# Only those spin symmetries are gathered, which were required in settings->scattering_states
sub convert_and_gather_cross_sections {
  my ($prefix, $r_par) = @_;
  my @lines = ();
  my @data = ();
  my $state = 0;
  my $title = "";

  # first convert all files to be readable by gnuplot
  for (my $ir = 0; $ir < $r_par->{'data'}->{'nir'}; $ir++) {
    foreach my $statespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$ir] == 1) {
        my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$ir];
        open(IN, "$prefix.$statespin.$str_ir");
        $state = 0;
        while (my $line = <IN>) {
          chomp($line);
          if ($line =~ /(CROSS SECTIONS.*\s)(\d+)\s*$/) {
            $title = $1;
            my $new_state = $2;
            if ($state < $new_state) {
              if ($state > 0) {
                my $file_name = "$r_par->{'dirs'}->{'data'}${bs}$prefix${bs}$prefix.$statespin.$str_ir.from_initial_state_$state.geom$r_par->{'data'}->{'igeom'}";
                &print_and_add_cross_sections($state, $file_name, "$title$state", \@lines, \@data);
              }
              for (my $i = 0; $i < scalar @lines; $i++) { $lines[$i] = ""; }
            }
            $state = $new_state;
          }
          elsif ($line =~ /^\s*I\s+E/) {
            if ($lines[0] && $lines[0] ne "") {
              $line =~ s/^\s*I\s+E\([^\)]*\)\s\s\s//;
              $lines[0] .= $line;
            }
            else {
              $line =~ s/^\s*I\s+/#   /;
              $lines[0] = $line;
            }
          }
          elsif ($line =~ /^\s*(\d+)\s+\d+\.\d+/) {
            my $i = $1;
            $line =~ s/D/E/gi;
            if ($lines[$i] && $lines[$i] ne "") {
              $line =~ s/^\s*\d+\s+[^\s]+//;
              $lines[$i] .= $line;
            }
            else {
              $line =~ s/^\s*\d+\s+/  /;
              $lines[$i] = $line;
            }
          }
        } # while (my $line = <IN>)
        if ($state > 0) {
          my $file_name = "$r_par->{'dirs'}->{'data'}${bs}$prefix${bs}$prefix.$statespin.$str_ir.from_initial_state_$state.geom$r_par->{'data'}->{'igeom'}";
          &print_and_add_cross_sections($state, $file_name, "$title$state", \@lines, \@data);
        }
        close(IN);
      } # if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$ir] == 1) {
    } # foreach my $statespin
  } # for (my $ir = 0; $ir < $r_par->{'data'}->{'nir'}; $ir++)
  for (my $is = 1; $is <= $state; $is++) {
    my $file_name = "$r_par->{'dirs'}->{'data'}${bs}$prefix${bs}$prefix.total.from_initial_state_$is.geom$r_par->{'data'}->{'igeom'}";
    open(OUT, ">$file_name");
    print OUT "# TOTAL $title$is\n";
    print OUT "$lines[0]";
    for (my $i = 0; $i < scalar @{$data[$is]}; $i++ ) {
      print OUT "  ".join("  ", map(sprintf("%14.6e", $_), @{$data[$is]->[$i]}))."\n";
    }
    close(OUT);
  }
  return 1;
}

sub print_and_add_cross_sections {
  my ($state, $filename, $title, $r_lines, $r_data) = @_;
  open(OUT, ">$filename");
  print OUT "# $title\n";
  for (my $i = 0; $i < scalar @{$r_lines}; $i++ ) {
    print OUT "$r_lines->[$i]\n";
  }
  close(OUT);
  if (!($r_data->[$state])) { $r_data->[$state] = []; }
  for (my $i = 1; $i < scalar @{$r_lines}; $i++) {
    $r_lines->[$i] =~ s/^\s*//;
    $r_lines->[$i] =~ s/\s*$//;
    my @values = split(/\s+/, $r_lines->[$i]);
    if ($r_data->[$state]->[$i]) {
      for (my $j = 1; $j < scalar @values; $j++) {
        $r_data->[$state]->[$i]->[$j] += $values[$j];
      }
    }
    else {
      $r_data->[$state]->[$i] = [];
      for (my $j = 0; $j < scalar @values; $j++) {
        $r_data->[$state]->[$i]->[$j] = $values[$j];
      }
    }
  }
  return 1;
}

# Save target energies from $r_par->{'data'}->{'target'}->{'states_all_geom'} into a specified file
# and create gnuplot file for all states
sub save_target_energies {
  my ($filename, $r_par) = @_;
  my $igeom = 0;

  # units
  my $str_r_unit = "a.u.";
  if    ($r_par->{'model'}->{'r_unit'} == 1) { $str_r_unit = "A"; }
#  my $str_e_unit = "a.u.";
#  if    ($r_par->{'model'}->{'e_unit'} == 2) { $str_e_unit = "eV"; }
#  elsif ($r_par->{'model'}->{'e_unit'} == 1) { $str_e_unit = "Rydberg"; }

  # below it is assumed that molecule is diatomic, change according to your species
  # initialize a gnuplot file
  open(GNUPLOT, ">$r_par->{'dirs'}->{'model'}${bs}$filename.gp");
  print GNUPLOT "r_unit  = 1.0 # to change units\n";
  print GNUPLOT "e_unit  = 1.0 # to change units\n";
  print GNUPLOT "e_shift = 0.0 # to shift energies\n\n";
  print GNUPLOT "set xlabel 'Internuclear distance ($str_r_unit)'\n";
  print GNUPLOT "set ylabel 'Energy (a.u.)'\n";
  print GNUPLOT "set title '".$r_par->{'model'}->{'molecule'}.", ".$r_par->{'model'}->{'model'}."'\n\n";
  print GNUPLOT "plot [:] \\\n";

  # initialize a file with energies
  open(OUTPUT, ">$r_par->{'dirs'}->{'model'}${bs}$filename");
  print OUTPUT "# Target energies / a.u.\n";
  # Labels
  print OUTPUT '#'.$r_par->{'data'}->{'geom_labels'};
  my $n = 1;
  my $gnuplot_str = "";
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
      for (my $j = 0; $j < $r_par->{'model'}->{'ntarget_states'}->{$statespin}->[$i]; $j++) {
        print OUTPUT "  ".sprintf("%16s", "$statespin.$irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i].".($j+1));
        $n++;
        $gnuplot_str .= "  '$filename' u (\$1*r_unit):((\$$n+e_shift)*e_unit) t '".($j+1)." $statespin.$irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i]' w l, \\\n";
      }
    }
  }
  print OUTPUT "\n";
  $gnuplot_str =~ s/, \\\n$/\n/g;
  print GNUPLOT $gnuplot_str;
  close(GNUPLOT);

  # Data
  foreach $r_geom (@{$r_par->{'data'}->{'geometries'}}) {
    print OUTPUT ' '.$r_geom->{'geometry'};
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'ntarget_states'}}) {
        for (my $j = 0; $j < $r_par->{'model'}->{'ntarget_states'}->{$statespin}->[$i]; $j++) {
          my $state = "$spin_multiplicity{$statespin}.$i.".($j+1);
          print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'target'}->{'states_all_geom'}->{$state}->[$igeom]);
        }
      }
    }
    print OUTPUT "\n";
    $igeom++;
  }
  close(OUTPUT);
  return 1;
}

# Save R-matrix energies from $r_par->{'data'}->{'scattering'}->{'states_all_geom'} into a specified file
sub save_rmatrix_energies {
  my ($filename, $r_par, $number) = @_;
  my $igeom = 0;

  open(OUTPUT, ">$r_par->{'dirs'}->{'model'}${bs}$filename");
  print OUTPUT "# R-matrix energies / a.u.\n";
  # Labels
  print OUTPUT '#'.$r_par->{'data'}->{'geom_labels'};
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] > 0) {
        for (my $j = 0; $j < $number; $j++) {
          print OUTPUT "  ".sprintf("%16s", "$statespin.$irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i].".($j+1));
        }
      }
    }
  }
  print OUTPUT "\n";
  # Data
  foreach $r_geom (@{$r_par->{'data'}->{'geometries'}}) {
    print OUTPUT ' '.$r_geom->{'geometry'};
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
        if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] > 0) {
          for (my $j = 0; $j < $number; $j++) {
            my $state = "$spin_multiplicity{$statespin}.$i.".($j+1);
            print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'scattering'}->{'states_all_geom'}->{$state}->[$igeom]);
          }
        }
      }
    }
    print OUTPUT "\n";
    $igeom++;
  }
  close(OUTPUT);
  return 1;
}

# Save target dipole moments from $r_par->{'data'}->{'target'}->{'dipole_all_geom'} into a specified file
sub save_dipole_moments {
  my ($filename, $r_par) = @_;
  my $igeom = 0;

  open(OUTPUT, ">$r_par->{'dirs'}->{'model'}${bs}$filename");
  print OUTPUT "# Dipole moments (critical value is 0.6393 a.u.\n";
  # Labels
  print OUTPUT '#'.$r_par->{'data'}->{'geom_labels'};
  print OUTPUT "  ".sprintf("%16s", "dm_x / a.u.")."  ".sprintf("%16s", "dm_y / a.u.")."  ".sprintf("%16s", "dm_z / a.u.")."  ".sprintf("%16s", "dm / a.u.")."\n";
  # Data
  foreach $r_geom (@{$r_par->{'data'}->{'geometries'}}) {
    my $dmx = $r_par->{'data'}->{'target'}->{'dipole_all_geom'}->{'x'}->[$igeom];
    my $dmy = $r_par->{'data'}->{'target'}->{'dipole_all_geom'}->{'y'}->[$igeom];
    my $dmz = $r_par->{'data'}->{'target'}->{'dipole_all_geom'}->{'z'}->[$igeom];
    print OUTPUT ' '.$r_geom->{'geometry'};
    print OUTPUT "  ".sprintf("%16.9f", $dmx);
    print OUTPUT "  ".sprintf("%16.9f", $dmy);
    print OUTPUT "  ".sprintf("%16.9f", $dmz);
    print OUTPUT "  ".sprintf("%16.9f", sqrt($dmx*$dmx + $dmy*$dmy + $dmz*$dmz))."\n";
    $igeom++;
  }
  close(OUTPUT);
  return 1;
}

# Save positions and widths of the resonance from $r_par->{'data'}->{'scattering'}->{'resonance_position'}
# and $r_par->{'data'}->{'scattering'}->{'resonance_width'}
sub save_resonance_positions_and_widths {
  my ($filename, $r_par) = @_;
  my $igeom = 0;

  open(OUTPUT, ">$r_par->{'dirs'}->{'model'}${bs}$filename");
  print OUTPUT "# A = anion energy / a.u., E = position / a.u., W = width / a.u.\n";
  # Labels
  print OUTPUT '#'.$r_par->{'data'}->{'geom_labels'};
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] > 0) {
        my $state = "$statespin.$irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i]";
        print OUTPUT "  ".sprintf("%16s", "$state.A")."  ".sprintf("%16s", "$state.E")."  ".sprintf("%16s", "$state.W");
      }
    }
  }
  print OUTPUT "\n";
  # Data
  foreach $r_geom (@{$r_par->{'data'}->{'geometries'}}) {
    print OUTPUT ' '.$r_geom->{'geometry'};
    for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
      foreach my $statespin (sort { $spin_multiplicity{$a} <=> $spin_multiplicity{$b} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
        if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] > 0) {
          my $state = "$spin_multiplicity{$statespin}.$i";
          print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'target'}->{'states_all_geom'}->{$r_par->{'data'}->{'target'}->{'ground_state'}}->[$igeom] +
                                              $r_par->{'data'}->{'scattering'}->{'resonance_position'}->{$state}->[$igeom]);
          print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'scattering'}->{'resonance_position'}->{$state}->[$igeom]);
          print OUTPUT "  ".sprintf("%16.9f", $r_par->{'data'}->{'scattering'}->{'resonance_width'}->{$state}->[$igeom]);
        }
      }
    }
    print OUTPUT "\n";
    $igeom++;
  }
  close(OUTPUT);
  return 1;
}

# Make gnuplot files to plot eigenphase sums
sub make_gnuplot_files_for_eigenphase_sums {
  my ($prefix, $r_par) = @_;
  my $igeom = 0;

  # we need a gnuplot file for each spin symmetry
  my $col = 2;
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] == 1) {
        my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i];
        open(GNUPLOT, ">$r_par->{'dirs'}->{'data'}${bs}$prefix.$statespin.$str_ir.gp");
        print GNUPLOT "dir = '$prefix'\n\n";
        print GNUPLOT "set title '".$r_par->{'model'}->{'molecule'}.", $statespin $str_ir'\n\n";
        my $str_unit = "a.u.";
        if    ($r_par->{'model'}->{'e_unit'} == 2) { $str_unit = "eV"; }
        elsif ($r_par->{'model'}->{'e_unit'} == 1) { $str_unit = "Rydberg"; }
        print GNUPLOT "set xlabel 'Energy ($str_unit)\n";
        print GNUPLOT "set ylabel 'Eigenphase sum (rad)\n\n";
        print GNUPLOT "e_unit = 1.0 \# to change energy units\n\n";
        print GNUPLOT "plot [:] \\\n";
        my $igeom = 1;
        my $gnuplot_str = "";
        foreach $r_geom (@{$r_par->{'data'}->{'geometries'}}) {
          $gnuplot_str .= "  dir.'${bs}$prefix.all.geom$igeom' u (\$1*e_unit):$col t '".$r_geom->{'gnuplot_desc'}."' w l, \\\n";
          $igeom++;
        }
        $gnuplot_str =~ s/, \\\n$/\n/g;
        print GNUPLOT $gnuplot_str."\n";
        close(GNUPLOT);
        $col++;
      }
    }
  }
  return 1;
}

# Make gnuplot files to plot total cross sections
sub make_gnuplot_files_for_cross_sections {
  my ($prefix, $r_par) = @_;
  my $igeom = 0;

  # we need a gnuplot file for each spin symmetry
  my $col = 2;
  for (my $i = 0; $i < $r_par->{'data'}->{'nir'}; $i++) {
    foreach my $statespin (sort { $spin_multiplicity{lc($a)} <=> $spin_multiplicity{lc($b)} } keys %{$r_par->{'model'}->{'scattering_states'}}) {
      if ($r_par->{'model'}->{'scattering_states'}->{$statespin}->[$i] == 1) {
        my $str_ir = $irred_repr{$r_par->{'model'}->{'symmetry'}}->[$i];
        open(GNUPLOT, ">$r_par->{'dirs'}->{'data'}${bs}$prefix.$statespin.$str_ir.gp");
        print GNUPLOT "dir = '$prefix'\n\n";
        print GNUPLOT "initial = 1 # initial electronic state\n";
        print GNUPLOT "final   = 1 # final electronic state (0 gives total cross section)\n\n";
        print GNUPLOT "set title '".$r_par->{'model'}->{'molecule'}.", $statespin $str_ir'\n\n";
        my $str_e_unit = "a.u.";
        if    ($r_par->{'model'}->{'e_unit'} == 2) { $str_e_unit = "eV"; }
        elsif ($r_par->{'model'}->{'e_unit'} == 1) { $str_e_unit = "Rydberg"; }
        my $str_x_unit = "a.u.";
        if    ($r_par->{'model'}->{'x_unit'} == 2) { $str_x_unit = "A^2"; }
        print GNUPLOT "set xlabel 'Energy ($str_e_unit)\n";
        print GNUPLOT "set ylabel 'Cross section ($str_x_unit)\n\n";
        print GNUPLOT "e_unit = 1.0 \# to change energy units\n";
        print GNUPLOT "x_unit = 1.0 \# to chenge cross-section units\n\n";
        print GNUPLOT "plot [:] \\\n";
        my $igeom = 1;
        my $gnuplot_str = "";
        foreach $r_geom (@{$r_par->{'data'}->{'geometries'}}) {
          $gnuplot_str .= "  dir.'${bs}$prefix.$statespin.$str_ir.from_initial_state_'.initial.'.geom$igeom'";
          $gnuplot_str .= " u (\$1*e_unit):(column(final+2)*x_unit) t '".$r_geom->{'gnuplot_desc'}."' w l, \\\n";
          $igeom++;
        }
        $gnuplot_str =~ s/, \\\n$/\n/g;
        print GNUPLOT $gnuplot_str."\n";
        close(GNUPLOT);
        $col++;
      }
    }
  }
  return 1;
}

1;

# configure CHF data
sub initialize_CHF_data {
  my (%parameters) = @_;

  my $model = %parameters{'model'};
  my $data = %parameters{'data'};
  my $run = %parameters{'run'};

  my @act_mos;
  my $n = 0;
  for (my $i = 0; $i < $data->{'nir'}; $i++) {
    for (my $j = 1; $j <= $data->{'orbitals'}->{'active'}->[$i]; $j++) {
      $n += 1;
      push(@act_mos, [$n, $i, $j]);
    }
  }

  my @vir_mos;
  $n = 0;
  for (my $i = 0; $i < $data->{'nir'}; $i++) {
    for (my $j = 1; $j <= $data->{'orbitals'}->{'virtual'}->[$i]; $j++) {
      push(@vir_mos, [$n, $i, $j]);
      $n += 1;
    }
  }

  my @chf_states;

  my %state = (
    'id', 0,
    'sym', 0,
    'multiplet', "",
    'config', 0,
    'no_csfs', 0,
  );

  my $na = scalar @act_mos; # number of active orbitals
  my $nv = scalar @vir_mos; # number of virtual orbitals
  my $ne;                   # number of electrons in shell

  if ($run->{'photoionization'} == 1) {

    for my $multiplet (keys %{$model->{'ntarget_states'}}) {
      for (my $i = 0; $i < $data->{'nir'}; $i++) {
        $model->{'ntarget_states'}->{$multiplet}->[$i] = 0;
      }
    }

    # generate SE-like states
    my $target_id = 0;
    for (my $k = $na; $k > 0; $k--){
      $target_id += 1;
      $state{'id'} = $target_id;
      $state{'multiplet'} = "doublet";
      my @config = ();
      for (@act_mos){
        $ne = 2;
        my ($nc, $sym, $j) = @$_;
        if ($nc eq $k){
          $state{'sym'} = $sym;
          $ne -= 1;
          $model->{'ntarget_states'}->{'doublet'}->[$sym] += 1;
        }
        push(@config, "$sym,$j,$ne");
      }
      $state{'config'} = \@config;
      push(@chf_states, {%state});
    }

    # generate SEP-like states
    for (@vir_mos){
      my ($vn, $vsym, $vj) = @$_;
      for (my $i = 1; $i <= $na; $i++){
        for (my $k = $i; $k <= $na; $k++){
          my $targ_sym = 0;
          $target_id += 1;
          $state{'id'} = $target_id;
          $state{'multiplet'} = "doublet";
          my @config = ();
          for (@act_mos){
            $ne = 2;
            my ($nc, $sym, $j) = @$_;
            if ($nc eq $k){ $targ_sym = $group_table[$targ_sym]->[$sym]; $ne-=1;}
            if ($nc eq $i){ $targ_sym = $group_table[$targ_sym]->[$sym]; $ne-=1;}
            push(@config, "$sym,$j,$ne");
          }
          $targ_sym = $group_table[$targ_sym]->[$vsym];
          $state{'sym'} = $targ_sym;
          push(@config, "$vsym,$vj,1");
          $state{'config'} = \@config;
          $model->{'ntarget_states'}->{'doublet'}->[$targ_sym] += 1;
          push(@chf_states, {%state});
        }
      }
    }
  }
  else{ # electron-scattering calculation

    # generate ground state
    # This part assumes ground state is always a singlet (this may not be true)
    # This part assumes ground state is always symmetric e.g. sym=0 (this may not be true)
    my $target_id = 1;
    $state{'id'} = $target_id;
    $state{'multiplet'} = 'singlet';
    $state{'targ_sym'} = 0;
    my @config = ();
    for (@act_mos){
      my ($nc, $sym, $j) = @$_;
      push(@config, "$sym,$j,2");
    }
    $state{'config'} = \@config;
    push(@chf_states, {%state});

    # generate SE-like states
    for my $multiplet (keys %{$model->{'ntarget_states'}}) {
      next if ($model->{'positron_flag'} == 1 && $multiplet eq "triplet");
      $state{'multiplet'} = $multiplet;

      for (@act_mos){
        my ($nact, $sym, $j) = @$_;
        for (@vir_mos){
          my ($vn, $vsym, $vj) = @$_;
          my $targ_sym = 0;
          $target_id += 1;
          $state{'id'} = $target_id;
          $state{'multiplet'} = $multiplet;
          my @config = ();
          for (@act_mos){
            my ($nc, $sym, $j) = @$_;
            my $ne = 2;
            if ($nc == $nact){
              $ne = 1;
            }
            push(@config, "$sym,$j,$ne");
          }
          $targ_sym = $group_table[$vsym]->[$sym];
          $state{'sym'} = $targ_sym;
          push(@config, "$vsym,$vj,1");
          $state{'config'} = \@config;
          $model->{'ntarget_states'}->{$multiplet}->[$targ_sym] += 1;
          push(@chf_states, {%state});
        }
      }
    }
  }


  # sort the target states by spin-space symmetry order
  @chf_states = sort { $a->{'multiplet'} cmp $b->{'multiplet'} } @chf_states;
  @chf_states = sort { $a->{'sym'} cmp $b->{'sym'} } @chf_states;

  # re-label state ids
  $target_id = 0;
  foreach $state (@chf_states) {
    $target_id += 1;
    $state->{'id'} = $target_id;
  }

  $model->{'ntarget_states_used'} = $target_id;
  $data->{'CHF'}->{'states'} = \@chf_states;
  $data->{'CHF'}->{'current_state'} = 0;

  return 1;
}
