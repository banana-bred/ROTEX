# Geometries can be specified either manually in the array 'geometries'
# by copying what is between start copy and end copy for each geometry
# or automatically as shown below the hash array %geometry
#
# Note that 'length_unit' should be set to the actual length unit used for values in 'atoms' array,
# later (in &generate_geometries) $model{'r_unit'} is set also to this unit for consistency reasons
#
# Note that there's a restriction on the orientation of the molecules in the scripts (although there isn't one in the codes themselves)
# Your molecular target should be oriented as follows for the corresponding point groups ('any orientation' still requires the axes to
# be along the cartesian axes):
# D2h: any orientation
# C2v: C2 axis along Z
# C2h: C2 axis along Z
# D2:  any orientation
# C2:  C2 axis along Z
# Cs:  Molecule in the YZ plane
# Ci:  any orientation
# C1:  any orientation

# -- equilibrium positions in Center of Charge frame
my $HY0 = 0.821231768657;
my $HZ0 = 0.003519205596;
my $OZ0 = 0.591463405455;
my $OHZ0 = $OZ0 - $HZ0;
my $O1X0 =  0.0; my $O1Y0 =  0.0;  my $O1Z0 =  $OZ0;
my $H2X0 =  0.0; my $H2Y0 =  $HY0; my $H2Z0 =  $HZ0;
my $H3X0 =  0.0; my $H3Y0 = -$HY0; my $H3Z0 =  $HZ0;

%geometry = (

  'suffix', ".EQ.CoM.lmax5",   # string added to model directory to distinguish runs
                              # with different geometry settings

  'geometry_labels',  "   R1       R2       Theta    ",         # labels used on the first line of output files
                                          # it should correspond to numbers given in 'geometries'->'description'
  'correct_cm',   1,                      # correct the center of mass to be at the origin
  'length_unit',  1,                      # 0 - atomic units, 1 - Angstroms

  'geometries', [
    { 'description',  "0.0",
      'gnuplot_desc', "0.0", # used in gnuplot files for keys
       # specify ALL atoms (even redundant with respect to symmetry elements)
       'atoms', [ [ "O", $O1X0,     $O1Y0,     $O1Z0],
                  [ "H", $H2X0,     $H2Y0,     $H2Z0],
                  [ "H", $H3X0,     $H3Y0,     $H3Z0] ]
     },
   ],

  # the following option can be used e.g. to continue an interupted run
  # all geometries are generated but codes will run only for specified geometries
  'start_at_geometry', 1,                # codes will run only for geometries with index >= than this number
  'stop_at_geometry',  0,               # this can be used to stop at certain geometry
                                         # if zero then codes will run for all geometries
);
