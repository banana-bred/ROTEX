dir = 'xsec'

initial = 1 # initial electronic state
final   = 1 # final electronic state (0 gives total cross section)

set title 'H2OX, singlet B1'

set xlabel 'Energy (eV)
set ylabel 'Cross section (a.u.)

e_unit = 1.0 # to change energy units
x_unit = 1.0 # to chenge cross-section units

plot [:] \
  dir.'/xsec.singlet.B1.from_initial_state_'.initial.'.geom1' u ($1*e_unit):(column(final+2)*x_unit) t '0.0' w l

