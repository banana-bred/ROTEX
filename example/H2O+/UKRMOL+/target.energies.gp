r_unit  = 1.0 # to change units
e_unit  = 1.0 # to change units
e_shift = 0.0 # to shift energies

set xlabel 'Internuclear distance (A)'
set ylabel 'Energy (a.u.)'
set title 'H2OX, CAS'

plot [:] \
  'target.energies' u ($1*r_unit):(($2+e_shift)*e_unit) t '1 doublet.B1' w l
