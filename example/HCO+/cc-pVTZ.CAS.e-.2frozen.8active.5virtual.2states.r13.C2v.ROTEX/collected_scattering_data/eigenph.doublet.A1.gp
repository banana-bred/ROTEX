dir = 'eigenph'

set title 'HCOX, doublet A1'

set xlabel 'Energy (eV)
set ylabel 'Eigenphase sum (rad)

e_unit = 1.0 # to change energy units

plot [:] \
  dir.'/eigenph.all.geom1' u ($1*e_unit):2 t '0.0' w l

