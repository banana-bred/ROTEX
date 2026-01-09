#let link2(label, dest) = link(dest)[#emph(underline(label))]

// make sure HTML renders equations
#let target = dictionary(std).at("target", default: () => "paged")
#show math.equation: it => context {
  if target() == "html" {
    // For inline equations, wrap in a box so they don't break paragraphs
    show: if it.block { it => it } else { box }
    html.frame(it)
  } else {
    it
  }
}

See #link2("docs", "https://banana-bred.github.io/ROTEX/index.html") for more detail.
Github markdown may have trouble rendering some things.

= ROTational EXcitation of ions by electron impact (ROTEX)
Calculates electron-impact rotational (de-)excitation cross sections for asymmetric top molecules, using
+ the _Coulomb-Born (CB) approximation_ and
+ _Multichannel Quantum Defect Theory (MQDT)_.
One or both can be used to determine electron-impact rotational excitation cross sections.
If both are used, then the so-called Born closure will be used for electron-impact transitions that appear in the MQDT and CB cross sections to account for partial waves that are not included in the MQDT treatment:
#align(center)[
  $σ("total") = σ("MQDT") + σ("TCB") - σ("PCB")$
]
- $σ("TCB")$: CB cross sections for $l=0-∞$
- $σ("PCB")$: CB cross sections for $l=0-l_"max"$
- $σ("MQDT")$: MQDT cross sections for $l=0-l_"max"$
The quantity $l_"max"$ is the largest $l$ value for the partial waves in the scattering calculations used to generate K-matrices used in the MQDT approach.

In principle, this also works for symmetric tops and linear rotors, but those special cases have not been explicitly programmed in.

=== Build dependencies
Neither of these are necessary of course, but it will be much easier to build with one of them.
- (_optional)_ The #link2("Fortran Package Manager (fpm)", "https://github.com/fortran-lang/fpm")
- (_optional_) #link2("GNU Make", "https://www.gnu.org/software/make/")
Currently, `make` gives more flexibility, as the user can choose which optional dependencies to exclude (because they are included by default), but `fpm` easily saves and automatically targets differently compiled versions of the code (e.g., "release" and "debug" profiles).

=== Dependencies
- A LAPACK implementation
- #link2("Bspline-fortran","https://github.com/jacobwilliams/bspline-fortran") for combining cross sections on the same energy grid
- #link2("WignerD","https://github.com/banana-bred/WignerD") for rotating rigid rotor eigenvectors
- (_optional_) #link2("CDMSreader", "https://github.com/banana-bred/CDMSreader") for reading in #link2("CDMS", "https://cdms.astro.uni-koeln.de/") data to use in the Coulomb-Born approximation
- (_optional_) #link2("Forbear", "https://github.com/szaghi/forbear") for progress bars
- (_optional_) #link2("OpenMP", "https://www.openmp.org/") for easy thread parallelization

=== Building with fpm
In the package directory, just run

`fpm build --profile release`

to download dependencies and build the project.
Fpm alway uses all dependencies because (as of writing this) it is not yet capable of using the C preprocessor to handle optional `USE` statements.

=== Building without fpm
In the package directory, just run

`make`

to install the executable to `bin/rotex`.
The dependency repositories are automatically assumed to be included in `external/` (see the Makefile for more details).
To omit the use of `CDMSreader`, `Forbear`, or `OpenMP`, just set the corresponding environmental variable to 0.
For example,

`make clean build USE_CDMSREADER=0 USE_FORBEAR=0 USE_OPENMP=0`

disables those three optional dependencies.
Be sure to change the `FC` variable if `gfortran` is not your compiler of choice (and adjust `FFLAGS` accordingly).
See the `Makefile` for more detail.

== Usage
An executable will be built by fpm which can be run directly, but fpm can also be invoked to run the code, which expects a namelist file to be given via stdin, i.e.

`fpm run --profile release < example/example.mqdtr2k.H2O+.namelist`

or if you built with `make`

`./bin/rotex < example/example.mqdtr2k.H2O+.namelist`

=== Compatibility with different operating systems
This has only been tested on Linux-based systems for now.
In principle, it should work on others but that is not a guarantee.

=== Structure of the input (namelist) files
Namelists are used to read in values to whole groups of variables that will control program execution.
A namelist group name is defined with the ampersand (&) and the namelist group is terminated with with a foward slash (/).
The program contains three namelists:
+ &control_namelist
+ &kmat_namelist
+ &coulomb_namelist
The namelist variables are declared in `src/rotex__reading.f`
The order of the namelist groups in the input file is unimportant.
Example namelists that should produce successful runes are given in the `example` directory.

==== #underline[&control_namelist]
The first and main namelist group is always used: *&control_namelist*.
It controls the main flow of the program and contains the following variables.

```
&control_namelist
  output_directory
    !! Output directory for the code
  Nmin =
    !! Smallest N (rotational quantum number) for excitation
  Nmax =
    !! Largest N (rotational quantum number) for excitation
  use_kmat =
    !! Use K-matrices to get cross sections ? If .true., use KMAT_NAMELIST variables
  use_CB =
    !! Use the Culomb-Born approximation to get cross sections ? If .true., use COULOMB_NAMELIST variables
  rotor_kind =
    !! The rotor kind of the target
    !!   "a": asymmetric top
    !!   "s": symmetric top (not implemented yet)
    !!   "l": symmetric top (not implemented yet)
  target_charge =
    !! The target charge (integer)
  spin_isomer_kind =
    !! The spin isomer kind to enforce for the molecule.
    !!   0*: none (spin symmetry may be self-enforced by the underlying calculations and dipole vector)
    !!   2 : Ka+Kc parity (odd or even)
  zaxis =
    !! The inertial axis along which the molecular z-axis is fixed. Values
    !! must be either "A", "B", or "C" (case independent). It is up to the user to ensure that
    !! the x, y, z axes are properly aligned with the A, B, C axes in the K-matrix calculations
    !! (if used) and the for the dipole moment vector. The following are the three possible cases:
    !!   ZAXIS | X   Y   Z
    !!   ------|----------
    !!   "A"   | B   C   A
    !!   "B"   | C   A   B
    !!   "C"   | A   B   C
  ABC =
    !! The A, B, and C rotationa constants (cm⁻¹) of the target molecule, such that
    !! A = ABC(1)
    !! B = ABC(2)
    !! C = ABC(3)
  xs_zero_threshold =
    !! Cross sections (cm²) that are below this threshold (for both excitation AND de-excitation) will
    !! be excluded from the list of transitions in each method (KMAT or CB). This is essentially a
    !! filter to reduce the number of files that are generated from K-matrix cross sections that have
    !! negligible magnitude. DEFAULT: 0.0, i.e. get the first matrix that we see
  ! -- 4th order centrifugal distortion
  add_cd4 =
    !! Whether to add 4th order centrifugal distortion. DEFAULT: .false.
  add_cd6 =
    !! Whether to add 6th order centrifugal distortion. DEFAULT: .false.

  ! -- 4th order centrifugal distortion constants (all DEFAULT 0.0)
  DN     = !! ΔN (ΔJ)
  DNK    = !! ΔNK (ΔJK)
  DK     = !! ΔK
  deltan = !! δN (δJ)
  deltak = !! δK
  ! -- 6th order centrifugal distortion constnats (all DEFAULT 0.0)
  HN    = !! HN (HJ)
  HNK   = !! HNK (HJK)
  HKN   = !! HKN (HKJ)
  HK    = !! HK
  etan  = !! ηN (ηJ)
  etank = !! ηNK (ηJK)
  etak  = !! ηK
/
```

==== #underline[&coulomb_namelist]
Next is *&coulomb_namelist*, which controls the execution of the code with respect to the CB approximation.
This can be run with only the data that is included in the namelist group, but can also use data from the CDMS to get Einstein A coefficients, which are used to determine transition dipole moments between rotational state, but can also use data from the CDMS to get Einstein A coefficients, which are used to determine transition dipole moments between rotational states.

```
! =============================================================================================
&coulomb_namelist
! =============================================================================================
  use_CDMS_einstA =
    !! Whether to use CDMS data to get Einstein A coefficients, which will be used
    !! in determining cross sections in the Coulomb-Born approximation (if that transition
    !! is available)
  CDMS_file =
    !! The file to use for reading the CDMS data if USE_CDMS_EINSTA is true.
  only_einsta =
    !! Whether to ONLY calculate Einstein A coefficients and ignore cross section calculations
    !! For the Coulomb-Born approximation. [DEFAULT .false.]
  do_dipole =
    !! Whether to use the Coulomb-Born approximation for the DIPOLE term of the
    !! interaction potential expansion over the multipoles [DEFAULT .true.]
  do_quadrupole =
    !! Whether to use the Coulomb-Born approximation for the QUADRUPOLE term of the
    !! interaction potential expansion over the multipoles [DEFAULT .false.]
  cartesian_dipole_moments =
    !! Array of dipole moments (Debye) in the order μx, μy, μz such that
    !!   μx = cartesian_dipole_moments(1)
    !!   μy = cartesian_dipole_moments(2)
    !!   μz = cartesian_dipole_moments(3)
  eta_thresh = 50 ! au
    !! The largest value of η' allowed for evaluating
    !! the hypergeometric functions ₂F₁(a,b;c;z)
  Ef =
    !! The largest energy (eV) above the lower state threshold t oconsider for performing
    !! cross section calculations. This is the endpoint of the electron energy grid for CB calculations
  nE =
    !! The number of energies for which to calculate cross sections in the CB approx
  do_xtrap =
    !! Whether to extrapolate the CB (de-)excitation cross sections to the excitation threshold
    !! assuming a 1/E dependence (this way, we don't need to calculate ₂F₁(a,b;c;z) with huge
    !! a,b,c)
  Ei_xtrap =
    !! The smallest energy above the excitation (upper state) threshold to which we extrapolate
    !! CB cross sections. This should not be 0.0. (units: eV)
  nE_xtrap =
    !! The number of extrapolation energies down to the excitation threshold. The total number
    !! of CB cross sections per transition will be approximately nE + nE_xtrap.
  analytic_total_cb =
    !! Whether to use the analytic formula to determine the Total CB cross sections.
    !!   .true.:  use the analytic formula
    !!   .false.: take the sum from l=0 to lmax_total and consider that as the "Total" converged
    !!     result
  lmax_partial =
    !! The partial waves to consider for the Partial CB cross sections.
  lmax_total =
    !! The partial waves to consider for the Total CB cross section if analytic_total_cb is .false.
    !! For a particular multipole expansion (quadrupole is not used, so this only needs to be specified
    !! if analytic_total_cb(1) is .false.)
/
```

==== #underline[&kmat_namelist]
Next is the *&kmat_namelist* <kmat_namelist> group.
It controls behavior of the code with respect to the MQDT implementation (K-matrices from scattering calculations -> rotational frame transformation -> closed-channel elimination -> cross sections).
This requires the output of electron-scattering calculations, which is detailed below.

```
! =============================================================================================
&kmat_namelist
! =============================================================================================
  kmat_dir =
    !! The directory containing the K-matrices that the code will read
  channels_dir =
    !! The directory containing the channel files that the code will read (unused if
    !! KMAT_OUTPUT_TYPE is "mqdtr2k")
  lmax_kmat =
    !! The max value of the partial wave orbital angular momentum quantum number l in the K-matrix files
  point_group =
    !! The Abelian point group for the K-matrix files:
    !!   untested: C1, C2, Ci, C2h, D2, D2h
    !!   tested:   C2v, Cs
  spinmults =
    !! Array of spin multiplicities 2S+1 to consider. K-matrix cross sections will be averaged
    !! over spin multiplicities if there are more than one
  kmat_output_type =
    !! The K-matrix format/type that we will read, for use by this code. Two possible values:
    !!   ukrmol+: the default output of the UKRmol+ codes
    !!   mqdtr2k: explained in the writeup
  real_spherical_harmonics =
    !! Whether the spherical harmonics basis used to generate the K-matrices are real-valued
    !!   .true.: real-valued Ylm (will be transformed to complex-valued Ylλ for frame transformation)
    !!   .false.: complex-valued Ylλ (nothing to do)
  num_egrid_segs =
    !! The number of segments in the total energy grid. Must be ≥ 1
  num_egrid =
    !! Array of number of energies per grid segment (length: NUM_EGRID_SEGS)
  egrid_segs =
    !! Array of the bounds (non-degenerate) of the TOTAL energy grid segments
    !! (length: NUM_EGRID_SEGS + 1, units: eV)
  egrid_spacing =
    !! How each energy grid segment is spaced. 'lin' for linear, 'log' for logarithmic
  kmat_energy_closest =
    !! The target energy (eV) at which to evaluate the K-matrix. In both cases of KMAT_OUTPUT_TYPE
    !! (which for this code should really be kmat_INPUT_type, but oh well), there may be
    !! several evaluation energies for the K-matrix, which means that one K-matrix exists per
    !! evaluation energy. The code tries to pick the first match that is closest to this energy
  kmat_energy_units_override =
    !! Force the code to interpret the K-matrix evaluation energy units as something other than
    !! the default for KMAT_OUTPUT_TYPE. Choices:
    !!   "h": hartree atomic units
    !!   "e": eV
    !!   "r": Rydberg
  channel_energy_units_override =
    !! Force the code to interpret the channel energy units as something other than
    !! the default for KMAT_OUTPUT_TYPE. Choices:
    !!   "h": hartree atomic units
    !!   "e": eV
    !!   "r": Rydberg
/
```

=== Structure of the K-matrix data
If using K-matrices for the MQDT approach, there are two critical options to pick between:

==== kmat_output_type: "ukrmol+"
This format uses data taken directly from the output of the UKRmol+ scattering suite @ukrmolx, run with the UKRmol scripts @ukrmolscripts.
It requires *channel* files and *K-matrix* files.
The channel files are expected to live in the directory determined by `channels_dir` and have the naming format

`channels.geom1.<<SPINMULT>>.<<IRREP>>`

because the code currently only expects one geometry for the target molecule (no vibration).
The K-matrices are expected to live in the directory specified by `kmat_dir` and have a similar naming format:

`K-matrix.geom1.<<SPINMULT>>.<<IRREP>>`

In both cases, `<<SPINMULT>>` is one of `singlet`, `doublet`.. and `<<IRREP>>` is one of the irreducible representations of the point group in which the calculations were performed.

- $"C"_s$ : `Ap` and `App` #h(4em) $"C"_(2v)$ : `A1`, `A2`, `B1`, and `B2`; #h(4em)...

==== kmat_output_type: "mqdtr2k"
This format is nonstandard and has very specific output with specific meaning.
The K-matrix files contain the channel data, so only `kmat_dir` needs to be defined.
It shall contain the K-matrix files with the naming pattern

`<<SPINMULT>><<IRREP>>.kmat`

where `<<IRREP>>` is the same as above but `<<SPINMULT>>` is the actual spin multiplicity and not the name, i.e.

- `1` (*not* `singlet`)
- `2` (*not* `doublet`)
- etc.

The line-by line structure is given as follows (values separated by some whitespace, unformatted read (`*`) can be used):

- header (no information to read)
- `nchan`
- `ichan   ielec   l   m   Echan   iq` $<-$ [nchan of these lines]
- `Ekmat   K(1,1)   K(2,1)   K(2,2)  ...  K(nchan,nchan)` [one line for each `Ekmat`]

where

#show table.cell: it => {
  if it.y == 0 {
    set text(weight: 700)
    strong(it)
  } else { it }
}
#table(
  columns: 2,
  [Variable], [Definition],
 [`nchan`], [number of electronic channels in this spin multiplicity and irrep] ,
 [`ichan`], [channel index] ,
 [`ielec`], [electronic state index for that channel] ,
 [`l`], [orbital angular momentum quantum number for that channel] ,
 [`m`], [projection of l on the molecular z-axis] ,
 [`Echan`], [channel energy (units given by default values or user choice; see #link2("kmat_namelist", <kmat_namelist>) above)] ,
 [`iq`], [
   The normalization scheme of the Coulomb $f$ and $g$ functions:
 - iq=4, the typical energy-normalized Coulomb functions
 - iq=0, the $f^η$, $g^η$ functions @hvizdos2023bound] ,
)

#bibliography("refs.bib", title: "References")
