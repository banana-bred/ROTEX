! ================================================================================================================================ !
module rotex__writing
  !! Procedures for writing data to disk

  use rotex__kinds,   only: dp
  use rotex__globals, only: G
  use rotex__system,  only: stderr, die, warn, mkdir

  implicit none (type, external)

  private

  public :: write_lifetimes_to_file
  public :: write_channels_to_file
  public :: write_CB_xs_to_file
  public :: write_smat_xs_to_file
  public :: write_total_xs_to_file
  public :: write_Smat_J_elems_to_file
  public :: write_elec_mat_elems_to_file

  character(*), parameter :: ENERGY_XS_WRITE_FMT = '(2X, 2E30.20)'

  interface write_xs_header
    module procedure :: write_xs_header_asymtop
    module procedure :: write_xs_header_symtop
    module procedure :: write_xs_header_linear
  end interface write_xs_header

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine write_lifetimes_to_file(N_states)
    !! Writes the states involved in the excitation and their lifetimes

    use rotex__types,     only: N_states_type
    use rotex__constants, only: au2ev, au2sec

    implicit none (type, external)

    type(N_states_type), intent(in) :: N_states(:)

    integer :: funit
    integer :: inlo, i_tau
    integer :: N, Ka, Kc
    integer :: nc
    real(dp) :: E, EinstA, lifetime
    character(:), allocatable :: filename
    character(:), allocatable :: fmt
    character(15) :: lifetime_char

    filename = G%OUTPUT_DIRECTORY // "lifetimes.dat"
    open(newunit = funit, file = filename)

    nc = maxval(N_states(:) % N) + 1

    ! -- write file header
    select case(G%ROTOR_KIND)
    case("a","A")
      write(funit, '("# ", 3A6, 2A16)') "N", "Ka", "Kc", "energy (meV)", "lifetime (s)"
    case("s","S")
      write(funit, '("# ", 2A6, 2A16)') "N", "K", "energy (meV)", "lifetime (s)"
    case("l","L")
      write(funit, '("# ", A6,  2A16)') "N", "energy (meV)", "lifetime (s)"
    end select

    do inlo = 1, size(N_states, 1)

      N = N_states(inlo) % N
      if(N .lt. G%NMIN) cycle
      if(N .gt. G%NMAX) cycle

      do i_tau = 1, 2*N+1

        Ka     = N_states(inlo) % Ka(i_tau)
        Kc     = N_states(inlo) % Kc(i_tau)
        E      = N_states(inlo) % eigenH % eigvals(i_tau) * au2ev * 1000
        EinstA = N_states(inlo) % einstA(i_tau)

        if(EinstA .eq. 0) then
          write(lifetime_char, '(A15)') "inf"
        else
          lifetime = 1/EinstA * au2sec
          write(lifetime_char, '(F15.6)') lifetime
          if(lifetime .lt. 1e-2 .OR. lifetime .gt. 9.9e5 ) write(lifetime_char, '(E15.6)') lifetime
        endif

        select case(G%ROTOR_KIND)
        case("a","A")
          write(funit, '(2X, 3I6)', advance = "no") N, Ka, Kc
        case("s","S")
          select case(G%SYMAXIS)
          case("a","A")
            write(funit, '(2X, 2I6)', advance = "no") N, Ka
          case("c","C")
            write(funit, '(2X, 2I6)', advance = "no") N, Kc
          case("b","B")
            write(funit, '(2X, I6, A6)', advance = "no") N, "X"
          end select
        case("l","L")
          write(funit, '(2X, I6)', advance = "no") N
        end select
        fmt    = '(X,F15.6)'
        if(E .ne. 0 .AND. E .lt. 1e-2) fmt  = '(X,E15.6)'
        write(funit, fmt, advance = "no") E
        write(funit, '(X, A)') lifetime_char

      enddo

    enddo

    close(funit)

  end subroutine write_lifetimes_to_file

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine write_CB_xs_to_file( &
      prefix                             &
    , output_directory                   &
    , E_el                                &
    , xs                                 &
    , init                               &
    , fin                                &
    , lmax                               &
    , xs_type                            &
    )
    !! Writes a Coulomb-Born cross section to a file whos name and file header
    !! carry information about the state symmetry
    use rotex__types,      only: N_states_type, asymtop_rot_channel_type
    use rotex__characters, only: add_trailing,  sub, sup
    use rotex__constants,  only: au2eV, au2cm
    use rotex__functions,  only: logrange

    implicit none (type, external)

    character(*),                      intent(in) :: prefix
      !! filename prefix
    character(*),                      intent(in) :: output_directory
      !! The directory in which output files are placed
    real(dp),                          intent(in) :: E_el(:)
      !! the scattering eneries in au
    real(dp),                          intent(in) :: xs(:)
      !! The excitation cross sections
    type(asymtop_rot_channel_type), intent(in) :: init, fin
      !! Initial and final states for this transition
    character(*),                      intent(in) :: lmax
      !! The max value of l for the CB cross sections ("inf" if total)
    character(*), intent(in), optional :: xs_type
      !! The kind of cross section that this is

    integer :: ie, iemin
    integer :: ni, nf, kai, kci, kaf, kcf, ne
    integer :: funit
    real(dp) :: Ei, Ef, dE
    character(:), allocatable :: state_name1
    character(:), allocatable :: state_name2
    character(:), allocatable :: filename
    character(:), allocatable :: prefix_local
    character(:), allocatable :: xs_type_

    prefix_local = trim(prefix)
    xs_type_ = trim(prefix) ; if(present(xs_type)) xs_type_ = xs_type
    call add_trailing(prefix_local, ".")

    ne = size(E_el, 1)

    ! -- filenames
    ni  = init % n
    nf  = fin  % n
    kai = init % Ka
    kci = init % Kc
    kaf = fin  % Ka
    kcf = fin  % Kc
    allocate(character(10) :: state_name1)
    allocate(character(10) :: state_name2)
    write(state_name1, '(I0, "_", I0, "_", I0)') ni,  kai,  kci
    write(state_name2, '(I0, "_", I0, "_", I0)') nf,  kaf,  kcf
    state_name1 = trim(state_name1)
    state_name2 = trim(state_name2)

    ! -- get lower and upper state energies
    Ei = init % E
    Ef = fin % E
    dE  = Ef - Ei

    ! -- write cross sections
    filename  = output_directory // prefix_local // state_name1 // "." // state_name2 // ".dat"
    open(newunit = funit,   file = filename)
    call write_xs_header(funit,  ni, kai, kci, nf, kaf, kcf, xs_type_, lmax)
    iemin = findloc(xs .gt. 0, .true., 1)
    do ie = iemin, ne
      write(funit, ENERGY_XS_WRITE_FMT) E_el(ie) * au2eV, xs(ie) * au2cm * au2cm
    enddo
    close(funit)

  end subroutine write_CB_xs_to_file

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine write_smat_xs_to_file( &
      prefix                               &
    , output_directory                     &
    , egrid_total                          &
    , transition                           &
    , exxs                                 &
    , dexxs                                &
    , lmax                                 &
    )
    !! Writes an S-matrix (+CB) cross section to a file whos name and file header
    !! carry information about the state symmetry
    use rotex__types,       only: rvector_type, asymtop_rot_transition_type,  asymtop_rot_channel_type
    use rotex__arrays,      only: size_check
    use rotex__characters,  only: add_trailing, sub, sup, i2c => int2char, state_label, lower
    use rotex__constants,   only: au2eV, au2cm, pi
    use rotex__functions,   only: logrange
    use rotex__channel_ops, only: assert_channel_validity

    implicit none (type, external)

    character(*), intent(in) :: prefix
      !! filename prefix
    character(*), intent(in) :: output_directory
      !! The directory in which output files are placed
    real(dp), intent(in) :: egrid_total(:)
      !! The grid of total energies in au
    type(asymtop_rot_transition_type), intent(in) :: transition
      !! The transition (excitation pair) to be consdered
    real(dp),  intent(in) :: exxs(:)
      !! Array of excitation cross sections for this/all spin multiplicities
    real(dp),  intent(in) :: dexxs(:)
      !! Array of de-excitation cross sections for all spin multiplicities
    integer, intent(in) :: lmax
      !! The max value of l for the K-matrices

    integer :: ie, iemin
    integer :: ne
    integer :: nlo, nup, kalo, kclo, kaup, kcup, ksymup, ksymlo, rchar
    integer :: funit_ex, funit_dex
    real(dp) :: Elo, Eup, Eel_ex, Eel_dex, sigmaup, sigmadown
    character(:), allocatable :: state_name1, state_name2
    character(:), allocatable :: filename_ex, filename_dex
    character(:), allocatable :: prefix_local
    type(asymtop_rot_channel_type) :: lo, up

    prefix_local = trim(prefix)
    call add_trailing(prefix_local, ".")

    ne = size(egrid_total, 1)
    call size_check(exxs,  ne, "EXXS")
    call size_check(dexxs, ne, "DEXXS")

    ! -- make sure the output directory exists before trying to create files therein
    call mkdir(output_directory)

    lo = transition % lo
    call assert_channel_validity(lo, G%ROTOR_KIND, "write_smat_xs_to_file (lo)")
    up = transition % up
    call assert_channel_validity(up, G%ROTOR_KIND, "write_smat_xs_to_file (up)")

    ! -- filenames
    nlo   = lo % n
    nup   = up % n
    kalo  = lo % ka
    kclo  = lo % kc
    kaup  = up % ka
    kcup  = up % kc
    rchar = lo % rchar

    ! -- get upper and lower channel energy
    Elo = lo % E
    Eup = up % E

    select case(G%ROTOR_KIND)
    case("a","A")
      state_name1 = state_label(nlo, kalo, kclo, sepstr="_")
      state_name2 = state_label(nup, kaup, kcup, sepstr="_")
    case("s","S")
      if(lo % rchar .ne. up % rchar .AND. (G%SYMTOP_REDUCE_PROJECTION .eqv. .false.)) then
        write(stderr, *) "lo: ", lo
        write(stderr, *) "up: ", up
        call die("Different rchars detected in this transition. Something has gone wrong !")
      endif
      ksymlo = merge(Kalo, Kclo, lower(G%SYMAXIS) .eq. "a")
      ksymup = merge(Kaup, Kcup, lower(G%SYMAXIS) .eq. "a")
      state_name1 = state_label(nlo, ksymlo, sepstr="_")
      state_name2 = state_label(nup, ksymup, sepstr="_")
    case("l","L")
      state_name1 = state_label(nlo)
      state_name2 = state_label(nup)
    end select
    ! -- write (de-)excitation data
    filename_ex  = output_directory // prefix_local // state_name1 // "." // state_name2 // ".dat"
    filename_dex = output_directory // prefix_local // state_name2 // "." // state_name1 // ".dat"

    iemin = findloc(exxs .gt. 0.0_dp, .true., 1)

    open(newunit = funit_ex,  file = filename_ex)
    open(newunit = funit_dex, file = filename_dex)

    select case(G%ROTOR_KIND)
    case("a","A")
      call write_xs_header(funit_ex,  Nlo, Kalo, Kclo, Nup, Kaup, Kcup, "S-matrix", i2c(lmax))
      call write_xs_header(funit_dex, Nup, Kaup, Kcup, Nlo, Kalo, Kclo, "S-matrix", i2c(lmax))
    case("s", "S")
      call write_xs_header(funit_ex,  Nlo, Ksymlo, Nup, Ksymup, rchar, "S-matrix", i2c(lmax))
      call write_xs_header(funit_dex, Nup, Ksymup, Nlo, Ksymlo, rchar, "S-matrix", i2c(lmax))
    case("l","L")
      call write_xs_header(funit_ex,  Nlo, Nup, "S-matrix", i2c(lmax))
      call write_xs_header(funit_dex, Nup, Nlo, "S-matrix", i2c(lmax))
    end select

    do ie = iemin, ne
      sigmaup   = exxs(ie)
      sigmadown = dexxs(ie)
      Eel_ex    = egrid_total(ie) - Elo
      Eel_dex   = egrid_total(ie) - Eup
      write(funit_ex,  ENERGY_XS_WRITE_FMT) Eel_ex  * au2eV, sigmaup   * au2cm*au2cm
      write(funit_dex, ENERGY_XS_WRITE_FMT) Eel_dex * au2eV, sigmadown * au2cm*au2cm
    enddo

    close(funit_ex)
    close(funit_dex)

  end subroutine write_smat_xs_to_file

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine write_total_xs_to_file( &
      prefix                                &
    , output_directory                      &
    , Eel_ex                                &
    , Eel_dex                               &
    , transition                            &
    , xs_xcite                              &
    , xs_dxcite                             &
    , lmax                                  &
    )
    !! Write the total cross-sections (S-matrix + CB correction) to disk for a single transition

    use rotex__types,      only: asymtop_rot_transition_type, asymtop_rot_channel_type
    use rotex__characters, only: add_trailing, i2c => int2char
    use rotex__arrays,     only: size_check
    use rotex__constants,  only: au2ev, au2cm

    implicit none (type, external)

    character(*), intent(in) :: prefix
      !! filename prefix
    character(*), intent(in) :: output_directory
      !! The directory in which output files are placed
    real(dp), intent(in) :: Eel_ex(:)
      !! The excitation electron energy grid in au
    real(dp), intent(in) :: Eel_dex(:)
      !! The de-excitation electron energy grid in au
    type(asymtop_rot_transition_type), intent(in) :: transition
      !! The transition to write to file
    real(dp),  intent(in) :: xs_xcite(:)
      !! Array of arrays of excitation cross sections for all spin multiplicities
    real(dp),  intent(in) :: xs_dxcite(:)
      !! Array of arrays of de-excitation cross sections for all spin multiplicities
    integer, intent(in) :: lmax
      !! The max value of l for the K/S-matrices

    integer :: ne, ie, iemin
    integer :: nlo, nup, kalo, kaup, kclo, kcup
    integer :: funit_ex, funit_dex
    real(dp) :: Elo,  Eup
    character(:), allocatable :: prefix_local, state_name1, state_name2,  filename_ex, filename_dex
    character(*), parameter :: XS_TYPE = "S-matrix + Coulomb-Born correction"
    type(asymtop_rot_channel_type) :: lo, up

    prefix_local = trim(prefix)
    call add_trailing(prefix_local, ".")

    ne = size(Eel_ex, 1)
    call size_check(Eel_dex,   ne, "EEL_DEX")
    call size_check(xs_xcite,  ne, "XS_XCITE")
    call size_check(xs_dxcite, ne, "XS_DXCITE")

    call mkdir(output_directory)

    lo = transition % lo
    up = transition % up

    ! -- filenames
    nlo  = lo % n
    nup  = up % n
    kalo = lo % ka
    kclo = lo % kc
    kaup = up % ka
    kcup = up % kc
    if(allocated(state_name1)) deallocate(state_name1)
    if(allocated(state_name2)) deallocate(state_name2)
    allocate(character(10) :: state_name1)
    allocate(character(10) :: state_name2)
    write(state_name1, '(I0, "_", I0, "_", I0)') Nlo,  Kalo,  Kclo
    write(state_name2, '(I0, "_", I0, "_", I0)') Nup,  Kaup,  Kcup
    state_name1 = trim(state_name1)
    state_name2 = trim(state_name2)

    Elo = lo % E
    Eup = up % E

    ! -- filenames
    filename_ex  = output_directory // prefix_local // state_name1 // "." // state_name2 // ".dat"
    filename_dex = output_directory // prefix_local // state_name2 // "." // state_name1 // ".dat"

    ! -- find first nonzero energy
    iemin = findloc(xs_xcite .gt. 0.0_dp, .true., 1)

    ! -- write data to file
    open(newunit = funit_ex,  file = filename_ex)
    open(newunit = funit_dex, file = filename_dex)
    call write_xs_header(funit_ex,  Nlo, Kalo, Kclo, Nup, Kaup, Kcup, XS_TYPE, i2c(lmax), "∞")
    call write_xs_header(funit_dex, Nup, Kaup, Kcup, Nlo, Kalo, Kclo, XS_TYPE, i2c(lmax), "∞")
    do ie = iemin, ne
      write(funit_ex,  ENERGY_XS_WRITE_FMT) Eel_ex(ie)  * au2ev, xs_xcite(ie)  * au2cm*au2cm
      write(funit_dex, ENERGY_XS_WRITE_FMT) Eel_dex(ie) * au2ev, xs_dxcite(ie) * au2cm*au2cm
    enddo
    close(funit_ex)
    close(funit_dex)

  end subroutine write_total_xs_to_file

  ! ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! pure elemental function xtrap_xs(xs1, E1, Etarg) result(res)
  !   !! Extrapolate a cross section as 1/E towards 0, given the cross section xs1 at
  !   !! energy E1
  !   use rotex__types, only: dp
  !   implicit none (type, external)
  !   real(dp), intent(in) :: xs1, E1, Etarg
  !   real(dp) :: res
  !   res = xs1*E1/Etarg
  ! end function xtrap_xs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine write_xs_header_linear(funit, N, Np, xs_type, lmax, lmax2)
    use rotex__characters, only: ndigits, i2c => int2char, state_label
    implicit none (type, external)
    integer, intent(in) :: funit
    integer, intent(in) :: N, Np
    character(:), allocatable :: nc, ncp
    character(*), intent(in) :: xs_type
    character(*), intent(in) :: lmax
    character(*), intent(in), optional :: lmax2
    ! -- determine how many characters to take up for N, Ka, and and values
    nc  = i2c(max(ndigits(N),  2) + 1)
    ncp = i2c(max(ndigits(Np), 2) + 1)
    write(funit, '("# Cross section type: ", A)') xs_type
    if(present(lmax2)) then
      write(funit, '("# S-matrix: l = 0 – ", A)') lmax
      write(funit, '("# Coulomb-Born correction: l = ", A, "..", A)') achar(ichar(lmax) + 1), lmax2
    else
      write(funit, '("# lmax: ", A)') lmax
    endif
    write(funit, '("# ", A' // nc  // ')',  advance = "no") "N"
    write(funit, '(2X,A)',                  advance = "no") "-->"
    write(funit, '(A'   // ncp // ')')                      "N'"
    write(funit, '(A)',    advance = "no") "# "
    write(funit, '(2X, A)',    advance = "no") state_label(N)
    write(funit, '(2X,A)', advance = "no") "-->"
    write(funit, '(1X, A)')                    state_label(Np)
    write(funit, '("# ", 2(A30))') "scattering energy (eV)", "cross section (cm²)"
  end subroutine write_xs_header_linear
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine write_xs_header_symtop(funit, N, K, Np, Kp, rchar, xs_type, lmax, lmax2)
    use rotex__characters,  only: ndigits, i2c => int2char, state_label
    use rotex__pointgroups, only: pg_nrot
    implicit none (type, external)
    integer, intent(in) :: funit
    integer, intent(in) :: N, K, Np, Kp, rchar
    character(:), allocatable :: nc, ncp
    character(*), intent(in) :: xs_type
    character(*), intent(in) :: lmax
    character(*), intent(in), optional :: lmax2
    ! -- determine how many characters to take up for N, Ka, and and values
    nc  = i2c(max(ndigits(N),  2) + 1)
    ncp = i2c(max(ndigits(Np), 2) + 1)
    write(funit, '("# Cross section type: ", A)') xs_type
    if(present(lmax2)) then
      write(funit, '("# S-matrix: l = 0 – ", A)') lmax
      write(funit, '("# Coulomb-Born correction: l = ", A, "..", A)') achar(ichar(lmax) + 1), lmax2
    else
      write(funit, '("# lmax: ", A)') lmax
    endif
    write(funit,   '("# The z-axis is aligned with the ", A, " axis")') G%SYMAXIS
    if(G%SYMTOP_REDUCE_PROJECTION .AND. modulo(K, pg_nrot(G%TARGET_POINT_GROUP)) .ne. 0) then
      write(funit, '("# Rotational character w.r.t. C2prime axis: n/a (±K folded)")')
    else
      write(funit, '("# Rotational character w.r.t. C2prime axis: ", I2)') rchar
    endif
    write(funit, '("# ", A' // nc  // ',",", A' // nc  // ')',  advance = "no") "N", "K"
    write(funit, '(2X,A)',                                         advance = "no") "-->"
    write(funit, '(A'   // ncp // ',",", A' // ncp // ')')                      "N'", "K'"
    write(funit, '(A)',    advance = "no") "# "
    write(funit, '(2X,A)',    advance = "no") state_label(N, K, sepstr=",  ")
    write(funit, '(2X,A)', advance = "no") "-->"
    write(funit, '(1X,A)')                    state_label(Np, Kp, sepstr=",  ")
    write(funit, '("# ", 2(A30))') "scattering energy (eV)", "cross section (cm²)"
  end subroutine write_xs_header_symtop
  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine write_xs_header_asymtop(funit, N, Ka, Kc, Np, Kap, Kcp, xs_type, lmax, lmax2)
    use rotex__characters, only: ndigits, i2c => int2char, state_label
    implicit none (type, external)
    integer, intent(in) :: funit
    integer, intent(in) :: N, Ka, Kc, Np, Kap, Kcp
    character(:), allocatable :: nc, ncp
    character(*), intent(in) :: xs_type
    character(*), intent(in) :: lmax
    character(*), intent(in), optional :: lmax2
    ! -- determine how many characters to take up for N, Ka, and and values
    nc  = i2c(max(ndigits(N),  2) + 1)
    ncp = i2c(max(ndigits(Np), 2) + 1)
    write(funit, '("# Cross section type: ", A)') xs_type
    if(present(lmax2)) then
      write(funit, '("# S-matrix: l = 0 – ", A)') lmax
      write(funit, '("# Coulomb-Born correction: l = ", A, "..", A)') achar(ichar(lmax) + 1), lmax2
    else
      write(funit, '("# lmax: ", A)') lmax
    endif
    write(funit,   '("# The z-axis is aligned with the ", A, " axis")') G%SYMAXIS
    write(funit, '("# ", 2(A' // nc  // ',","), A' // nc  // ')',  advance = "no") "N", "Ka", "Kc"
    write(funit, '(2X,A)',                                         advance = "no") "-->"
    write(funit, '(2(A'   // ncp // ',","), A' // ncp // ')')                      "N'", "Ka'", "Kc'"
    write(funit, '(A)',    advance = "no") "# "
    write(funit, '(2X,A)',    advance = "no") state_label(N, Ka, Kc, sepstr=",  ")
    write(funit, '(2X,A)', advance = "no") "-->"
    write(funit, '(1X,A)')                    state_label(Np, Kap, Kcp, sepstr=",  ")
    write(funit, '("# ", 2(A30))') "scattering energy (eV)", "cross section (cm²)"
  end subroutine write_xs_header_asymtop

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine write_channels_to_file( &
        filename                            &
      , jmin                                &
      , jmax                                &
      , n_states                            &
      , channels_l                          &
      , channels_l_j                        &
    )
    !! Write rotational channel info to file

    use rotex__types, only: asymtop_rot_channel_l_type, asymtop_rot_channel_l_vector_type &
                          , n_states_type, asymtop_rot_channel_type
    use rotex__channel_ops, only: operator(.eq.)
    use rotex__constants, only: au2ev
    use rotex__symmetry, only: spin_symmetry
    implicit none (type, external)
    character(*), intent(in) :: filename
      !! File to which we write channels
    integer, intent(in) :: jmin, jmax
      !! Total angular momentum mim/max
    type(n_states_type), intent(in) :: n_states(:)
    type(asymtop_rot_channel_l_type),        intent(in) :: channels_l(:)
      !! Rotational channels
    type(asymtop_rot_channel_l_vector_type), intent(in) :: channels_l_j(jmin:jmax)
      !! Rotational channels for each J
    type(asymtop_rot_channel_type) :: channel_without_l
    logical, allocatable :: jtest(:)
    integer  :: funit, nchans, j
    integer  :: in, n, itau, ka, kc, sym, ksym
    integer, allocatable :: lvals(:), jvals(:)
    real(dp) :: e
    nchans = size(channels_l, 1)
    open(newunit = funit, file = filename)
    select case(G%ROTOR_KIND)
    case("a","A")
      write(funit, '("# ", 3(A7), A15, A7, 3X, A4, A7)') "N", "Ka", "Kc", "E (meV)", "sym", "l", "J"
    case("s","S")
      write(funit, '("# ", 2(A7), A15, A7, 3X, A4, A7)') "N", "K", "E (meV)", "sym", "l", "J"
    case("l","L")
      write(funit, '("# ", (A7), A15, A7, 3X, A4, A7)') "N", "E (meV)", "sym", "l", "J"
    end select
    do in=1, size(n_states, 1)
      n  = n_states(in) % n
      do itau=1, 2*n+1

        ! -- for symmsteric tops, one of these is not allocated. for now, just print 0
        ka = 0 ; if(allocated(n_states(in) % ka)) ka = n_states(in) % ka(itau)
        kc = 0 ; if(allocated(n_states(in) % kc)) kc = n_states(in) % kc(itau)

        if(any(G%ROTOR_KIND .eq. ["s","S"])) then
          select case(G%SYMAXIS)
          case("a","A")
            Ksym = Ka
          case("c","C")
            Ksym = Kc
          end select
        endif

        e  = n_states(in) % eigenh % eigvals(itau)
        sym = spin_symmetry(n, ka, kc)
        channel_without_l = asymtop_rot_channel_type(nelec=1, e=e, n=n, ka=ka, kc=kc, sym=sym)

        ! -- which l values are inlcuded ?
        lvals = channels_l % l
        lvals = pack(lvals, channels_l(:) .eq. channel_without_l)

        ! -- which J values are included ?
        jtest = [(any(channels_l_j(j) % channels(:) .eq. channel_without_l), j=jmin, jmax)]
        jvals = pack([(j, j=jmin, jmax)], jtest)

        ! -- N, Ka, Kc, E
        select case(G%ROTOR_KIND)
        case("a","A")
          write(funit, '(2X, 3(I7), F15.8, I7)', advance="no") n, ka,   kc, e*au2ev*1000, sym
        case("s","S")
          write(funit, '(2X, 2(I7), F15.8, I7)', advance="no") n, ksym,     e*au2ev*1000, sym
        case("l","L")
          write(funit, '(2X, I7, F15.8, I7)',    advance="no") n,           e*au2ev*1000, sym
        end select

        ! -- l
        write(funit, '(I4)', advance="no") minval(lvals)
        if(size(lvals, 1) .gt. 2) then
          write(funit, '(A2,I0)', advance="no") "..", maxval(lvals)
        elseif(size(lvals, 1) .eq. 2) then
          write(funit, '(A1,I0)', advance="no") ",", maxval(lvals)
        endif

        ! -- J
        write(funit, '(I4)', advance="no") minval(jvals)
        if(size(jvals, 1) .gt. 2) then
          write(funit, '(A2,I0)', advance="no") "..", maxval(jvals)
        elseif(size(jvals, 1) .eq. 2) then
          write(funit, '(A1,I0)', advance="no") ",", maxval(jvals)
        endif
        write(funit, *)

      enddo
    enddo
    close(funit)
  end subroutine write_channels_to_file

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine write_elec_mat_elems_to_file( &
        kmat_eval_energies                        &
      , elec_channels                             &
      , eigenphases                               &
      , eigenphases_unwrapped                     &
      , sinmat                                    &
      , cosmat                                    &
      , spinmult_name                             &
    )
    !! Writes elements of K(E), sin(E), cos(E) to the corresponding files/directories determined by
    !! `matname`.

    use rotex__types,      only: elec_channel_type
    use rotex__constants,  only: au2ev
    use rotex__characters, only: i2c => int2char
    use rotex__arrays,     only: size_check

    implicit none (type, external)

    real(dp),                intent(in) :: Kmat_eval_energies(:)
      !! nE-element array of matrix evaluation energie
    type(elec_channel_type), intent(in) :: elec_channels(:)
      !! n-element array of electronic channels
    real(dp),                intent(in) :: eigenphases(:,:)
      !! Electronic eigenphases n×nE before identification and branch correction
    real(dp),                intent(in) :: eigenphases_unwrapped(:,:)
      !! Electronic eigenphases n×nE after identification and branch correction
    real(dp),                intent(in) :: sinmat(:,:,:), cosmat(:,:,:)
      !! Sine and cosine matrice n×n×nE
    character(*),            intent(in) :: spinmult_name

    integer :: n, ne, i, j, ie, nwrite
    integer :: funit_sine, funit_cosine, funit_eigenphases, funit_eigenphases2
    character(:), allocatable :: eigenphases_file, eigenphases_file2, sin_file, cos_file, phases_dir

    n  = size(eigenphases, 1)
    ne = size(eigenphases, 2)
    nwrite = (n*(n+1)) / 2 ! -- symemtric matrices, write only a triangle

    call size_check(elec_channels, n,   "electronic channels")
    call size_check(sinmat, [n, n, ne], "SIN MATRIX")
    call size_check(cosmat, [n, n, ne], "COS MATRIX")
    call size_check(eigenphases, shape(eigenphases_unwrapped), "EIGENPHASES_UNWRAPPED")

    phases_dir       = G%OUTPUT_DIRECTORY // "phases/"
    call mkdir(phases_dir    // spinmult_name)
    call mkdir(phases_dir    // spinmult_name)
    call mkdir(phases_dir    // spinmult_name)
    eigenphases_file = phases_dir    // spinmult_name // "/eigenphases.dat"
    eigenphases_file2 = phases_dir   // spinmult_name // "/eigenphases_smooth.dat"
    sin_file         = phases_dir    // spinmult_name // "/sine_elements.dat"
    cos_file         = phases_dir    // spinmult_name // "/cosine_elements.dat"

    open(newunit=funit_eigenphases,  file=eigenphases_file)
    open(newunit=funit_eigenphases2, file=eigenphases_file2)
    open(newunit=funit_sine,         file=sin_file)
    open(newunit=funit_cosine,       file=cos_file)

    ! -- channel header
    write(funit_eigenphases,  '("# ", '//i2c(n)//'(I0,",",I0,",",I0,2X))') &
      ( elec_channels(i) % nelec                                          &
      , elec_channels(i) % l                                              &
      , elec_channels(i) % ml, i=1,n)
    write(funit_eigenphases2, '("# ", '//i2c(n)//'(I0,",",I0,",",I0,2X))') &
      ( elec_channels(i) % nelec                                          &
      , elec_channels(i) % l                                              &
      , elec_channels(i) % ml, i=1,n)
    write(funit_sine, '("# ", '//i2c(nwrite)//'(I0,",",I0,",",I0," <-> "I0,",",I0,",",I0,2X))') &
      ((elec_channels(i) % nelec                                                               &
      , elec_channels(i) % l                                                                   &
      , elec_channels(i) % ml                                                                  &
      , elec_channels(j) % nelec                                                               &
      , elec_channels(j) % l                                                                   &
      , elec_channels(j) % ml                                                                  &
      , i=1, j), j=1,n)
    write(funit_cosine, '("# ", '//i2c(nwrite)//'(I0,",",I0,",",I0," <-> "I0,",",I0,",",I0,2X))') &
      ((elec_channels(i) % nelec                                                               &
      , elec_channels(i) % l                                                                   &
      , elec_channels(i) % ml                                                                  &
      , elec_channels(j) % nelec                                                               &
      , elec_channels(j) % l                                                                   &
      , elec_channels(j) % ml                                                                  &
      , i=1, j), j=1,n)

    ! -- data
    do ie=1, ne
      write(funit_eigenphases,  '(ES15.7,X)', advance='no') kmat_eval_energies(ie)*au2ev
      write(funit_eigenphases2, '(ES15.7,X)', advance='no') kmat_eval_energies(ie)*au2ev
      write(funit_sine,         '(ES15.7,X)', advance='no') kmat_eval_energies(ie)*au2ev
      write(funit_cosine,       '(ES15.7,X)', advance='no') kmat_eval_energies(ie)*au2ev
      do i=1,n
        write(funit_eigenphases,  '(ES15.7,X)', advance='no') eigenphases(i, ie)
        write(funit_eigenphases2, '(ES15.7,X)', advance='no') eigenphases_unwrapped(i, ie)
        do j=1,i
          write(funit_sine,         '(ES15.7,X)', advance='no') sinmat(i, j, ie)
          write(funit_cosine,       '(ES15.7,X)', advance='no') cosmat(i, j, ie)
        enddo
      enddo
      write(funit_eigenphases,*)
      write(funit_eigenphases2,*)
      write(funit_sine,*)
      write(funit_cosine,*)
    enddo

    close(funit_eigenphases)
    close(funit_eigenphases2)
    close(funit_sine)
    close(funit_cosine)

  end subroutine write_elec_mat_elems_to_file

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine write_Smat_J_elems_to_file(spinmult, J, kmat_eval_energies, channels, smat_flat)
    !! Writes elements of S^J(E) to disk for inspection

    use rotex__types,      only: asymtop_rot_channel_l_type
    use rotex__constants,  only: au2ev
    use rotex__globals,    only: spinmult_names
    use rotex__characters, only: i2c => int2char

    implicit none (type, external)

    integer, intent(in) :: spinmult
      !! The current spin multiplicity
    integer,                                 intent(in) :: J
      !! The current value of J
    real(dp),                                intent(in) :: Kmat_eval_energies(:)
      !! The evaluation energy grid
    type(asymtop_rot_channel_l_type), intent(in) :: channels(:)
      !! The rotational channels per J
    complex(dp),                      intent(in) :: smat_flat(:,:)
      !! The current flattened S^J sub block

    integer :: ie, i, i1, i2, k, ne, nchans, funitr, funitc, funitch
    character(:), allocatable :: fnamer, fnamec, fnamech

    ne = size(kmat_eval_energies, 1)

    call mkdir(G%OUTPUT_DIRECTORY // "/smat_J" // spinmult_names(spinmult))
    call mkdir(G%OUTPUT_DIRECTORY // "channels")

    nchans = size(channels, 1)

    fnamer  = G%OUTPUT_DIRECTORY // "/smat_J" // spinmult_names(spinmult) // "/S_J"   // trim(adjustl(i2c(J))) // "_real.dat"
    fnamec  = G%OUTPUT_DIRECTORY // "/smat_J" // spinmult_names(spinmult) // "/S_J"   // trim(adjustl(i2c(J))) // "_cplx.dat"
    fnamech = G%OUTPUT_DIRECTORY // "channels/channels_J" // trim(adjustl(i2c(J))) // ".txt"

    ! -- channels
    open(newunit=funitch,  file=fnamec,  status="replace", action="write")
    write(funitch, '("# J = ", I0)')
    write(funitch, '("# ")', advance="no")
    write(funitch, '(A6)', advance="no") "idx"
    write(funitch, '(2X, A2)', advance="no") "n"
    write(funitch, '(2X, A6)', advance="no") "N"
    write(funitch, '(2X, A6)', advance="no") "Ka"
    write(funitch, '(2X, A6)', advance="no") "Kc"
    write(funitch, '(2X, A2)', advance="no") "l"
    write(funitch, '(2X, A22)') "E (meV)"
    do i=1, nchans
      write(funitch, '(2X, I6, 2X, I2, 3(2X, I6), 2X, I2, 2X, ES22.14)') &
          i                                                            &
        , channels(i)%nelec                              &
        , channels(i)%N                                  &
        , channels(i)%Ka                                 &
        , channels(i)%Kc                                 &
        , channels(i)%l                                  &
        , channels(i)%E
    enddo
    close(funitch)

    ! -- S
    open(newunit=funitr,  file=fnamer,  status="replace", action="write")
    open(newunit=funitc,  file=fnamec,  status="replace", action="write")
    write(funitr, '("# J = ", I0)') J             ; write(funitc, '("# J = ", I0)') J
    write(funitr, '("# ")', advance="no")         ; write(funitc, '("# ")', advance="no")
    write(funitr, '(A15)', advance="no") "E (eV)" ; write(funitc, '(A15)') "E (eV)"
    write(funitr, '(" S(E)..")')                  ; write(funitr, '(" S(E)..")')
    do ie=1, ne
      write(funitr, '(ES17.7)', advance="no") Kmat_eval_energies(ie) * au2ev
      write(funitc, '(ES17.7)', advance="no") Kmat_eval_energies(ie) * au2ev
      k=0
      do i1=1, nchans ; do i2=1, i1
        k=k+1
        write(funitr, '(2X,ES22.12)', advance="no") smat_flat(k,ie)%re
        write(funitc, '(2X,ES22.12)', advance="no") smat_flat(k,ie)%im
      enddo ; enddo
      write(funitr, *)
      write(funitc, *)
    enddo

    close(funitr)
    close(funitc)

  end subroutine write_Smat_J_elems_to_file

! ================================================================================================================================ !
end module rotex__writing
