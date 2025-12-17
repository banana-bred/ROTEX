! ================================================================================================================================ !
module rotex__writing
  !! Procedures for writing data to disk

  implicit none

  private

  public :: write_lifetimes_to_file
  public :: write_channels_to_file
  public :: write_CB_xs_to_file
  public :: write_smat_xs_to_file
  public :: write_total_xs_to_file

  character(*), parameter :: ENERGY_XS_WRITE_FMT = '(2X, 2E30.20)'

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine write_lifetimes_to_file(output_directory, N_min, N_max, N_states, zaxis)
    !! Writes the states involved in the excitation and their lifetimes

    use rotex__types,     only: dp, N_states_type
    use rotex__constants, only: au2ev, au2sec

    implicit none

    character(*), intent(in) :: output_directory
      !! The directory in which output files are placed
    integer, intent(in) :: N_min, N_max
      !! Minimum and maximum value of N for which lifetimes were evaluated
    type(N_states_type), intent(in) :: N_states(:)
    character(1),        intent(in) :: zaxis
      !! The array of states

    integer :: funit
    integer :: inlo, i_tau
    integer :: N, Ka, Kc
    integer :: nc
    real(dp) :: E, EinstA, lifetime
    character(:), allocatable :: filename
    character(:), allocatable :: fmt
    character(15) :: lifetime_char

    filename = output_directory // "lifetimes.dat"
    open(newunit = funit, file = filename)

    nc = maxval(N_states(:) % N) + 1

    ! -- write file header
    write(funit, '("# The z-axis is aligned with the ", A, " axis")') zaxis
    write(funit, '("# ", 3A6, 2A16)') "N", "Ka", "Kc", "energy (meV)", "lifetime (s)"

    do inlo = 1, size(N_states, 1)

      N = N_states(inlo) % N
      if(N .lt. N_min) cycle
      if(N .gt. N_max) cycle

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

        write(funit, '(2X, 3I6)', advance = "no") N, Ka, Kc
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
    , zaxis                              &
    , E_el                                &
    , xs                                 &
    , init                               &
    , fin                                &
    , lmax                               &
    , xs_type                            &
    )
    !! Writes a Coulomb-Born cross section to a file whos name and file header
    !! carry information about the state symmetry
    use rotex__types,      only: dp, N_states_type, asymtop_rot_channel_type
    use rotex__characters, only: add_trailing,  sub, sup
    use rotex__constants,  only: au2eV, au2cm
    use rotex__functions,  only: logrange

    implicit none

    character(*),                      intent(in) :: prefix
      !! filename prefix
    character(*),                      intent(in) :: output_directory
      !! The directory in which output files are placed
    character(1),                      intent(in) :: zaxis
      !! The A, B, or C axis that lies along z
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
    call write_xs_header(funit,   zaxis, ni, kai, kci, nf, kaf, kcf, xs_type_, lmax)
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
    , zaxis                                &
    , egrid_total                          &
    , transition                           &
    , exxs                                 &
    , dexxs                                &
    , lmax                                 &
    )
    !! Writes an S-matrix (+CB) cross section to a file whos name and file header
    !! carry information about the state symmetry
    use rotex__types,      only: dp, rvector_type, asymtop_rot_transition_type,  asymtop_rot_channel_type
    use rotex__arrays,     only: size_check
    use rotex__system,     only: die, warn, mkdir
    use rotex__characters, only: add_trailing, sub, sup, i2c => int2char
    use rotex__constants,  only: au2eV, au2cm, pi
    use rotex__functions,  only: logrange

    implicit none

    character(*), intent(in) :: prefix
      !! filename prefix
    character(*), intent(in) :: output_directory
      !! The directory in which output files are placed
    character(1), intent(in) :: zaxis
      !! The A, B, or C axis that lies along z
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
    integer :: nlo, nup, kalo, kclo, kaup, kcup
    integer :: funit_ex, funit_dex
    real(dp) :: Elo, Eup, Eel_ex, Eel_dex, sigmaup, sigmadown
    character(:), allocatable :: state_name1
    character(:), allocatable :: state_name2
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

    ! -- get upper and lower channel energy
    Elo = lo % E
    Eup = up % E

    ! -- write (de-)excitation data
    filename_ex  = output_directory // prefix_local // state_name1 // "." // state_name2 // ".dat"
    filename_dex = output_directory // prefix_local // state_name2 // "." // state_name1 // ".dat"

    iemin = findloc(exxs .gt. 0.0_dp, .true., 1)

    open(newunit = funit_ex,  file = filename_ex)
    open(newunit = funit_dex, file = filename_dex)

    call write_xs_header(funit_ex,  zaxis, Nlo, Kalo, Kclo, Nup, Kaup, Kcup, "S-matrix", i2c(lmax))
    call write_xs_header(funit_dex, zaxis, Nup, Kaup, Kcup, Nlo, Kalo, Kclo, "S-matrix", i2c(lmax))

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
    , zaxis                                 &
    , Eel_ex                                &
    , Eel_dex                               &
    , transition                            &
    , xs_xcite                              &
    , xs_dxcite                             &
    , lmax                                  &
    )
    !! Write the total cross-sections (S-matrix + CB correction) to disk for a single transition

    use rotex__kinds,      only: dp
    use rotex__types,      only: asymtop_rot_transition_type, asymtop_rot_channel_type
    use rotex__characters, only: add_trailing, i2c => int2char
    use rotex__arrays,     only: size_check
    use rotex__system,     only: mkdir
    use rotex__constants,  only: au2ev, au2cm

    implicit none

    character(*), intent(in) :: prefix
      !! filename prefix
    character(*), intent(in) :: output_directory
      !! The directory in which output files are placed
    character(1), intent(in) :: zaxis
      !! The A, B, or C axis that lies along z
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
    call write_xs_header(funit_ex,  zaxis, Nlo, Kalo, Kclo, Nup, Kaup, Kcup, XS_TYPE, i2c(lmax), "∞")
    call write_xs_header(funit_dex, zaxis, Nup, Kaup, Kcup, Nlo, Kalo, Kclo, XS_TYPE, i2c(lmax), "∞")
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
  !   implicit none
  !   real(dp), intent(in) :: xs1, E1, Etarg
  !   real(dp) :: res
  !   res = xs1*E1/Etarg
  ! end function xtrap_xs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine write_xs_header(funit, axis, N, Ka, Kc, Np, kaup, kcup, xs_type, lmax, lmax2)
    use rotex__characters, only: ndigits, i2c => int2char
    implicit none
    integer, intent(in) :: funit
    character(*), intent(in) :: axis
    integer, intent(in) :: N, Ka, Kc, Np, kaup, kcup
    character(:), allocatable :: fmt, fmtp
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
    write(funit,   '("# The z-axis is aligned with the ", A, " axis")') axis
    write(funit, '("# ", 2(A' // nc  // ',","), A' // nc  // ')',  advance = "no") "N", "Ka", "Kc"
    write(funit, '(2X,A)',                                                advance = "no") "-->"
    write(funit, '(2(A'   // ncp // ',","), A' // ncp // ')')                    "N", "Ka", "Kc"
    fmt  = '(2(I' // nc  // ', ","), I' // nc  // ',)'
    fmtp = '(2(I' // ncp // ', ","), I' // ncp // ',)'
    write(funit, '(A)',     advance = "no") "# "
    write(funit, fmt,       advance = "no") N, Ka, Kc
    write(funit, '(2X,A)', advance = "no") "-->"
    write(funit, fmtp)                      Np, kaup, kcup
    write(funit, '("# ", 2(A30))') "scattering energy (eV)", "cross section (cm²)"
  end subroutine write_xs_header

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine write_channels_to_file(filename, jmin, jmax, n_states, channels_l, channels_l_j, spin_isomer_kind, symaxis)
    !! Write rotational channel info to file
    use rotex__types, only: dp, asymtop_rot_channel_l_type, asymtop_rot_channel_l_vector_type &
                          , n_states_type, asymtop_rot_channel_type, operator(.eq.)
    use rotex__constants, only: au2ev
    use rotex__symmetry, only: spin_symmetry
    implicit none
    character(*), intent(in) :: filename
      !! File to which we write channels
    integer, intent(in) :: jmin, jmax
      !! Total angular momentum mim/max
    integer, intent(in) :: spin_isomer_kind
      !! Kind of nuclear spin to preserve
    character(1), intent(in) :: symaxis
      !! Symmetry axis for nuclear spin
    type(n_states_type), intent(in) :: n_states(:)
    type(asymtop_rot_channel_l_type),        intent(in) :: channels_l(:)
      !! Rotational channels
    type(asymtop_rot_channel_l_vector_type), intent(in) :: channels_l_j(jmin:jmax)
      !! Rotational channels for each J
    type(asymtop_rot_channel_type) :: channel_without_l
    logical, allocatable :: jtest(:)
    integer  :: funit, nchans, j
    integer  :: in, n, itau, ka, kc, sym
    integer, allocatable :: lvals(:), jvals(:)
    real(dp) :: e
    nchans = size(channels_l, 1)
    open(newunit = funit, file = filename)
    write(funit, '("# ", 3(A7), A15, A7, 3X, A4, A7)') "N", "Ka", "Kc", "E (meV)", "sym", "l", "J"
    do in=1, size(n_states, 1)
      n  = n_states(in) % n
      do itau=1, 2*n+1
        ka = n_states(in) % ka(itau)
        kc = n_states(in) % kc(itau)
        e  = n_states(in) % eigenh % eigvals(itau)
        sym = spin_symmetry(n, ka, kc, spin_isomer_kind, symaxis)
        channel_without_l = asymtop_rot_channel_type(nelec=1, e=e, n=n, ka=ka, kc=kc, sym=sym)
        ! -- which l values are inlcuded ?
        lvals = channels_l % l
        lvals = pack(lvals, channels_l(:) .eq. channel_without_l)
        ! -- which J values are included ?
        jtest = [(any(channels_l_j(j) % channels(:) .eq. channel_without_l), j=jmin, jmax)]
        jvals = pack([(j, j=jmin, jmax)], jtest)
        ! -- N, Ka, Kc, E
        write(funit, '(2X, 3(I7), F15.8, I7)', advance="no") n, ka, kc, e*au2ev*1000, sym
        ! -- l
        write(funit, '(I4)', advance="no") minval(lvals)
        if(size(lvals, 1) .gt. 2) then
          write(funit, '(A2,I0)', advance="no") "..", maxval(lvals)
        elseif(size(lvals, 1) .eq. 2) then
          write(funit, '(A1,I0)', advance="no") ",", maxval(lvals)
        endif
        ! -- j
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

! ================================================================================================================================ !
end module rotex__writing
