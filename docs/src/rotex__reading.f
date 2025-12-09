! ================================================================================================================================ !
module rotex__reading
  !! Contains procedures used in reading data (K-matrices and namelist data)
  use rotex__constants, only: UKRMOLX, MQDTR2K

  implicit none

  private

  public :: read_kmats
  public :: read_namelists

! ================================================================================================================================ !
contains
! ================================================================================================================================ !


  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine read_kmats( kmat_dir          &
                              , channels_dir      &
                              , point_group       &
                              , spinmult          &
                              , kmat_lmax         &
                              , Kmat              &
                              , elec_channels     &
                              , channel_E_units   &
                              , kmat_eval_E_units &
                              , kmat_output_type  &
                              , kmat_e_closest    )
    !! Reads in a K-matrix from a file with a very particular file format given by kmat_output_type

    use rotex__types,      only: dp, elec_channel_type, permsort_channels, rvector_type, ivector_type
    use rotex__utils,      only: read_blank
    use rotex__arrays,     only: append, is_symmetric, realloc
    use rotex__system,     only: die, stdout, stderr
    use rotex__symmetry,   only: group_size, irrep_name
    use rotex__constants,  only: au2ev, IOSTAT_END, IOSTAT_OK, spinmult_names, DEFAULT_INT
    use rotex__characters, only: int2char

    implicit none

    character(*), intent(in) :: kmat_dir
      !! Directory in which the files containing the K-matrices
    character(*), intent(in) :: channels_dir
      !! Directory in which the channel data are located
    character(*), intent(in) :: point_group
      !! The point group of the calculations
    integer, intent(in) :: spinmult
      !! The spin multiplicity (2S+1) of the system (target + e⁻)
    integer, intent(in) :: kmat_lmax
      !! The max value of l in the electronic partial wave basis
    real(dp), intent(out), allocatable :: Kmat(:,:)
      !! K(i, j)
    type(elec_channel_type), intent(out), allocatable :: elec_channels(:)
      !! The channel basis of the K-matrix: \(n,l,λ\) (the code calls λ \(m_l\))
    character(1), intent(in) :: channel_E_units
      !! The units of the channel energies in the Kmat file. Options are :
      !!  - "r" for Rydberg, "h" for hartree, "e" for eV
    character(1), intent(in) :: kmat_eval_E_units
      !! The units of the energy at which the K-matrix was evaluated in the Kmat file. Options are :
      !!  - "h" for hartree
      !!  - "e" for eV
      !!  - "r" for Rydberg
    real(dp), intent(in) :: kmat_e_closest
      !! Evaluate the K-matrix that is closest to this energy
    character(*), intent(in) :: kmat_output_type
      !! The kind of K-matrix output to read

    logical :: skip_this_irrep
    integer :: i, j, ichan, iflat, i1, i2
    integer :: irrep
    integer :: nirreps
    integer :: nchans_total
    integer, allocatable :: nchans_irrep(:), idx(:)
    character(:), allocatable :: kmat_filename, channels_filename, filename, irrepname

    real(dp), allocatable :: channel_e_convert, kmat_e_convert

    type(elec_channel_type), allocatable :: elec_channels_this_irrep(:)
    type(ivector_type), allocatable :: index_map(:)
    type(rvector_type), allocatable :: kmat_flat(:)

    nirreps = group_size(point_group)
    allocate(nchans_irrep(nirreps))

    write(stdout, '(A)')     "----------------------------------------"
    write(stdout, '(A, I0)') "Spin multiplicity: ", spinmult
    write(stdout, '(A)')     "⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻"

    ! -- set the energy conversion scheme
    select case(channel_E_units)
    case("h") ; channel_e_convert = 1._dp
    case("e") ; channel_e_convert = 1._dp / au2ev
    case("r") ; channel_e_convert = 1._dp / 2._dp
    case default
      call die("Unable to determine the K-matrix channel energy units. Current value : " // channel_E_units // ". Supported&
        & values are 'r'ydberg, 'e'lectron-volt, and 'h'artree.")
    end select
    select case(Kmat_eval_E_units)
    case("h") ; kmat_e_convert = 1._dp
    case("e") ; kmat_e_convert = 1._dp / au2ev
    case("r") ; kmat_e_convert = 1._dp / 2._dp
    case default
      call die("Unable to determine the K-matrix evaluation energy units. Current value : " // kmat_eval_E_units // ". Supported&
        & values are 'r'ydberg, 'e'lectron-volt, and 'h'artree.")
    end select

    allocate(index_map(nirreps))
    allocate(kmat_flat(nirreps))

    nchans_total = 0

    write(stdout, '("Point group: ", A)') point_group
    ! -- read the K-matrices and electronic channels
    irrep_loop_kmats: do irrep = 1, nirreps
      irrepname = irrep_name(irrep, point_group)
      write(stdout, '(2X, "Irrep: ", A)') irrepname
      filename = kmat_dir // int2char(spinmult) // irrepname // ".kmat"

      select case(kmat_output_type)
      case(UKRMOLX)

        channels_filename = channels_dir // "channels.geom1." // spinmult_names(spinmult) // "." // irrepname
        kmat_filename     = kmat_dir     // "K-matrix.geom1." // spinmult_names(spinmult) // "." // irrepname

        call get_flat_kmat_and_channels_ukrmolx( &
            channels_filename                    &
          , kmat_filename                        &
          , kmat_flat(irrep)%vec                 &
          , channel_E_convert                    &
          , kmat_e_convert                       &
          , kmat_e_closest                       &
          , elec_channels_this_irrep             &
          , nchans_irrep(irrep)                  &
          , skip_this_irrep)

      case(MQDTR2K)

        ! -- the kmat file is expected to have the channels
        kmat_filename = kmat_dir // int2char(spinmult) // irrep_name(irrep, point_group) // ".kmat"

        call get_flat_kmat_and_channels_mqdtr2k(   &
            kmat_filename                          &
          , kmat_flat(irrep)%vec                   &
          , channel_E_convert                      &
          , kmat_e_convert                         &
          , kmat_e_closest                         &
          , elec_channels_this_irrep               &
          , nchans_irrep(irrep)                    &
          , skip_this_irrep&
        )

      case default
        call die("KMAT_OUTPUT_TYPE must be "//UKRMOLX//" or "//MQDTR2K)
      end select

      ! -- if we need to focus on channel parity
      if(point_group .eq. "cs") call fill_parity_array_this_irrep(kmat_lmax, elec_channels_this_irrep, irrep, point_group)

      call append(elec_channels, elec_channels_this_irrep)

      if(skip_this_irrep .eqv. .true.) cycle irrep_loop_kmats

      ! -- total number of channels across all irreps
      nchans_total = nchans_total + nchans_irrep(irrep)

      ! -- index map: this irrep → full basis of channels
      allocate(index_map(irrep)%vec(nchans_irrep(irrep)))
      do concurrent( ichan=1:nchans_irrep(irrep) )
        index_map(irrep)%vec(ichan) = find_global_channel_index(elec_channels_this_irrep(ichan), elec_channels)
      enddo

    enddo irrep_loop_kmats

    deallocate(elec_channels_this_irrep)

    allocate(Kmat(nchans_total, nchans_total))
    Kmat = 0

    ! -- flattened K-matrices -> K-matrix for all irreps
    do irrep = 1, nirreps
      ! -- flattened K-matrix for this irrep -> full symmetric K-matrix for all irreps
      iflat = 0
      do i = 1, nchans_irrep(irrep)
        i1 = index_map(irrep)%vec(i)
        do j = i, nchans_irrep(irrep)
          iflat = iflat + 1
          i2 = index_map(irrep)%vec(j)
          Kmat(i1, i2) = kmat_flat(irrep)%vec(iflat)
          Kmat(i2, i1) = kmat_flat(irrep)%vec(iflat)
        enddo
      enddo
    enddo

    deallocate(index_map, kmat_flat)

    ! -- at this point the K-matrix is block-diagonal w.r.t. irrep, as it should be.
    !    Sort the channels by quantum number, and permute the elements of K accordingly
    call permsort_channels(elec_channels, idx)
    Kmat = Kmat(idx, idx)
    call print_channels(elec_channels, stdout)

    if(is_symmetric(Kmat) .eqv. .true.) return

    write(stderr, '("maxval(abs(Kmat - transpose(Kmat))): ", E20.10)') maxval(abs(Kmat - transpose(Kmat)))
    call die("The K-matrix is not symmeric !")

  end subroutine read_kmats

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine fill_parity_array_this_irrep(lmax_kmat, elec_channels, irrep, point_group)
    !! Fill the array M_PARITY, stored in the module ROTEX__SYMMETRY, that will later be accessed
    !! by the rotational frame transformation to determine which values of m from -lmax_kmat to lmax_kmat
    !! correspond to even and odd combinations of partial waves
    use rotex__types,      only: elec_channel_type
    use rotex__system,     only: die
    use rotex__symmetry,   only: Ap, App, m_parity, even, odd
    use rotex__characters, only: i2c => int2char
    implicit none
    integer, intent(in) :: lmax_kmat
      !! Max value of l for the electronic channels
    type(elec_channel_type), intent(in) :: elec_channels(:)
      !! The electronic channels for the current irrep
    integer, intent(in) :: irrep
      !! The current irrep index
    character(*), intent(in) :: point_group
      !! The point group for the scattering calculations
    integer :: ichan, m
    if(point_group .ne. "cs") call die("Attempting to fill m_parity array for a point group&
      & other than Cs: " // point_group)
    if(allocated(m_parity) .eqv. .false.) allocate(m_parity(-lmax_kmat:lmax_kmat), source = 0)
    do ichan=1, size(elec_channels, 1)
      m = elec_channels(ichan) % ml
      if(m_parity(m) .ne. 0) cycle ! skip if set, but this probably should not happen
      select case(irrep)
      case(Ap)
        m_parity(m) = even
      case(App)
        m_parity(m) = odd
      case default
        call die("Somehow, irrep ("//i2c(irrep)//") is not one of the valid values for the point group " // point_group)
      end select
    enddo
  end subroutine fill_parity_array_this_irrep

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine get_flat_kmat_and_channels_ukrmolx( &
        channels_filename                        &
      , kmat_filename                            &
      , kmat_flat                                &
      , channel_e_convert                        &
      , kmat_e_convert                           &
      , kmat_e_closest                           &
      , elec_channels_this_irrep                 &
      , nchans_this_irrep                        &
      , skip_this_irrep)
    !! Return the flattened (1D) K-matrix that is closest to the desired evaluation energy
    !! given by kmat_e_closest

    use rotex__constants, only: IOSTAT_END, IOSTAT_OK, au2ev
    use rotex__kinds,     only: dp
    use rotex__arrays,    only: realloc, size_check, append
    use rotex__types,     only: elec_channel_type
    use rotex__utils,     only: read_blank
    use rotex__system,    only: stdout, stderr, die

    implicit none

    character(*), intent(in) :: channels_filename
    character(*), intent(in) :: kmat_filename
    real(dp), intent(inout), allocatable :: kmat_flat(:)
    real(dp), intent(in) :: channel_e_convert
    real(dp), intent(in) :: kmat_e_convert
    real(dp), intent(in) :: kmat_e_closest
    type(elec_channel_type), intent(out), allocatable :: elec_channels_this_irrep(:)
    integer, intent(out) :: nchans_this_irrep
    logical, intent(out) :: skip_this_irrep

    integer,      parameter :: UKRMOL_KMAT_ELEMENTS_PER_LINE = 4
    character(8), parameter :: UKRMOL_KMAT_ELEMENTS_FMT = '(D20.13)'

    integer  :: iostat, ne, i, ichan, ie, ie_closest, l, ml, nelec, nchans, iflat
    integer  :: funit, nchans_max, nskip, nskip_header, iline, icol, nlines, nchans_flat
    real(dp) :: E
    real(dp), allocatable :: kmat_energies(:)
    type(elec_channel_type) :: chan

    skip_this_irrep = .false.

    ! -- check file existence
    inquire(file = kmat_filename, iostat = iostat)
    if(iostat .ne. IOSTAT_OK) then
      skip_this_irrep = .true.
      return
    endif
    inquire(file = channels_filename, iostat = iostat)
    if(iostat .ne. IOSTAT_OK) call die("K-matrix is here for this irrep, but the channel file is not !")

    ! -- note that the K-matrices are expressed in the basis of OPEN channels, so the
    !    included channels changes when thresholds are crossed. This changes the number
    !    of channels and obviously the size of the resulting K-matrix.
    ! -- Count the number of energies/K-matrices
    open(newunit = funit, file = kmat_filename)
    ne = 0
    call read_blank(funit, 3)
    read(funit, *) i, i, i, nchans_max, nskip_header
    call read_blank(funit, nskip_header)
    nskip_header = nskip_header + 4 ! for the next re-reads
    do
      read(funit, *, iostat = iostat) nchans, i, nchans_flat, E
      if(iostat .eq. IOSTAT_END) exit
      ne = ne + 1
      nskip = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind = dp))
      call read_blank(funit, nskip)
    enddo

    call realloc(kmat_flat, nchans_flat)
    kmat_flat = 0._dp

    ! -- read in the K-matrix energies
    rewind(funit)
    call realloc(kmat_energies, ne)
    call read_blank(funit, nskip_header)
    ie = 0
    do
      read(funit, *, iostat = iostat) nchans, i, nchans_flat, E
      if(iostat .eq. IOSTAT_END) exit
      ie = ie + 1
      kmat_energies(ie) = E * kmat_e_convert
      nskip = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind=dp))
      call read_blank(funit, nskip)
    enddo

    ! -- find the lowest energy
    ie_closest = minloc(abs(kmat_energies - kmat_e_closest), 1)
    if(ie_closest .lt. 1) call die("Somehow, IE_CLOSEST returned a non-positive integer !")
    write(stdout, '(4X, "User requested K-matrix at ", E20.10, " eV")') kmat_e_closest            * au2ev
    write(stdout, '(7X, "Found Kmatrix at energy ",    E20.10, " eV")') kmat_energies(ie_closest) * au2ev

    ! -- Now, actually go and read that K-matrix
    rewind(funit)
    call read_blank(funit, nskip_header)
    ie = 0
    kmat_read_loop: do
      read(funit, *, iostat = iostat) nchans, i, nchans_flat, E
      if(iostat .eq. IOSTAT_END) then
        call die("Reach end of ")
        write(stderr, '("Number of K-matrices/energies: ", I0)') ne
        write(stderr, '("Target K-matrix energy: ", E20.10)') kmat_e_closest * au2ev
        write(stderr, '("Closest available K-matrix is number ", I0)') ie_closest
        call die("Could not find the K-matrix that is closest to the given target energy before EOF")
      endif
      ie = ie + 1
      if(ie .ne. ie_closest) then
        nskip = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind=dp))
        call read_blank(funit, nskip)
      else
        nlines = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind=dp))
        iflat = 0
        ! -- iterate through the lines and columns of the flattened K-matrix
        !    in the file
        do iline = 1, nlines
          do icol = 1, UKRMOL_KMAT_ELEMENTS_PER_LINE
            iflat = iflat + 1
            read(funit, UKRMOL_KMAT_ELEMENTS_FMT, advance = 'no') kmat_flat(iflat)
            if(iflat .eq. nchans_flat) exit kmat_read_loop
          enddo
          read(funit, *)
        enddo
        write(stderr, '("NLINES: ", I0)') nlines
        write(stderr, '("UKRMOL_KMAT_ELEMENTS_PER_LINE: ", I0)') UKRMOL_KMAT_ELEMENTS_PER_LINE
        write(stderr, '("NCHANS: ", I0)') nchans
        write(stderr, '("NCHANS_FLAT: ", I0)') nchans_flat
        call die("Improper K-matrix read loop exit. ")
      endif
    enddo kmat_read_loop

    close(funit)

    nchans_this_irrep = nchans
    if(nchans .lt. 1) then
      write(stderr, '("NCHANS: ", I0)') nchans
      call die("The K-matrix cannot have less than one channel !")
    endif
    call realloc(elec_channels_this_irrep, nchans_this_irrep)

    ! -- time to read the channels
    open(newunit = funit, file = channels_filename)
    write(stdout, '(A)') "Reading channels file at " // channels_filename

    call read_blank(funit, 2)
    ! -- read number of electronic states and number of electronic channels
    read(funit, *) nelec, i, i, nchans_max
    if(nchans_max .lt. 1) then
      write(stderr, '("NCHANS_MAX: ")') nchans_max
      call die("Total number of channels computed for this irrep is < 1 !")
    endif

    ! -- skip electronic state lines
    call read_blank(funit, nelec + 1)

    ! -- channel read
    !    There is no iq from UKRmol+. Coulomb f and g have the default normalization q=1
    do ichan = 1, nchans_this_irrep

      read(funit, *) i, nelec, l, ml, E

      ! -- convert channel energy to atomic units
      E = E * channel_E_convert

      chan = elec_channel_type(nelec=nelec, l=l, ml=ml, E=E)
      elec_channels_this_irrep(ichan) = chan

    enddo

    close(funit)

  end subroutine get_flat_kmat_and_channels_ukrmolx

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine get_flat_kmat_and_channels_mqdtr2k( &
        kmat_filename                            &
      , kmat_flat                                &
      , channel_e_convert                        &
      , kmat_e_convert                           &
      , kmat_e_closest                           &
      , elec_channels_this_irrep                 &
      , nchans_this_irrep                        &
      , skip_this_irrep)
    !! Return the flattened (1D) K-matrix that is closest to the desired evaluation energy
    !! given by kmat_e_closest, where the K-matrix channels format is that of MQDTR2K

    use rotex__constants, only: IOSTAT_END, IOSTAT_OK, au2ev
    use rotex__kinds,     only: dp
    use rotex__arrays,    only: realloc, size_check
    use rotex__types,     only: elec_channel_type
    use rotex__utils,     only: read_blank
    use rotex__system,    only: stdout, die

    implicit none

    character(*), intent(in) :: kmat_filename
    real(dp), intent(inout), allocatable :: kmat_flat(:)
    real(dp), intent(in) :: channel_e_convert
    real(dp), intent(in) :: kmat_e_convert
    real(dp), intent(in) :: kmat_e_closest
    type(elec_channel_type), intent(out), allocatable :: elec_channels_this_irrep(:)
    integer, intent(out) :: nchans_this_irrep
    logical, intent(out) :: skip_this_irrep

    integer  :: iostat, ne, funit, ne_skip, i, ichan, ie, ie_closest, l, ml, iq, nelec
    integer  :: nchans_irrep_flat
    real(dp) :: E
    real(dp), allocatable :: kmat_energies(:)

    skip_this_irrep = .false.

    inquire(file = kmat_filename, iostat = iostat)
    if(iostat .ne. IOSTAT_OK) then
      skip_this_irrep = .true.
      return
    endif

    open(newunit = funit, file = kmat_filename)

    ! -- get number of channels, skip to K-matrices, determine number of energies/K-matrices to read
    ne = 0
    call read_blank(funit)
    read(funit, *) nchans_this_irrep
    nchans_irrep_flat = nchans_this_irrep * (nchans_this_irrep+1) / 2
    call read_blank(funit, nchans_this_irrep)
    do
      ! -- energies are the first number, ignore the rest
      read(funit, *, iostat = iostat) E
      if(iostat .eq. IOSTAT_END) exit
      ne = ne + 1
    enddo
    rewind(funit)
    ! -- skip header and channels; read K-matrix evaluation energies
    call realloc(kmat_energies, ne)
    call read_blank(funit, 2+nchans_this_irrep)
    ie = 0
    do
      ! -- energies are the first number, ignore the rest
      read(funit, *, iostat = iostat) E
      if(iostat .eq. IOSTAT_END) exit
      ie = ie + 1
      kmat_energies(ie) = E * kmat_e_convert
    enddo

    ! -- find the lowest energy
    ie_closest = minloc(abs(kmat_energies - kmat_e_closest), 1)
    if(ie_closest .lt. 1) call die("Somehow, IE_CLOSEST returned a non-positive integer !")
    ne_skip = ie_closest - 1
    write(stdout, '(4X, "User requested K-matrix at ", E20.10, " eV")') kmat_e_closest            * au2ev
    write(stdout, '(7X, "Found Kmatrix at energy ",    E20.10, " eV")') kmat_energies(ie_closest) * au2ev

    ! -- skip header, read channels for this irrep
    rewind(funit)
    call read_blank(funit, 2)
    call realloc(elec_channels_this_irrep, nchans_this_irrep)
    do ichan=1, nchans_this_irrep
      read(funit, *) i, nelec, l, ml, E, iq
      ! -- convert channel energy to atomic units
      E = E * channel_E_convert
      elec_channels_this_irrep(ichan) = elec_channel_type(nelec = nelec, l = l, ml = ml, iq = iq, E = E)
    enddo

    ! -- skip to the desired K-matrix
    call read_blank(funit, ne_skip)

    ! -- read the evaluation energy and the flattened K-matrix. We're only interested in the
    !    K-matrix now that we've identified the evaluation energy
    call realloc(kmat_flat, nchans_irrep_flat)
    kmat_flat = 0.0_dp
    read(funit, *) E, (kmat_flat(i), i=1, nchans_irrep_flat)

    close(funit)

  end subroutine get_flat_kmat_and_channels_mqdtr2k

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_channels(channels, funit)
    use rotex__kinds,  only: dp
    use rotex__types,  only: elec_channel_type
    use rotex__system, only: stdout
    implicit none
    type(elec_channel_type), intent(in) :: channels(:)
    integer, intent(in), optional :: funit
    integer :: nelec, l, ml, iq, funit_, ichan, nchans
    real(dp) :: E
    funit_ = stdout ; if(present(funit)) funit_ = funit
    nchans = size(channels, 1)
    write(funit_, *)
    write(funit_, *)
    write(funit_, '(A)') "Electronic channels:"
    write(funit_, '(A)') "⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻⁻"
    write(funit_, '(8X, 4A4, A15)') "n", "l", "ml", "iq", "E (eV)"
    do ichan=1, nchans
      nelec = channels(ichan) % nelec
      l     = channels(ichan) % l
      ml    = channels(ichan) % ml
      iq    = channels(ichan) % iq
      E     = channels(ichan) % E
      write(funit_, '(4X, 5I4, E15.7)') ichan, nelec, l, ml, iq, E
    enddo
    write(funit_, *)
  end subroutine print_channels

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function find_global_channel_index(channel, channels) result(i)
    use rotex__types,  only: operator(.ne.), elec_channel_type
    use rotex__system, only: die
    implicit none
    type(elec_channel_type), intent(in) :: channel
    type(elec_channel_type), intent(in) :: channels(:)
    integer :: i, n
    n = size(channels, 1)
    do i = 1, n
      if(channel .ne. channels(i)) cycle
      return
    enddo
    call die("Failed to find the given channel !")
  end function find_global_channel_index

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine read_namelists(cfg)
    !! Reads user parameters and puts them into the config derived type
    use rotex__types,      only: dp, config_type, cd4_type, cd6_type
    use rotex__arrays,     only: append, remove_value
    use rotex__system,     only: stdin, stdout, ds => directory_separator, die
    use rotex__constants,  only: au2invcm, au2ev, macheps => macheps_dp, au2cm, au2deb, DEFAULT_CHAR1&
                               , UKRMOLX, MQDTR2K
    use rotex__characters, only: add_trailing, to_lower

    implicit none

    type(config_type), intent(out) :: cfg

    integer, parameter :: DEFAULT_INT = huge(1)

    ! -- namelist: control
    integer :: Nmin
    integer :: Nmax
    logical :: use_kmat
    logical :: use_CB
    integer :: spin_isomer_kind = 0
    character(:), allocatable :: output_directory
    character(1) :: rotor_kind = DEFAULT_CHAR1
    character(1) :: zaxis
    real(dp) :: abc(3) = 0.0_dp
    integer :: target_charge = DEFAULT_INT
    logical :: add_cd4 = .false.
    logical :: add_cd6 = .false.
    ! -- cd4
    real(dp) :: dn     = 0.0_dp
    real(dp) :: dnk    = 0.0_dp
    real(dp) :: dk     = 0.0_dp
    real(dp) :: deltan = 0.0_dp
    real(dp) :: deltak = 0.0_dp
    ! -- cd6
    real(dp) :: hn     = 0.0_dp
    real(dp) :: hnk    = 0.0_dp
    real(dp) :: hkn    = 0.0_dp
    real(dp) :: hk     = 0.0_dp
    real(dp) :: etan   = 0.0_dp
    real(dp) :: etank  = 0.0_dp
    real(dp) :: etak   = 0.0_dp

    ! -- namelist: kmat
    logical :: real_spherical_harmonics = .true.
    integer :: lmax_kmat = DEFAULT_INT
    integer :: num_egrid_segs
    integer, allocatable :: num_egrid(:)
    integer, allocatable :: spinmults(:)
    real(dp), allocatable :: egrid_segs(:)
    character(1) :: channel_energy_units_override = DEFAULT_CHAR1
    character(1) :: kmat_energy_units_override    = DEFAULT_CHAR1
    character(3) :: egrid_spacing
    character(7) :: kmat_output_type = "======="
    character(:), allocatable :: point_group
    character(:), allocatable :: kmat_dir
    character(:), allocatable :: channels_dir

    ! -- namelist: coulomb
    logical :: use_CDMS_einstA = .false.
    logical :: only_einsta = .false.
    logical :: do_xtrap = .false.
    logical :: do_dipole = .true.
    logical :: do_quadrupole = .false.
    logical :: analytic_total_cb(2)
    integer :: nE = DEFAULT_INT
    integer :: nE_xtrap = DEFAULT_INT
    integer :: lmax_partial = DEFAULT_INT
    integer :: lmax_total = DEFAULT_INT
    real(dp) :: Ei = 0.0_dp
    real(dp) :: Ef = 0.0_dp
    real(dp) :: Ei_xtrap = 0.0_dp
    real(dp) :: cartesian_dipole_moments(3)
    ! real(dp) :: cartesian_quadrupole_moments(6)
    real(dp) :: xs_zero_threshold = 0.0_dp   ! include all cross sections by default
    real(dp) :: kmat_energy_closest = 0.0_dp ! just take the first one
    character(:), allocatable :: CDMS_file

    namelist / control_namelist /                 &
      !! Contains parameters and values that are necessary to run the program
                         output_directory         &
                       , spin_isomer_kind         &
                       , nmin                     &
                       , nmax                     &
                       , use_kmat                 &
                       , use_cb                   &
                       , zaxis                    &
                       , rotor_kind               &
                       , target_charge            &
                       , abc                      &
                       , add_cd4                  &
                       , add_cd6                  &
                       , dn, dnk, dk, deltan, deltak &
                       , hn, hnk, hkn, hk, etan, etank, etak &
                       , xs_zero_threshold

    namelist / kmat_namelist /                      &
      !! Parameters regarding the K-matrces used for (de-excitation)
                      kmat_dir                      &
                    , channels_dir                  &
                    , lmax_kmat                     &
                    , point_group                   &
                    , num_egrid_segs                &
                    , num_egrid                     &
                    , egrid_segs                    &
                    , egrid_spacing                 &
                    , spinmults                     &
                    , kmat_output_type              &
                    , kmat_energy_closest           &
                    , real_spherical_harmonics      &
                    , channel_energy_units_override &
                    , kmat_energy_units_override

    namelist / coulomb_namelist /                     &
      !! Parameters regarding the Coulomb-Born approximation
      !! used for (de-)excitation
                         use_CDMS_einstA              &
                       , only_einsta                  &
                       , cdms_file                    &
                       , ei                           &
                       , ef                           &
                       , ne                           &
                       , ne_xtrap                     &
                       , do_xtrap                     &
                       , ei_xtrap                     &
                       , cartesian_dipole_moments     &
                       ! , cartesian_quadrupole_moments &
                       , do_dipole                    &
                       , do_quadrupole                &
                       , analytic_total_cb            &
                       , lmax_partial                 &
                       , lmax_total

    !!!!!!!!!!!!!!!!!!!!!! CONTROL_NAMELIST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(character(1000) :: output_directory)
    ! -- read
    read(stdin, control_namelist)
    ! -- check defaults
    if(Nmin          .eq. DEFAULT_INT)   call die("Must specify NMIN in CONTROL_NAMELIST")
    if(Nmax          .eq. DEFAULT_INT)   call die("Must specify NMAX in CONTROL_NAMELIST")
    if(rotor_kind    .eq. DEFAULT_CHAR1) call die("Must specify ROTOR_KIND in CONTROL_NAMELIST")
    if(target_charge .eq. DEFAULT_INT)   call die("Must specify TARGET_CHARGE in CONTROL_NAMELIST")
    if(ZAXIS         .eq. DEFAULT_CHAR1) call die("Must specify ZAXIS in CONTROL_NAMELIST")
    if(any(ABC       .eq. 0.0_dp))       call die("Must specify nonzero rotational constants ABC in CONTROL_NAMELIST")
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! KMAT_NAMELIST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(use_kmat .eqv. .true.) then
      ! -- prepare for reading
      rewind(stdin)
      allocate(character(1000) :: kmat_dir)
      allocate(character(1000) :: channels_dir)
      allocate(num_egrid(100))
      allocate(egrid_segs(101))
      allocate(character(10)   :: point_group)
      allocate(spinmults(10))
      spinmults = DEFAULT_INT
      kmat_dir(1:1)     = DEFAULT_CHAR1
      channels_dir(1:1) = DEFAULT_CHAR1
      kmat_output_type(1:1) = DEFAULT_CHAR1
      ! -- read
      read(stdin, kmat_namelist)
      ! -- trim arrays
      num_egrid  = num_egrid(1:num_egrid_segs)
      egrid_segs = egrid_segs(1:num_egrid_segs+1) / au2ev
      ! -- normalize characters
      call to_lower(egrid_spacing)
      call to_lower(kmat_output_type)
      select case(kmat_output_type)
      case(UKRMOLX, MQDTR2K)
        continue
      case default
        call die("KMAT_OUTPUT_TYPE in KMAT_NAMELIST must be one of " // UKRMOLX // " or " // MQDTR2K)
      end select
      select case(egrid_spacing)
        case("lin", "log") ; continue
        case default ; call die("EGRID_SPACING (" // egrid_spacing // ") must be LIN or LOG in KMAT_NAMELIST")
      end select
      ! -- check defaults
      if(kmat_dir(1:1) .eq. DEFAULT_CHAR1) call die("Must specify KMAT_DIR in KMAT_NAMELIST")
      if(lmax_kmat .eq. DEFAULT_INT .OR. lmax_kmat .lt. 0) &
        call die("LMAX_KMAT in KMAT_NAMELIST must be defined and be non-negative")
      call remove_value(spinmults, DEFAULT_INT)
      if(any(spinmults .lt. 1)) call die("SPINMULTS in KMAT_NAMELIST cannot have values that are < 1")
      if(     channels_dir(1:1) .eq. DEFAULT_CHAR1 &
        .AND. kmat_output_type  .eq. UKRMOLX) call die("Must specify CHANNELS_DIR in KMAT_NAMELIST with&
          & KMAT_OUPUT_TYPE = " // UKRMOLX)
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!! COULOMB_NAMELIST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(use_CB .eqv. .true.) then
      allocate(character(1000) :: cdms_file)
      ! -- prepare for reading
      cdms_file(1:1) = DEFAULT_CHAR1
      rewind(stdin)
      ! -- read
      read(stdin, coulomb_namelist)
      ! -- check values
      if(use_CDMS_einstA) then
        if(CDMS_file(1:1) .eq. DEFAULT_CHAR1) call die("Must define CDMS_FILE in COULOMB_NAMELIST&
          & when USE_CDMS_EINSTA is .TRUE.")
      endif
      if(do_dipole .eqv. .false.) call die("DO_DIPOLE in COULOMB_NAMELIST should not be set to .FALSE.; nothing would be done")
      if(do_quadrupole) call die("DO_QUADRUPOLE in COULOMB_NAMELIST should not be set to true; it is not implemented")
      if(Ei .eq. 0.0_dp) call die("EI in COULOMB_NAMELIST must be defined and be positive")
      if(Ef .eq. 0.0_dp) call die("EF in COULOMB_NAMELIST must be defined and be positive")
      if(Ef .le. Ei) call die("EF ≤ EI in COULOMB_NAMELIST is not allowed")
      if(nE .eq. DEFAULT_INT .OR. ne .le. 0) call die("NE in COULOMB_NAMELIST must be defined and positive")
      if(lmax_partial .eq. DEFAULT_INT) call die("LMAX_PARTIAL in COULOMB_NAMELIST must be defined and nonnegative")
      if(lmax_total .eq. DEFAULT_INT .AND. (analytic_total_cb(1) .eqv. .false.)) &
        call die("LMAX_TOTAL in COULOMB_NAMELIST must be defined and nonnegative if ANALYTIC_TOTAL_CB(1) is .false.")
      if(do_xtrap) then
        if(ne_xtrap .eq. DEFAULT_INT .OR. ne_xtrap .le. 0) then
          call die("NE_XTRAP in COULOMB_NAMELIST must be defined and positive if DO_XTRAP is .TRUE.")
        endif
        if(Ei_xtrap .eq. 0.0_dp .OR. Ei_xtrap .gt. Ei) then
          call die("EI_XTRAP in COULOMB_NAMELIST must be defined, nonzero, and less than EI if DO_XTRAP is .TRUE.")
        endif
      endif
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call to_lower(zaxis)

    ! -- convert to atomic units
    Ei                = Ei                / au2ev
    Ef                = Ef                / au2ev
    Ei_xtrap          = Ei_xtrap          / au2ev
    ABC(:)            = ABC(:)            / au2invcm
    xs_zero_threshold = xs_zero_threshold / (au2cm*au2cm)

    ! -- convert to lower case
    call to_lower(rotor_kind)
    call to_lower(point_group)
    call to_lower(kmat_energy_units_override)
    call to_lower(channel_energy_units_override)

    ! -- remove spaces
    point_group      = trim(point_group)
    kmat_dir         = trim(kmat_dir)
    channels_dir     = trim(channels_dir)
    output_directory = trim(output_directory)

    ! -- add trailing directory separator to directories if needed, make directories as needed
    call add_trailing(output_directory, ds)
    call add_trailing(kmat_dir,         ds)
    call add_trailing(channels_dir,     ds)

    write(stdout, '(A)') "--------------------------------------------------------------------------------------------------------"
    write(stdout, *)
    write(stdout, control_namelist)
    write(stdout, *)
    if(use_kmat .eqv. .true.) then
      write(stdout, kmat_namelist)
      write(stdout, *)
    endif
    if(use_CB .eqv. .true.) then
      write(stdout, coulomb_namelist)
      write(stdout, *)
    endif
    write(stdout, '(A)') "--------------------------------------------------------------------------------------------------------"
    write(stdout, *)

    ! -- checks
    if(Nmin .gt. Nmax) call die("Nmin > Nmax not allowed")
    if(target_charge .eq. DEFAULT_INT) call die("Must set the charge of the target in namelist CONTROL !")
    if(target_charge .eq. 0) call die("Neutral targets not programmed yet !")

    ! -- namelist: control
    cfg%nmin                     = nmin
    cfg%nmax                     = nmax
    cfg%use_kmat                 = use_kmat
    cfg%use_cb                   = use_cb
    cfg%spin_isomer_kind         = spin_isomer_kind
    cfg%output_directory         = output_directory
    cfg%rotor_kind               = rotor_kind
    cfg%zaxis                    = zaxis
    cfg%abc                      = abc(:)
    cfg%target_charge            = target_charge
    cfg%add_cd4                  = add_cd4
    cfg%add_cd6                  = add_cd6
    cfg%xs_zero_threshold        = xs_zero_threshold
    if(add_cd4 .eqv. .true.) then
      if(any([dn,dnk,dk,deltan,deltak] .eq. 0.0_dp)) &
        call die("Please specify all 4th-order centrifugal distortion parameters !")
      dn      = dn     / au2invcm
      dnk     = dnk    / au2invcm
      dk      = dk     / au2invcm
      deltan  = deltan / au2invcm
      deltak  = deltak / au2invcm
      cfg%cd4 = cd4_type(dn = dn, dnk = dnk, dk = dk, deltan = deltan, deltak = deltak)
    endif
    if(add_cd6 .eqv. .true.) then
      if(add_cd4 .eqv. .false.) call die("Don't add the sextic correction while omitting the quartic correction !")
      if(any([hn,hnk,hkn,hk,etan,etank,etak] .eq. 0.0_dp)) &
        call die("Please specify all 6th-order centrifugal distortion parameters !")
      hn    = hn    / au2invcm
      hnk   = hnk   / au2invcm
      hkn   = hkn   / au2invcm
      hk    = hk    / au2invcm
      etan  = etan  / au2invcm
      etank = etank / au2invcm
      etak  = etak  / au2invcm
      cfg%cd6 = cd6_type(hn = hn, hnk = hnk, hkn = hkn, hk = hk, etan = etan, etank = etank, etak = etak)
    endif

    ! -- namelist: kmat
    if(use_kmat .eqv. .true.) then
      if(kmat_output_type .eq. "=======") then
        call die("Must speficy KMAT_OUTPUT_TYPE. It should be one of "// UKRMOLX //" or "// MQDTR2K)
      elseif(all(kmat_output_type .ne. [UKRMOLX, MQDTR2K])) then
        call die("Poorly specified KMAT_OUTPUT_TYPE. It should be one of "// UKRMOLX //" or "// MQDTR2K)
      endif
      cfg%kmat_dir                      = kmat_dir
      cfg%channels_dir                  = channels_dir
      cfg%lmax_kmat                     = lmax_kmat
      cfg%point_group                   = point_group
      cfg%spinmults                     = spinmults(:)
      cfg%num_egrid_segs                = num_egrid_segs
      cfg%num_egrid                     = num_egrid(:)
      cfg%egrid_segs                    = egrid_segs(:)
      cfg%egrid_spacing                 = egrid_spacing
      cfg%real_spherical_harmonics      = real_spherical_harmonics
      cfg%kmat_energy_closest           = kmat_energy_closest / au2ev
      cfg%kmat_output_type              = kmat_output_type
      cfg%kmat_energy_units_override    = kmat_energy_units_override
      cfg%channel_energy_units_override = channel_energy_units_override
    endif

    ! -- namelist: coulomb
    if(use_CB .eqv. .true.) then
      if(do_quadrupole) call die("DO_QUADRUPOLE exists as an option, but I'm yet confident in its&
        & implementation. Remove this call if you want and see what happens.")
#ifndef USE_CDMSREADER
      if(use_cdms_einsta .eqv. .true.) call die("User requested use of CDMS data, but the code is&
        & not compiled with that capability. Build with 'USE_CDMSREADER=1' to change this.")
#endif
      cfg%use_cdms_einsta              = use_cdms_einsta
      cfg%analytic_total_cb            = analytic_total_cb(:)
      cfg%ei                           = ei
      cfg%ef                           = ef
      cfg%ne                           = ne
      cfg%ne_xtrap                     = ne_xtrap
      cfg%ei_xtrap                     = ei_xtrap
      cfg%do_xtrap                     = do_xtrap
      cfg%do_dipole                    = do_dipole
      cfg%do_quadrupole                = do_quadrupole
      cfg%lmax_partial                 = lmax_partial
      cfg%lmax_total                   = lmax_total
      cfg%cartesian_dipole_moments     = cartesian_dipole_moments(:)     / au2deb
      ! cfg%cartesian_quadrupole_moments = cartesian_quadrupole_moments(:) !/ au2deb
      cfg%cdms_file                    = cdms_file
      cfg%only_einsta                  = only_einsta
    endif

  end subroutine read_namelists

! ================================================================================================================================ !
end module rotex__reading
! ================================================================================================================================ !
