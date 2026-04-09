! ================================================================================================================================ !
module rotex__reading
  !! Contains procedures used in reading data (K-matrices and namelist data)
  use rotex__globals,   only: G, UKRMOLX, MQDTR2K
  use rotex__system, only: stdout, stderr, die

  implicit none (type, external)

  private

  public :: read_kmats

! ================================================================================================================================ !
contains
! ================================================================================================================================ !


  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine read_kmats(                    &
                                kmat               &
                              , kmat_eval_energies &
                              , spinmult           &
                              , elec_channels      &
                              , channel_E_units    &
                              , kmat_eval_E_units  &
    )
    !! Reads in a K-matrix from a file with a very particular file format given by kmat_output_type

    use rotex__kinds,      only: dp
    use rotex__types,      only: elec_channel_type, rmatrix_type, ivector_type, r3rarr_type
    use rotex__channel_ops, only: permsort_channels
    use rotex__utils,      only: read_blank
    use rotex__arrays,     only: append, is_symmetric, realloc
    use rotex__system,     only: die, stdout, stderr, IOSTAT_END, IOSTAT_OK
    use rotex__symmetry,   only: group_size, irrep_name
    use rotex__constants,  only: au2ev
    use rotex__globals,    only: spinmult_names, DEFAULT_INT
    use rotex__characters, only: int2char

    implicit none (type, external)

    real(dp), intent(out), allocatable :: Kmat(:,:,:)
      !! K-matrix: nchan × nchan × ne
    real(dp), intent(out), allocatable :: kmat_eval_energies(:)
      !! Evaluation energies of the K-matrix
    integer, intent(in) :: spinmult
      !! The current spin multiplicity
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

    logical :: skip_this_irrep
    integer :: i, j, ichan, i1, i2
    integer :: irrep
    integer :: nirreps
    integer :: ne
    integer :: nchans_total
    integer, allocatable :: nchans_irrep(:), idx(:)
    character(:), allocatable :: kmat_filename, channels_filename, filename, irrepname

    real(dp) :: channel_e_convert, kmat_e_convert
    real(dp), allocatable :: kmat_eval_energies_local(:)

    type(elec_channel_type), allocatable :: elec_channels_this_irrep(:)
    type(ivector_type), allocatable :: index_map(:)
    type(r3rarr_type), allocatable :: kmat_irrep(:)

    nirreps = group_size(G%POINT_GROUP)
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
    allocate(kmat_irrep(nirreps))

    nchans_total = 0

    write(stdout, '("Point group: ", A)') G%POINT_GROUP

    ! -- read the K-matrices and electronic channels
    irrep_loop_kmats: do irrep = 1, nirreps
      irrepname = irrep_name(irrep, G%POINT_GROUP)
      write(stdout, '(2X, "Irrep: ", A)') irrepname
      filename = G%KMAT_DIR // int2char(spinmult) // irrepname // ".kmat"

      select case(G%KMAT_OUTPUT_TYPE)
      case(UKRMOLX)

        channels_filename = G%CHANNELS_DIR // "channels.geom1." // spinmult_names(spinmult) // "." // irrepname
        kmat_filename     = G%KMAT_DIR     // "K-matrix.geom1." // spinmult_names(spinmult) // "." // irrepname

        call get_kmat_and_channels_ukrmolx( &
            channels_filename               &
          , kmat_filename                   &
          , kmat_irrep(irrep)%r3arr         &
          , kmat_eval_energies_local        &
          , channel_E_convert               &
          , kmat_e_convert                  &
          , elec_channels_this_irrep        &
          , nchans_irrep(irrep)             &
          , skip_this_irrep)

      case(MQDTR2K)

        ! -- the kmat file is expected to have the channels in this format, so no need for a channels file
        kmat_filename = G%KMAT_DIR // int2char(spinmult) // irrep_name(irrep, G%POINT_GROUP) // ".kmat"

        call get_kmat_and_channels_mqdtr2k( &
            kmat_filename                   &
          , kmat_irrep(irrep)%r3arr     &
          , kmat_eval_energies_local        &
          , channel_E_convert               &
          , kmat_e_convert                  &
          , elec_channels_this_irrep        &
          , nchans_irrep(irrep)             &
          , skip_this_irrep                 &
        )

      case default
        call die("KMAT_OUTPUT_TYPE must be "//UKRMOLX//" or "//MQDTR2K)
      end select

      ! ! -- if we need to focus on channel parity
      ! if(G%POINT_GROUP .eq. "cs") call fill_parity_array_this_irrep_cs(G%LMAX_KMAT, elec_channels_this_irrep, irrep, G%POINT_GROUP)

      call append(elec_channels, elec_channels_this_irrep)

      if(skip_this_irrep .eqv. .true.) cycle irrep_loop_kmats

      ! -- total number of channels across all irreps
      nchans_total = nchans_total + nchans_irrep(irrep)
      ne = size(kmat_irrep(irrep)%r3arr, 3)

      ! -- index map: this irrep → full basis of channels
      allocate(index_map(irrep)%vec(nchans_irrep(irrep)))
      do concurrent( ichan=1:nchans_irrep(irrep) )
        index_map(irrep)%vec(ichan) = find_global_channel_index(elec_channels_this_irrep(ichan), elec_channels)
      enddo

    enddo irrep_loop_kmats

    call move_alloc(kmat_eval_energies_local, kmat_eval_energies)

    deallocate(elec_channels_this_irrep)

    ! -- K-matrix per irrep -> total K-matrix
    allocate(kmat(nchans_total, nchans_total, ne), source=0.0_dp)

    call fill_total_kmat(kmat_irrep, index_map, nchans_irrep, nchans_total, ne, kmat)

    deallocate(kmat_irrep)
    deallocate(index_map)

    ! -- this is probably not necessary, but it's a nice order in which to have channels
    call permsort_channels(elec_channels, idx)
    kmat(:,:,:) = kmat(idx,idx,:)

    if(G%print_elec_channels) call print_channels(elec_channels, stdout)

  end subroutine read_kmats

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine fill_total_kmat( &
        kmat_irrep                 &
      , index_map                  &
      , nchans_irrep               &
      , nchans_total               &
      , ne                         &
      , kmat                       &
    )
    !! Using index_map, put the elements of kmat_irrep into the full kmat
    use rotex__kinds,    only: dp
    use rotex__types,    only: r3rarr_type, ivector_type
    use rotex__arrays,   only: ij2k, size_check
    use rotex__symmetry, only: group_size
    implicit none (type, external)
    type(r3rarr_type), intent(in)   :: kmat_irrep(:)
      !! Array of K-matrices for each irrep
    type(ivector_type), intent(in)  :: index_map(:)
      !! Irrep channel -> total channel index map
    integer,            intent(in)  :: nchans_irrep(:)
      !! Number of channels per irrep
    integer,            intent(in)  :: nchans_total
      !! Number of channels across all irreps
    integer,            intent(in)  :: ne
      !! Number of K-matrix evaluation energies
    real(dp),           intent(out) :: kmat(:,:,:)
      !! The K-matrices for all energies. nchan × nchan × nE
    integer :: iloc, jloc, kloc, itot, jtot, ktot, irrep, nirreps, ie
    call size_check(kmat, [nchans_total, nchans_total, ne], "KMAT")
    nirreps = group_size(G%POINT_GROUP)
    do irrep=1,nirreps
      kloc = 0
      do concurrent(iloc=1:nchans_irrep(irrep), jloc=1:nchans_irrep(irrep), ie=1:ne)
        itot = index_map(irrep)%vec(iloc)
        jtot = index_map(irrep)%vec(jloc)
        kmat(itot, jtot, ie) = kmat_irrep(irrep)%r3arr(iloc, jloc, ie)
      enddo
    enddo
  end subroutine fill_total_kmat

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  ! subroutine fill_parity_array_this_irrep_cs(lmax_kmat, elec_channels, irrep, point_group)
  !   !! Fill the array M_PARITY, stored in the module ROTEX__SYMMETRY, that will later be accessed
  !   !! by the rotational frame transformation to determine which values of m from -lmax_kmat to lmax_kmat
  !   !! correspond to even and odd combinations of partial waves. The irrep corresponds to total irrep,
  !   !! so we need to factor out the irrep of the electronic state. This routine
  !   !! is only expected to be called for Cs symmetry
  !   use rotex__types,      only: elec_channel_type
  !   use rotex__system,     only: die
  !   use rotex__symmetry,   only: Ap, App, m_parity, even, odd, elecstate_parity_set, elecstate_parity
  !   use rotex__characters, only: i2c => int2char
  !   implicit none (type, external)
  !   integer, intent(in) :: lmax_kmat
  !     !! Max value of l for the electronic channels
  !   type(elec_channel_type), intent(in) :: elec_channels(:)
  !     !! The electronic channels for the current irrep
  !   integer, intent(in) :: irrep
  !     !! The current irrep index
  !   character(*), intent(in) :: point_group
  !     !! The point group for the scattering calculations
  !   integer :: ichan, m
  !   integer :: gs_parity
  !     !! Parity of the ground state
  !   if(point_group .ne. "cs") call die("Attempting to fill m_parity array for a point group&
  !     & other than Cs: " // point_group)
  !   if(allocated(m_parity) .eqv. .false.) allocate(m_parity(-lmax_kmat:lmax_kmat), source = 0)
  !   ! -- determine the parity of the electronic state
  !   if(elecstate_parity_set .eqv. .false.) then
  !     if(any(elec_channels%ml .eq. 0) .eqv. .true.) then
  !       ! -- m=0 behaves as A' (even) and therefore the parity of the electronic state
  !       !    will be the parity of the total channel if m=0 is included:
  !       !    Γtot = Γelec × Γm = Γelec × A' = Γelec
  !       elecstate_parity = merge(even, odd, irrep .eq. Ap)
  !       elecstate_parity_set = .true.
  !     else
  !       ! -- having no m=0 (A') channels means that this electronic state has
  !       !    the opposite parity of total channel:
  !       !      Γtot = Γelec × Γm = Γelec × A''
  !       !    If Γtot is A' (even), Γelec is A'' (odd) and vice versa
  !       elecstate_parity = merge(odd, even, irrep .eq. Ap)
  !       elecstate_parity_set = .true.
  !     endif
  !   endif
  !   do ichan=1, size(elec_channels, 1)
  !     m = elec_channels(ichan) % ml
  !     if(m_parity(m) .ne. 0) cycle ! skip if set, but this probably should not happen
  !     select case(irrep)
  !     case(Ap)
  !       m_parity(m) = even
  !     case(App)
  !       m_parity(m) = odd
  !     case default
  !       call die("Somehow, irrep ("//i2c(irrep)//") is not one of the valid values for the point group " // point_group)
  !     end select
  !   enddo
  ! end subroutine fill_parity_array_this_irrep_cs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine get_kmat_and_channels_ukrmolx( &
        channels_filename                        &
      , kmat_filename                            &
      , kmat                                     &
      , kmat_energies                            &
      , channel_e_convert                        &
      , kmat_e_convert                           &
      , elec_channels_this_irrep                 &
      , nchans_this_irrep                        &
      , skip_this_irrep)
    !! Return the K-matrix that is closest to the desired evaluation energy
    !! given by G%KMAT_ENERGY_CLOSEST
    use rotex__constants, only: au2ev
    use rotex__kinds,     only: dp
    use rotex__arrays,    only: realloc, size_check, append, unpackmat
    use rotex__types,     only: elec_channel_type
    use rotex__utils,     only: read_blank
    use rotex__system,    only: stdout, stderr, die, IOSTAT_OK, IOSTAT_END

    implicit none (type, external)

    character(*), intent(in) :: channels_filename
      !! Where to read channels
    character(*), intent(in) :: kmat_filename
      !! Where to read K-matrices
    real(dp), intent(out), allocatable :: kmat(:,:,:)
      !! K-matrix for each energy: nchan × nchan × nE_include
    real(dp), allocatable, intent(inout) :: kmat_energies(:)
      !! Array of K-matrix evaluation energies that we want
    real(dp), intent(in) :: channel_e_convert
      !! Energy type of channel energy. Should be replaced by a global
    real(dp), intent(in) :: kmat_e_convert
      !! Energy type of K-matrix evaluation energy. Should be replaced by a global
    type(elec_channel_type), intent(out), allocatable :: elec_channels_this_irrep(:)
      !! Channels for this irrep
    integer, intent(out) :: nchans_this_irrep
      !! Number of channels for this irrep
    logical, intent(out) :: skip_this_irrep
      !! Whether to skip this irrep (e.g., no channels)

    integer,      parameter :: UKRMOL_KMAT_ELEMENTS_PER_LINE = 4
    character(8), parameter :: UKRMOL_KMAT_ELEMENTS_FMT = '(D20.13)'

    integer  :: iostat, ne, ne_inlcude, i, ichan, ie, ie_closest, l, ml, nelec, nchans, iflat
    integer  :: ie_include, ne_include, iemin, iemax
    integer  :: funit, nchans_max, nskip, nskip_header, iline, icol, nlines, nchans_flat
    real(dp) :: E
    real(dp), allocatable :: kmat_flat(:, :)
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

    ! -- total number of evaluation energies in file
    ne = 0

    ! -- initial number of evaluation energies to retain for the K-matrix
    !    0: count up for EDFT
    !    1: does not change for non-EDFT
    ne_include = merge(0, 1, G%EDFT)

    call read_blank(funit, 3)
    read(funit, *) i, i, i, nchans_max, nskip_header
    call read_blank(funit, nskip_header)
    nskip_header = nskip_header + 4 ! for the next re-reads

    ! -- count number of energies and the number of energies that we want to include if EDFT
    do
      read(funit, *, iostat = iostat) nchans, i, nchans_flat, E
      if(iostat .eq. IOSTAT_END) exit
      E = E*kmat_e_convert ! store in au
      ne = ne + 1
      nskip = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind = dp))
      call read_blank(funit, nskip)
      ! -- don't worry about ne_include if we just want one
      if(G%EDFT .eqv. .false.) cycle
      if(E .lt. G%KMAT_EI .OR. (E .gt. G%KMAT_EF .AND. G%KMAT_EF .gt. 0._dp)) cycle
      ne_include = ne_include + 1
    enddo

    call realloc(kmat_flat, nchans_flat, ne_include)
    kmat_flat = 0._dp

    ! -- read in the K-matrix energies (including those we won't use)
    rewind(funit)
    call realloc(kmat_energies, ne) ! this will be trimmed down to included energies later
    call read_blank(funit, nskip_header)
    ie = 0
    nchans_this_irrep = 0
    do
      read(funit, *, iostat = iostat) nchans, i, nchans_flat, E
      if(iostat .eq. IOSTAT_END) exit
      ie = ie + 1
      E = E*kmat_e_convert ! store in au
      kmat_energies(ie) = E
      nskip = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind=dp))
      call read_blank(funit, nskip)
      ! -- take the smallest valid K-matrix and set its number of channels as the number of channels
      if(nchans_this_irrep .ne. 0) cycle
      ! -- take first matrix satisfying E > KMAT_EI if EDFT
      if(G%EDFT) then
        if(E .lt. G%KMAT_EI) cycle
      endif
      ! -- just take the first matrix otherwise
      nchans_this_irrep = nchans
    enddo

    ! -- inform user of energy selection
    if(G%EDFT) then
      ! -- minimum and maximum evaluation energies to consider for EDFT
      iemin = findloc(kmat_energies .ge. G%KMAT_EI, .true., 1)
      iemax = merge(ne, findloc(kmat_energies .le. G%KMAT_EF, .true., 1), G%KMAT_EF .eq. 0._dp)
      write(stdout, '(4x, "user requested k-matrix energies between ", es10.3, " and ", es10.3, " ev")') &
        kmat_energies(iemin)*au2ev, kmat_energies(iemax)*au2ev
    else
      ! -- find the lowest energy if energy independent
      ie_closest = minloc(abs(kmat_energies - G%KMAT_ENERGY_CLOSEST), 1)
      iemin = ie_closest
      iemax = ie_closest
      if(ie_closest .lt. 1) call die("Somehow, IE_CLOSEST returned a non-positive integer !")
      write(stdout, '(4X, "User requested K-matrix at ", E20.10, " eV")') G%KMAT_ENERGY_CLOSEST          * au2ev
      write(stdout, '(7X, "Found Kmatrix at energy ",    E20.10, " eV")') kmat_energies(ie_closest) * au2ev
    endif

    ! -- Now, actually go and read that (those) K-matrix (K-matrices)
    rewind(funit)
    call read_blank(funit, nskip_header)
    ie = 0
    ie_include = 0

    kmat_read: if(G%EDFT) then

      ! -- energy dependent read
      do

        read(funit, *, iostat = iostat) nchans, i, nchans_flat, E
        if(iostat .eq. IOSTAT_END) exit kmat_read

        ie = ie + 1
        ! -- skip this K-matrix if it's too low in E
        if(ie .lt. iemin) then
          nskip = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind=dp))
          call read_blank(funit, nskip)
          cycle
        endif
        ! -- skip the remaining K-matrices if we've reached our max energy
        if(ie .gt. iemax) exit kmat_read

        ! -- at this point, we have found a K-matrix that we want to add.
        ie_include = ie_include + 1
        ! -- iterate through the lines and columns of the flattened K-matrix
        !    in the file
        iflat = 0
        nlines = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind=dp))
        do iline = 1, nlines
          do icol = 1, UKRMOL_KMAT_ELEMENTS_PER_LINE
            iflat = iflat + 1

            ! -- ignore the rest of the elements if the number of channels has increased
            if(iflat .gt. size(kmat_flat, 1)) cycle

            read(funit, UKRMOL_KMAT_ELEMENTS_FMT, advance = 'no') kmat_flat(iflat, ie_include)
            if(iflat .eq. nchans_flat) cycle
          enddo
          read(funit, *)
        enddo

      enddo

    else

      ! -- energy independent read
      do

        read(funit, *, iostat = iostat) nchans, i, nchans_flat, E
        if(iostat .eq. IOSTAT_END) then
          write(stderr, '("Number of K-matrices/energies: ", I0)') ne
          write(stderr, '("Target K-matrix energy: ", E20.10)') G%KMAT_ENERGY_CLOSEST * au2ev
          write(stderr, '("Closest available K-matrix is number ", I0)') ie_closest
          call die("Could not find the K-matrix that is closest to the given target energy before EOF")
        endif

        ie = ie + 1

        ! -- find the corresponding K-matrix
        if(ie .ne. ie_closest) then
          ! -- skip
          nskip = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind=dp))
          call read_blank(funit, nskip)
        else
          ! -- found it
          nlines = ceiling(nchans_flat / real(UKRMOL_KMAT_ELEMENTS_PER_LINE, kind=dp))
          iflat = 0
          ! -- iterate through the lines and columns of the flattened K-matrix
          !    in the file
          do iline = 1, nlines
            do icol = 1, UKRMOL_KMAT_ELEMENTS_PER_LINE
              iflat = iflat + 1
              read(funit, UKRMOL_KMAT_ELEMENTS_FMT, advance = 'no') kmat_flat(iflat, 1)
              if(iflat .eq. nchans_flat) exit kmat_read
            enddo
            read(funit, *)
          enddo

          write(stderr, '("NLINES: ", I0)') nlines
          write(stderr, '("UKRMOL_KMAT_ELEMENTS_PER_LINE: ", I0)') UKRMOL_KMAT_ELEMENTS_PER_LINE
          write(stderr, '("NCHANS: ", I0)') nchans
          write(stderr, '("NCHANS_FLAT: ", I0)') nchans_flat

          call die("Improper K-matrix read loop exit. ")

        endif
      enddo
    endif kmat_read

    close(funit)

    ! -- filter K-matrix energies to keep only those that we want to include
    kmat_energies = kmat_energies(iemin:iemax)

    if(nchans_this_irrep .lt. 1) then
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

    ! -- flat Kmat -> Kmat
    allocate(Kmat(nchans_this_irrep, nchans_this_irrep, ne_include), source=0.0_dp)
    do ie=1, ne_include
      call unpackmat(kmat_flat(:,ie), kmat(:,:,ie))
    enddo

  end subroutine get_kmat_and_channels_ukrmolx

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine get_kmat_and_channels_mqdtr2k( &
        kmat_filename                       &
      , kmat                                &
      , kmat_energies                       &
      , channel_e_convert                   &
      , kmat_e_convert                      &
      , elec_channels_this_irrep            &
      , nchans_this_irrep                   &
      , skip_this_irrep)
    !! Return the K-matrix that is closest to the desired evaluation energy
    !! given by G%KMAT_ENERGY_CLOSEST where the K-matrix channels format is that of MQDTR2K

    use rotex__constants, only: au2ev
    use rotex__kinds,     only: dp
    use rotex__arrays,    only: realloc, size_check, unpackmat
    use rotex__types,     only: elec_channel_type
    use rotex__utils,     only: read_blank
    use rotex__system,    only: stdout, die, IOSTAT_END, IOSTAT_OK

    implicit none (type, external)

    character(*), intent(in) :: kmat_filename
      !! Where to read the channels and K-matrices
    real(dp), intent(out), allocatable :: kmat(:,:,:)
      !! K-matrices: nchan × nchan × nE_include
    real(dp), intent(out), allocatable :: kmat_energies(:)
      !! Array of K-matrix evaluation energies to include
    real(dp), intent(in) :: channel_e_convert
      !! Channel energy conversion units. Should be replaced by a global maybe
    real(dp), intent(in) :: kmat_e_convert
      !! K-matrix evaluation energy conversion units. Should be replaced by a global maybe
    type(elec_channel_type), intent(out), allocatable :: elec_channels_this_irrep(:)
      !! Array of electronic channels for this irrep
    integer, intent(out) :: nchans_this_irrep
      !! Number of electronic channels this irrep
    logical, intent(out) :: skip_this_irrep
      !! Whether to skip this irrep (e.g., no channels for this irrep)

    integer  :: iostat, ne, funit, ne_skip, i, ichan, ie, ie_closest, l, ml, iq, nelec
    integer  :: ie_include, ne_include, iemin, iemax
    integer  :: nchans_irrep_flat
    real(dp) :: E
    real(dp), allocatable :: kmat_flat(:,:)

    skip_this_irrep = .false.

    inquire(file = kmat_filename, iostat = iostat)
    if(iostat .ne. IOSTAT_OK) then
      skip_this_irrep = .true.
      return
    endif

    open(newunit = funit, file = kmat_filename)

    ! -- get number of channels, skip to K-matrices, determine number of energies/K-matrices to read
    ne = 0
    ne_include = 0
    call read_blank(funit)
    read(funit, *) nchans_this_irrep
    if(nchans_this_irrep .lt. 1) then
      write(stderr, '("NCHANS_THIS_IRREP: I0")') nchans_this_irrep
      call die("READ <1 channels for this irrep")
    endif
    nchans_irrep_flat = (nchans_this_irrep * (nchans_this_irrep+1)) / 2
    call read_blank(funit, nchans_this_irrep)
    do
      ! -- energies are the first number, ignore the rest
      read(funit, *, iostat = iostat) E
      if(iostat .eq. IOSTAT_END) exit
      ne = ne + 1
      if(G%EDFT .eqv. .false.) cycle
      ! -- if EDFT, keep track of number of energies to include
      E = E*kmat_e_convert
      if(E .lt. G%KMAT_EI) cycle
      if(E .gt. G%KMAT_EF .AND. G%KMAT_EF .gt. 0._dp) cycle
      ne_include = ne_include + 1
    enddo

    if(G%EDFT) then
      if(ne_include .lt. 1) then
        write(stderr, '("NE_INCLUDE: ", I0)') NE_INCLUDE
        call die("NE_INCLUDE < 1")
      endif
    else
      ne_include = 1
    endif

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

    ! -- min and max energies to consider, inform user of energy selection
    if(G%EDFT) then
      iemin = findloc(kmat_energies .ge. G%KMAT_EI, .true., 1)
      iemax = merge(ne, findloc(kmat_energies .le. G%KMAT_EF, .true., 1), G%KMAT_EF .eq. 0._dp)
      write(stdout, '(4x, "user requested k-matrix energies between ", es10.3, " and ", es10.3, " ev")') &
        kmat_energies(iemin)*au2ev, kmat_energies(iemax)*au2ev
      ne_skip = iemin - 1
    else
      iemin = 1
      iemax = ne
      ! -- find the closest energy
      ie_closest = minloc(abs(kmat_energies - G%KMAT_ENERGY_CLOSEST), 1)
      if(ie_closest .lt. 1) call die("Somehow, IE_CLOSEST returned a non-positive integer !")
      ne_skip = ie_closest - 1
      write(stdout, '(4X, "User requested K-matrix at ", E20.10, " eV")') G%KMAT_ENERGY_CLOSEST            * au2ev
      write(stdout, '(7X, "Found Kmatrix at energy ",    E20.10, " eV")') kmat_energies(ie_closest) * au2ev
    endif

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
    call realloc(kmat_flat, nchans_irrep_flat, ne_include)
    kmat_flat = 0.0_dp

    ! -- read the K-matrix/K-matrices
    if(G%EDFT) then
      do ie_include = 1, ne_include
        read(funit, *) E, (kmat_flat(i, ie_include), i=1, nchans_irrep_flat)
      enddo
    else
      read(funit, *) E, (kmat_flat(i, 1), i=1, nchans_irrep_flat)
    endif
    close(funit)

    ! -- flat Kmat -> Kmat
    allocate(kmat(nchans_this_irrep, nchans_this_irrep, ne_include), source=0._dp)
    do concurrent (ie=1:ne_include)
      call unpackmat(kmat_flat(:,ie), kmat(:,:,ie))
    enddo

  end subroutine get_kmat_and_channels_mqdtr2k

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine print_channels(channels, funit)
    use rotex__kinds,  only: dp
    use rotex__types,  only: elec_channel_type
    use rotex__system, only: stdout
    implicit none (type, external)
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
    use rotex__types,  only: elec_channel_type
    use rotex__channel_ops,  only: operator(.ne.)
    use rotex__system, only: die
    implicit none (type, external)
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

! ================================================================================================================================ !
end module rotex__reading
! ================================================================================================================================ !
