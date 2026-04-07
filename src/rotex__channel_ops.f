! ================================================================================================================================ !
module rotex__channel_ops
  !! Various channel operators
  use rotex__kinds, only: dp
  use rotex__types, only: channel_type, elec_channel_type &
                        , asymtop_rot_channel_type        &
                        , asymtop_rot_channel_l_type      &
                        , asymtop_rot_transition_type

  implicit none (type, external)

  private

  public :: operator(.eq.)
  public :: operator(.ne.)
  public :: operator(.isin.)
  public :: assignment(=)
  public :: findloc_transitions
  public :: sort_channels
  public :: permsort_channels
  public :: sort_channels_by_energy
  public :: trim_channel_l
  public :: get_channel_index

  interface operator(.eq.)
    module procedure :: channel_iseq
    module procedure :: transition_iseq
  end interface operator(.eq.)

  interface operator(.ne.)
    module procedure :: channel_isne
  end interface operator(.ne.)

  interface operator(.isin.)
    module procedure channel_isin
    module procedure transition_isin
  end interface operator(.isin.)

  interface assignment(=)
    module procedure :: channel_set_eq
  end interface assignment(=)

  interface sort_channels
    module procedure :: sort_elec_channels
  end interface sort_channels

  interface permsort_channels
    module procedure :: permsort_elec_channels
  end interface permsort_channels

  interface findloc_transitions
    module procedure :: findloc_transitions_scl
    module procedure :: findloc_transitions_arr
  end interface findloc_transitions

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure elemental module subroutine channel_set_eq(channel_out, channel_in)
    !! Sets the channel left equal to the channel right
    use rotex__system, only: die
    implicit none (type, external)
    class(channel_type), intent(out) :: channel_out
    class(channel_type), intent(in)  :: channel_in

    ! -- start at base class
    channel_out % nelec = channel_in % nelec
    channel_out % E     = channel_in % E

    select type (left => channel_out)

    ! -- return if this is just a base channel
    class default
      return

    ! -- electronic channels
    type is (elec_channel_type)
      select type (right => channel_in)
      type is (elec_channel_type)
        left % l  = right % l
        left % ml = right % ml
        left % iq = right % iq
      class default
        call die("Trying to assign a non-electronic channel to an electronic channel.")
      end select

    ! -- rotational channels + l
    type is (asymtop_rot_channel_l_type)
      select type (right => channel_in)
      type is (asymtop_rot_channel_l_type)
        left % N   = right % N
        left % Ka  = right % Ka
        left % Kc  = right % Kc
        left % l   = right % l
        left % iq  = right % iq
        left % sym = right % sym
      class default
        call die("Trying to assign rotation+l channel to a different kind of channel.")
      end select

    ! -- rotational channels (no l)
    type is (asymtop_rot_channel_type)
      select type (right => channel_in)
      type is (asymtop_rot_channel_type)
        left % N   = right % N
        left % Ka  = right % Ka
        left % Kc  = right % Kc
        left % sym = right % sym
      class default
        call die("Trying to assign rotational channel (no l) to a different kind of channel.")
      end select

    end select
  end subroutine channel_set_eq

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module function channel_isin(channel, channels) result(res)
    !! Check if a channel is in the array channels
    implicit none (type, external)
    class(channel_type), intent(in) :: channel, channels(:)
    logical :: res
    res = .true.
    if(any(channel .eq. channels)) return
    res = .false.
  end function channel_isin
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module function transition_isin(transition, transitions) result(res)
    !! Check if a transition is in the array transitions
    implicit none (type, external)
    type(asymtop_rot_transition_type), intent(in) :: transition, transitions(:)
    logical :: res
    res = .true.
    if(any(transition .eq. transitions)) return
    res = .false.
  end function transition_isin


  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental module function transition_iseq(transition1, transition2) result(res)
    !! Test if two transitions are equal
    implicit none (type, external)
    type(asymtop_rot_transition_type), intent(in) :: transition1, transition2
    logical :: res
    res = .false.
    if(transition1 % lo .ne. transition2 % lo) return
    if(transition1 % up .ne. transition2 % up) return
    res = .true.
  end function transition_iseq

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental module function channel_iseq(channel1, channel2) result(res)
    !! Test for channel equality on the basis of their quantum numbers only.
    !!   Asymmetric tops: (N,Ka,Kc)
    !!   Symmetric tops: (N,0,Kc) or (N,Ka,0)
    !!   Linear tops: (N,0,0)
    !! This will also look at l,λ and sym if they are present; these do not depend on top kind
    use rotex__system,  only: die
    implicit none (type, external)
    class(channel_type), intent(in) :: channel1, channel2
    logical :: res
    integer :: Ksym1, Ksym2

    res = .false.

    if(channel1 % nelec .ne. channel2 % nelec)  return

    select type(channel1)

    ! -- electronic channel equality comparison
    type is (elec_channel_type)
      select type(channel2)
      type is (elec_channel_type)
        if(channel1 % l  .ne. channel2 % l)   return
        if(channel1 % ml .ne.  channel2 % ml) return
      class default
        call die("Cannot compare an electronic channel to a different channel")
      end select

    ! -- rotational channel equality comparison (with l)
    type is (asymtop_rot_channel_l_type)
      select type(channel2)
      type is (asymtop_rot_channel_l_type)
        if(channel1 % N  .ne. channel2 % N)  return
        if(channel1 % Ka .ne. channel2 % Ka) return
        if(channel1 % Kc .ne. channel2 % Kc) return
        if(channel1 % l  .ne. channel2 % l)  return
        if(channel1 % sym.ne. channel2 % sym)return
      type is (asymtop_rot_channel_type)
        if(channel1 % N  .ne. channel2 % N)  return
        if(channel1 % Ka .ne. channel2 % Ka) return
        if(channel1 % Kc .ne. channel2 % Kc) return
        if(channel1 % sym.ne. channel2 % sym)return
      class default
        call die("Cannot compare a rotational channel (with l) to a non-rotational channel")
      end select

    ! -- rotational channel equality comparison (no l)
    type is (asymtop_rot_channel_type)
      select type(channel2)
      class is (asymtop_rot_channel_type)
        if(channel1 % N  .ne. channel2 % N)  return
        if(channel1 % Ka .ne. channel2 % Ka) return
        if(channel1 % Kc .ne. channel2 % Kc) return
        if(channel1 % sym.ne. channel2 % sym)return
      class default
        call die("Cannot compare a rotational channel (no l) to a non-rotational channel")
      end select

    end select

    res = .true.

  end function channel_iseq
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental module function channel_isne(channel1, channel2) result(res)
    implicit none (type, external)
    class(channel_type), intent(in) :: channel1, channel2
    logical :: res
    res = .not. (channel1 .eq. channel2)
  end function channel_isne

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module function findloc_transitions_scl(targ, search) result(idxtarg)
    !! Find the first index for targs that map to the elements in search.
    !! Return 0 if there is no such mapping.
    implicit none (type, external)
    type(asymtop_rot_transition_type), intent(in) :: targ, search(:)
    integer :: idxtarg
      !! TARG -> SEARCH mapping
    integer, parameter :: IDX_NOT_FOUND = 0
    integer :: isearch, itarg, nsearch
    integer, allocatable :: idxsearch(:)
    logical, allocatable :: mask(:)
    nsearch = size(search, 1)
    idxtarg = IDX_NOT_FOUND
    mask = search .eq. targ
    idxsearch = pack([(isearch, isearch=1, nsearch)], mask)
    idxtarg = idxsearch(1) ! take first match
  end function findloc_transitions_scl
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module function findloc_transitions_arr(targs, search) result(idxtarg)
    !! Find the indices for each element in targs that map to the elements in search.
    !! Return 0 if there is no such mapping.
    implicit none (type, external)
    type(asymtop_rot_transition_type), intent(in) :: targs(:), search(:)
    integer, allocatable :: idxtarg(:)
      !! TARGS -> SEARCH mapping
    integer, parameter :: IDX_NOT_FOUND = 0
    integer :: isearch, itarg, ntargs, nsearch
    integer, allocatable :: idxsearch(:)
    logical, allocatable :: mask(:)
    ntargs  = size(targs, 1)
    nsearch = size(search, 1)
    allocate(idxtarg(ntargs), source = IDX_NOT_FOUND)
    do itarg=1, ntargs
      mask = search .eq. targs(itarg)
      if(all(mask .eqv. .false.)) cycle ! cycle if no match
      idxsearch = pack([(isearch, isearch=1, nsearch)], mask)
      idxtarg(itarg) = idxsearch(1) ! take first match
    enddo
  end function findloc_transitions_arr

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure module function get_channel_index(channels, channel, reverse) result(i)
    !! Return the first index i where channel .eq. channels(i) is .true., or
    !! the last index is reverse is .true.
    use rotex__system, only: die
    implicit none (type, external)
    class(channel_type), intent(in) :: channel, channels(:)
    logical, intent(in), optional :: reverse
    integer :: i, istart, iend, istep
    logical :: reverse_local
    reverse_local = .false. ; if(present(reverse)) reverse_local = reverse
    if(reverse_local .eqv. .true.) then
      istart = lbound(channels, 1)
      iend   = ubound(channels, 1)
      istep  = 1
    else
      istart = ubound(channels, 1)
      iend   = lbound(channels, 1)
      istep  = -1
    endif
    do i = istart, iend, istep
      if(channel .ne. channels(i)) cycle
      return
    enddo
    call die("Failed finding a channel in the channel array")
  end function get_channel_index

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure elemental module function trim_channel_l(channel_with_l) result(channel_without_l)
    !! Given a rotational channel with the l quantum number, return the equvalent channel without the l quantum number
    implicit none (type, external)
    type(asymtop_rot_channel_l_type), intent(in) :: channel_with_l
    type(asymtop_rot_channel_type) :: channel_without_l
    integer :: nelec, N, Ka, Kc, sym
    real(dp) :: E
    nelec = channel_with_l % nelec
    N     = channel_with_l % N
    Ka    = channel_with_l % Ka
    Kc    = channel_with_l % Kc
    sym   = channel_with_l % sym
    E     = channel_with_l % E
    channel_without_l = asymtop_rot_channel_type(nelec = nelec, N = N, Ka = Ka, Kc = Kc, E = E, sym = sym)
  end function trim_channel_l

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine permsort_elec_channels(elec_channels, idx)
    !! Sorts the inout array based on the quantum numbers of the channels, and returns
    !! the permutation array that would produce the same output
    implicit none (type, external)
    type(elec_channel_type), intent(inout) :: elec_channels(:)
    integer, allocatable :: idx(:)
    integer :: ichan, jchan, nchan
    nchan = size(elec_channels, 1)
    idx = [(ichan,ichan=1,nchan)]
    ! -- sort by nelec
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ge. elec_channels(ichan) % nelec) cycle
      call swap_channels(elec_channels, ichan, jchan)
      call swap_ints(idx, ichan, jchan)
    enddo ; enddo
    ! -- sort by l
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ne. elec_channels(ichan) % nelec) cycle
      if(elec_channels(jchan) % l .ge. elec_channels(ichan) % l) cycle
      call swap_channels(elec_channels, ichan, jchan)
      call swap_ints(idx, ichan, jchan)
    enddo ; enddo
    ! -- sort by λ (ml)
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ne. elec_channels(ichan) % nelec) cycle
      if(elec_channels(jchan) % l .ne. elec_channels(ichan) % l) cycle
      if(elec_channels(jchan) % ml .ge. elec_channels(ichan) % ml) cycle
      call swap_channels(elec_channels, ichan, jchan)
      call swap_ints(idx, ichan, jchan)
    enddo ; enddo
  end subroutine permsort_elec_channels

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine sort_elec_channels(elec_channels)
    !! Sorts the inout array based on the quantum numbers of the channels
    implicit none (type, external)
    type(elec_channel_type), intent(inout) :: elec_channels(:)
    integer :: ichan, jchan, nchan
    nchan = size(elec_channels, 1)
    ! -- sort by nelec
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ge. elec_channels(ichan) % nelec) cycle
      call swap_channels(elec_channels, ichan, jchan)
    enddo ; enddo
    ! -- sort by l
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ne. elec_channels(ichan) % nelec) cycle
      if(elec_channels(jchan) % l .ge. elec_channels(ichan) % l) cycle
      call swap_channels(elec_channels, ichan, jchan)
    enddo ; enddo
    ! -- sort by λ (ml)
    do ichan = 1, nchan ; do jchan = ichan+1, nchan
      if(elec_channels(jchan) % nelec .ne. elec_channels(ichan) % nelec) cycle
      if(elec_channels(jchan) % l .ne. elec_channels(ichan) % l) cycle
      if(elec_channels(jchan) % ml .ge. elec_channels(ichan) % ml) cycle
      call swap_channels(elec_channels, ichan, jchan)
    enddo ; enddo
  end subroutine sort_elec_channels

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine sort_channels_by_energy(channels)
    !! Bubble sort the array of channels such that the channel energies are in ascending order
    implicit none (type, external)
    class(channel_type), intent(inout) :: channels(:)
    logical :: swapped
    integer :: nchans
    integer :: j, i
    nchans = size(channels, 1)
    do j=nchans-1, 1, -1
      swapped = .false.
      do i=1,j
        if(channels(i) % E .le. channels(i+1) % E) cycle
        call swap_channels(channels, i, i+1)
        swapped = .true.
      enddo
      if(.not. swapped) exit
    enddo
  end subroutine sort_channels_by_energy

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine swap_ints(arr, ichan, jchan)
    !! Swaps the array elements at arr(ichan) and arr(jchan)
    use rotex__system, only: die
    implicit none (type, external)
    integer, intent(inout) :: arr(:)
    integer, intent(in) :: ichan, jchan
    integer :: tmp
    if(ichan .lt. lbound(arr, 1) .OR. ichan .gt. ubound(arr, 1)) call die("Trying to swap ints,&
      & but one of the indices exceeds the bounds of the int array")
    if(jchan .lt. lbound(arr, 1) .OR. jchan .gt. ubound(arr, 1)) call die("Trying to swap ints,&
      & but one of the indices exceeds the bounds of the int array")
    tmp = arr(ichan)
    arr(ichan) = arr(jchan)
    arr(jchan) = tmp
  end subroutine swap_ints

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  impure module subroutine swap_channels(channels, ichan, jchan)
    !! Swaps the channels at channels(ichan) and channels(jchan)
    use rotex__system, only: die
    implicit none (type, external)
    class(channel_type), intent(inout) :: channels(:)
    integer, intent(in) :: ichan, jchan
    class(channel_type), allocatable :: tmp
    if(ichan .lt. lbound(channels, 1) .OR. ichan .gt. ubound(channels, 1)) call die("Trying to swap channels,&
      & but one of the indices exceeds the bounds of the channel array")
    if(jchan .lt. lbound(channels, 1) .OR. jchan .gt. ubound(channels, 1)) call die("Trying to swap channels,&
      & but one of the indices exceeds the bounds of the channel array")
    allocate(tmp, source = channels(ichan))
    tmp = channels(ichan)
    channels(ichan) = channels(jchan)
    channels(jchan) = tmp
    deallocate(tmp)
  end subroutine swap_channels

! ================================================================================================================================ !
end module rotex__channel_ops
! ================================================================================================================================ !
