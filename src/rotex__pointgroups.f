! ================================================================================================================================ !
module rotex__pointgroups
  !! Point groups related stuff

  use rotex__system,     only: die, stderr
  use rotex__types,      only: pg_info_type
  use rotex__characters, only: lower

  implicit none (type, external)

  private

  public :: is_supported_pg
  public :: is_abelian_pg
  public :: group_size
  public :: pg_nrot
  public :: pg_nrot_needed
  public :: req_scat_pg
  public :: is_subgroup
  public :: abelian_pg_list
  public :: supported_pg_list
  public :: get_group_irreps
  public :: write_pg_table
  public :: irrep_name

  integer,      parameter, public :: PG_NROT_UNSUPPORTED = -1
  character(*), parameter, public :: PG_UNSUPPORTED = "NAPG" !! Unsupported point group

  integer, parameter, public :: even = 1
  integer, parameter, public :: odd  =-1

  type(pg_info_type), parameter, public :: PG_TABLE(*) = [                                                                   &
      !            name  nelem nrot nrot_needed req_scat is_abelian   irreps
      pg_info_type("c1 ",  1,   1,    1,        "c1 ",   .true.,  ["A  ", "   ", "   ", "   ", "   ", "   ", "   ", "   "]), &
      pg_info_type("cs ",  2,   1,    1,        "cs ",   .true.,  ["Ap ", "App", "   ", "   ", "   ", "   ", "   ", "   "]), &
      pg_info_type("ci ",  2,   1,    1,        "ci ",   .true.,  ["Ag ", "Au ", "   ", "   ", "   ", "   ", "   ", "   "]), &
      pg_info_type("c2 ",  2,   2,    1,        "c2 ",   .true.,  ["A  ", "B  ", "   ", "   ", "   ", "   ", "   ", "   "]), &
      pg_info_type("c2v",  4,   2,    1,        "c2v",   .true.,  ["A1 ", "B1 ", "B2 ", "A2 ", "   ", "   ", "   ", "   "]), &
      pg_info_type("c2h",  4,   2,    1,        "c2h",   .true.,  ["Ag ", "Au ", "Bg ", "Bu ", "   ", "   ", "   ", "   "]), &
      pg_info_type("d2 ",  4,   2,    1,        "d2 ",   .true.,  ["A  ", "B1 ", "B2 ", "B3 ", "   ", "   ", "   ", "   "]), &
      pg_info_type("d2h",  8,   2,    1,        "d2h",   .true.,  ["Ag ", "Au ", "B1g", "B1u", "B2g", "B2u", "B3g", "B3u"]), &
      pg_info_type("c3v",  6,   3,    3,        "cs ",   .false., ["A1 ", "A2 ", "E  ", "   ", "   ", "   ", "   ", "   "]), &
      pg_info_type("d3h", 12,   3,    3,        "c2v",   .false., ["A1p", "A2p", "Ep ", "A1q", "A2q", "Eq ", "   ", "   "]), &
      ! -- tbh we don't really need to non-abelian point group irreps bc we never look em up
      pg_info_type("d4h", 16,   4,    4,        "d2h",   .false., ["   ","   ","   ","   ","   ","   ","   ","   "]), &
      pg_info_type("d5h", 20,   5,    5,        "c2v",   .false., ["   ","   ","   ","   ","   ","   ","   ","   "]), &
      pg_info_type("d6h", 24,   6,    3,        "d2h",   .false., ["   ","   ","   ","   ","   ","   ","   ","   "])  &
  ]
  character(*), parameter :: PG_HEAD_FMT = '(A10, 1X, 3(A10, 1X), A11, 1X, A7, 3X, A )'
  character(*), parameter :: PG_ITEM_FMT = '(A10, 1X, 3(I10, 1X), A11, 1X, L7, 3X, 8(A3, 1X) )'

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function find_pg(pg) result(idx)
    !! Returns index of point group in the table, -1 if not found
    implicit none(type, external)
    character(*), intent(in) :: pg
    integer :: idx
    character(:), allocatable :: pg_lo
    pg_lo = trim(lower(pg))
    do idx = 1, size(PG_TABLE)
      if(pg_lo .ne. trim(PG_TABLE(idx)%name)) cycle
      return
    enddo
    idx=-1
  end function find_pg

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_supported_pg(pg) result(res)
    !! Is PG in the table ?
    implicit none(type, external)
    character(*), intent(in) :: pg
    logical :: res
    res = find_pg(pg) .gt. 0
  end function

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_abelian_pg(pg) result(res)
    !! Is PG an Abelian point group ?
    implicit none(type, external)
    character(*),  intent(in) :: pg
    logical :: res
    integer :: idx
    idx = find_pg(pg)
    res = idx .gt. 0
    if(res .eqv. .false.) return
    res = PG_TABLE(idx) % is_abelian
  end function is_abelian_pg

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function group_size(pg) result(n)
    !! Returns number of elements in the group PG, or -1 if PG accepted
    implicit none(type, external)
    character(*),  intent(in) :: pg
    integer :: n
    integer  :: idx
    idx = find_pg(pg)
    if(idx .lt. 1) then
      n = -1
      return
    endif
    n = PG_TABLE(idx)%nelem
  end function group_size

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function pg_nrot(pg) result(nrot)
    !! Returns the order of PG's principal cyclic rotation axis (1 for groups without rotation),
    !! or the  value of PG_NROT_UNSUPPORTED for unknown PGs
    implicit none(type, external)
    character(*), intent(in) :: pg
    integer :: nrot
    integer :: idx
    idx = find_pg(pg)
    if(idx .lt. 1) then
      nrot = PG_NROT_UNSUPPORTED
      return
    endif
    nrot = PG_TABLE(idx)%nrot
  end function pg_nrot

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function pg_nrot_needed(pg) result(nrot)
    !! Return the order of the cyclic rotation recovered by K-matrix projection from PG's max
    !! Abelian subgroup. 1 if `pg` is Abelian. PG_NROT_UNSUPPORTED for unknown PGs.
    implicit none(type, external)
    character(*), intent(in) :: pg
    integer :: nrot
    integer :: idx
    idx = find_pg(pg)
    if(idx .lt. 1) then
      nrot = PG_NROT_UNSUPPORTED
      return
    endif
    nrot = PG_TABLE(idx)%nrot_needed
  end function pg_nrot_needed

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function req_scat_pg(pg) result(res)
    !! Max Abelian subgroup in which the K-matrix must be computed to recover PG symmetry.
    !! Returns PG_UNSUPPORTED for unknown PGs
    implicit none(type, external)
    character(*), intent(in) :: pg
    character(:), allocatable :: res
    integer :: idx
    idx = find_pg(pg)
    if(idx .lt. 1) then
      res = PG_UNSUPPORTED
      return
    endif
    res = trim(PG_TABLE(idx)%req_scat_pg)
  end function req_scat_pg

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function is_subgroup(pg1, pg2) result(res)
    !! Test if PG1 ⊆ PG2 within the supported PG set. Returns false if either
    !! PG is not in PG_TABLE, or if no subgroup relation is encoded for the pair
    implicit none(type, external)
    character(*), intent(in) :: pg1, pg2
    logical :: res
    character(:), allocatable :: p1, p2

    p1 = trim(lower(pg1))
    p2 = trim(lower(pg2))

    if(p1 .eq. p2) then
      res = .true.
      return
    endif

    ! -- the people's subgroup (subrgroup of every other group)
    if(p1 .eq. "c1") then
      res = is_supported_pg(p2)
      return
    endif

    select case(p2)
    case("c2v") ; res = any(p1 .eq. ["cs", "c2"])
    case("c2h") ; res = any(p1 .eq. ["cs", "ci", "c2"])
    case("d2")  ; res =     p1 .eq. "c2"
    case("d2h") ; res = any(p1 .eq. ["cs ", "ci ", "c2 ", "c2v", "c2h", "d2 "])
    case("c3v") ; res =     p1 .eq. "cs"
    case("d3h") ; res = any(p1 .eq. ["cs ", "c2 ", "c2v"])
    case("d4h") ; res = any(p1 .eq. ["cs ", "ci ", "c2 ", "c2v", "c2h", "d2 ", "d2h"])
    case("d5h") ; res = any(p1 .eq. ["cs ", "c2 ", "c2v"])
    case("d6h") ; res = any(p1 .eq. ["cs ", "ci ", "c2 ", "c2v", "c2h", "d2 ", "d2h"])
    case default ; res = .false.
    end select
  end function is_subgroup

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function abelian_pg_list() result(list)
    !! Returns an array of names of all Abelian PGs in PG_TABLE, in table order.
    implicit none(type, external)
    character(:), allocatable :: list(:)
    integer :: i, j, nabelian
    nabelian = count(PG_TABLE(:)%is_abelian)
    allocate(character(3) :: list(nabelian))
    j = 0
    do i = 1, size(PG_TABLE, 1)
      if(PG_TABLE(i)%is_abelian .eqv. .false.) cycle
      j = j + 1
      list(j) = PG_TABLE(i)%name
    enddo
  end function abelian_pg_list

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function supported_pg_list() result(list)
    !! Returns an array of names of all supported PGs in PG_TABLE, in table order.
    character(:), allocatable :: list(:)
    integer :: i
    allocate(character(3) :: list(size(PG_TABLE, 1)))
    do i = 1, size(PG_TABLE)
      list(i) = PG_TABLE(i)%name
    enddo
  end function supported_pg_list

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine get_group_irreps(pg, irreps)
    implicit none(type, external)
    character(*),              intent(in)  :: pg
    character(:), allocatable, intent(out) :: irreps(:)
    integer :: idx, n_irreps
    idx = find_pg(pg)
    if(idx .lt. 1) then
      call write_pg_table(stderr)
      call die("Unknown point group: " // pg)
    endif
    n_irreps = PG_TABLE(idx) % nelem
    allocate(character(3) :: irreps(n_irreps))
    irreps(:) = PG_TABLE(idx)%irreps(1:n_irreps)
  end subroutine

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  function irrep_name(irrep, pg) result(output)
    integer,      intent(in) :: irrep
    character(*), intent(in) :: pg
    character(:), allocatable :: output
    integer :: idx
    idx = find_pg(pg)
    if(idx .lt. 1) then
      call write_pg_table(stderr)
      call die("Unknown point group: " // pg)
    endif
    output = trim(PG_TABLE(idx)%irreps(irrep))
  end function

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine write_pg_table(funit)
    !! Write the global point group table to funit
    implicit none(type, external)
    integer, intent(in) :: funit
    integer :: i
    write(funit, *)
    write(funit, *)
    write(funit, PG_HEAD_FMT) "NAME", "NELEM", "NROT", "NROT_NEED", "REQ_SCAT_PG", "ABELIAN", "IRREPS.."
    write(funit, *)
    do i=1, size(PG_TABLE, 1)
      write(funit, PG_ITEM_FMT)   &
        PG_TABLE(i) % name        &
      , PG_TABLE(i) % nelem       &
      , PG_TABLE(i) % nrot        &
      , PG_TABLE(i) % nrot_needed &
      , PG_TABLE(i) % req_scat_pg &
      , PG_TABLE(i) % is_abelian  &
      , PG_TABLE(i) % irreps(:)
    enddo
    write(funit, *)
  end subroutine write_pg_table

! ================================================================================================================================ !
end module rotex__pointgroups
! ================================================================================================================================ !
