! ================================================================================================================================ !
module rotex__symmetry
  !! All things related to symmetry that aren't point just groups

  use rotex__globals,     only: G
  use rotex__system,      only: die
  use rotex__pointgroups, only: pg_nrot

  implicit none (type, external)

  private

  public :: spin_symmetry
  public :: is_spin_allowed
  public :: is_spin_forbidden
  public :: symtop_rotstate_is_allowed
  public :: rotstate_is_allowed
  public :: sigma_v_class

  interface is_spin_allowed
    module procedure :: is_spin_allowed_chan
    module procedure :: is_spin_allowed_qnums
  end interface is_spin_allowed

  interface is_spin_forbidden
    module procedure :: is_spin_forbidden_chan
    module procedure :: is_spin_forbidden_qnums
  end interface is_spin_forbidden

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function spin_symmetry(n, ka, kc) result(res)
    !! Returns the nuclear spin symmetry class of the current N, Ka, Kc state.
    !! Two states with the same return value of this function can interconvert via electron collisions.
    !! Two states with different return values cannot. If G%ENFORCE_SPIN_ISOMER is .false., then
    !! this function always returns 0. Nuclear spin symmetry should still be recoverable numerically,
    !! if the scattering calculations were performed in highest possible Abelian point group, but this
    !! function can be used with G%ENFORCE_SPIN_ISOMER=.true. to force symmetry class separation if this
    !! is not enough.
    !!
    !!   Linear rotors: N mod 2 (centrosymmetry; linear molecules that are, e.g., D∞h )
    !!   Asymmetric rotors: K mod pg_nrot(target_pg)
    !!     K is Ka, Ka+Kc, or Kc depending on if G%SYMAXIS is A, B, or C
    !!   Symmetric rotors: K mod pg_nrot(target_pg)
    !!     K is Ka or Kc depending on if G%SYMAXIS is A or C (B disallowed)

    implicit none (type, external)

    integer, intent(in) :: n, ka, kc
    integer :: res
    integer :: ksym

    if(G%ENFORCE_SPIN_ISOMER .eqv. .false.) then
      res = 0
      return
    endif

    select case(G%ROTOR_KIND)
    case("l", "L")

      res = modulo(n,2)

    case("s","S","a","A")

      select case(G%SYMAXIS)
      case("a","A") ; ksym = ka
      case("b","B") ; ksym = ka+kc
      case("c","C") ; ksym = kc
      case default
        ksym = 0
      end select

      res = sigma_v_class(ksym, pg_nrot(G%TARGET_POINT_GROUP))

    case default

      res = 0

    end select

  end function spin_symmetry

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function is_spin_allowed_qnums(nlo, kalo, kclo, nup, kaup, kcup) result(res)
    !! Determine if the transition Nlo,Kalo,Kclo -> Nup,Kaup,Kcup is allowed by nuclear spin symmetry
    !! selection rules
    implicit none (type, external)
    integer,      intent(in) :: nlo, kalo, kclo, nup, kaup, kcup
    logical :: res
    res = spin_symmetry(nlo, kalo, kclo) .eq. spin_symmetry(nup, kaup, kcup)
  end function is_spin_allowed_qnums
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental function is_spin_allowed_chan(channel1, channel2) result(res)
    !! Test if two rotational channels respect ortho/para symmetry
    use rotex__types, only: asymtop_rot_channel_type
    implicit none (type, external)
    type(asymtop_rot_channel_type), intent(in) :: channel1, channel2
    logical :: res
    res = .false.
    if(channel1 % sym .ne. channel2 % sym) return
    res = .true.
  end function is_spin_allowed_chan

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function is_spin_forbidden_qnums(nlo, kalo, kclo, nup, kaup, kcup) result(res)
    !! Determine if the transition Nlo,Kalo,Kclo -> Nup,Kaup,Kcup is forbidden by nuclear spin symmetry
    !! selection rules
    implicit none (type, external)
    integer,      intent(in) :: nlo, kalo, kclo, nup, kaup, kcup
    logical :: res
    res = .not. is_spin_allowed_qnums(nlo, kalo, kclo, nup, kaup, kcup)
  end function is_spin_forbidden_qnums
  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental function is_spin_forbidden_chan(channel1, channel2) result(res)
    !! Test if two rotational channels respect ortho/para symmetry
    use rotex__types, only: asymtop_rot_channel_type
    implicit none (type, external)
    type(asymtop_rot_channel_type), intent(in) :: channel1, channel2
    logical :: res
    res = .true.
    if(channel1 % sym .ne. channel2 % sym) return
    res = .false.
  end function is_spin_forbidden_chan

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental function symtop_rotstate_is_allowed(N, K, rchar) result(res)
    !! Test whether the rotational state (N,K) is allowed, based on G%FORBIDDEN_STATES_KIND
    !!   kind = 0: all states allowed
    !!   kind = 1: all states allowed, except for even N, K=0 states. This is a legacy
    !!             check, predating wangify_symtop_eigvecs. It should probably be disabled.
    !!             The correct version is kind=2; this is basically half of kind=1
    !!   kind = 2: based on character. States are forbidden if it transforms as A₁ w.r.t. both
    !!             C_n (K = 0 mod n) and C₂' (rchar = +1).
    use rotex__functions,   only: isodd
    use rotex__pointgroups, only: pg_nrot
    implicit none (type, external)
    integer, intent(in) :: N, K, rchar
    logical :: res
    integer :: nrot
    select case(G%FORBIDDEN_STATES_KIND)
    case(0) ; res = .true.
    case(1) ; res = isodd(N) .OR. K .ne. 0
    case(2)
      nrot = pg_nrot(G%TARGET_POINT_GROUP)
      res  = .true.
      if(nrot .lt. 2) return
      if(modulo(K, nrot) .ne. 0) return ! <-- E-class: both characters exist; fine
      res = rchar .ne. +1               ! <-- A-class: drop A1 member
    case default
      call die("Unexpected G%FORBIDDEN_STATES_KIND. Must be one of 0,1,2")
    end select
  end function symtop_rotstate_is_allowed

  ! ------------------------------------------------------------------------------------------------------------------------------- !
  pure elemental function rotstate_is_allowed(N, Ka, Kc, rchar) result(res)
    !! Test whether the rotational state (N,K) is allowed, check special cases
    use rotex__functions, only: isodd
    implicit none (type, external)
    integer, intent(in) :: N, Ka, Kc, rchar
    logical :: res
    integer :: Ksym
    select case(G%FORBIDDEN_STATES_KIND)

    case(0)

      res = .true.

    case(1,2)

      select case(G%ROTOR_KIND)
      case("s")

        select case(G%SYMAXIS)
        case("a") ; Ksym = Ka
        case("b") ; call die("Ksym = B is ambiguous; pick A or C")
        case("c") ; Ksym = Kc
        end select

        res = symtop_rotstate_is_allowed(N, Ksym, rchar)

      case("a", "l")

        call die("G%FORBIDDEN_STATES_KIND ≠ 0 is not meaningful for a non-symmetric-top rotor")

      end select

    case default

      call die("Unexpected G%FORBIDDEN_STATES_KIND. Must be one of 0,1,2")

    end select
  end function rotstate_is_allowed

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure elemental function sigma_v_class(m, n) result(class)
    !! Returns the σ_v equivalence class of a projection quantum number m under an n-fold principal axis.
    !!   n<2: no axis, returns 0
    !!   n=2,3: {0,1,2}
    !!   n=4,5: {0,1,2,3}
    !!   etc.
    implicit none(type, external)
    integer, intent(in) :: m, n
    integer :: class
    class = 0
    if(n .lt. 2) return
    class = modulo(m, n)
    class = min(class, n - class)
  end function sigma_v_class

! ================================================================================================================================ !
end module rotex__symmetry
! ================================================================================================================================ !
