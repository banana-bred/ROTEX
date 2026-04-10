! ================================================================================================================================ !
module rotex__CBXS
  !! Routines to calculate cross sections in the Coulomb-Born approximation

  use rotex__globals, only: G

  implicit none (type, external)

  private

  public :: get_einsta_only
  public :: get_CB_xs_asym
  public :: xtrapolate_cb_xs

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_einsta_only( &
        einsta                       &
      , nlo                          &
      , nup                          &
      , elo                          &
      , eup                          &
      , eigveclo                     &
      , eigvecup                     &
      , dipole_moments               &
      , quadrupole_moments           &
      )
    !! Calculate only the Einstein A coefficeints for a transition

    use rotex__kinds,      only: dp
    use rotex__system,     only: die
    use rotex__constants,  only: invc
    use rotex__functions,  only: istriangle

    implicit none (type, external)

    integer, intent(in) :: nlo
      !! the angular momentum quantum number \(N\)
    integer, intent(in) :: nup
      !! the angular momentum quantum number \(N'\)
    real(dp), intent(in) :: elo
      !! The energy (hartrees) of the initial state
    real(dp), intent(in) :: eup
      !! The energy (hartrees) of the final state
    complex(dp), intent(in), allocatable :: eigveclo(:)
      !! Eigenvector of the initial state in the basis of symmetric top wavefunctions; the coefficients
      !! \(c^{N,\tau}_K)
    complex(dp), intent(in), allocatable :: eigvecup(:)
      !! Eigenvector of the final state in the basis of symmetric top wavefunctions; the coefficients
      !! \(c^{N',\tau'}_{K})
    real(dp), intent(inout) :: einsta
      !! The Einstein coefficient for the transition
    complex(dp), intent(in), target :: dipole_moments(3)
      !! The spherical dipole moments
    complex(dp), intent(in), target :: quadrupole_moments(5)
      !! The spherical quadrupole moments

    integer :: lambda
    real(dp) :: omega, am_summation
    complex(dp), pointer  :: multipole_moments(:)

    omega = eup - elo

    multipoles: do lambda = 1, 1
      am_summation = 0
      select case(lambda)
      case(1)
        if(G%DO_DIPOLE .eqv. .false.) cycle multipoles
        multipole_moments => dipole_moments
      case(2)
        if(G%DO_QUADRUPOLE .eqv. .false.) cycle multipoles
        multipole_moments => quadrupole_moments
      case default
        call die("Somehow, λ ≠ 1 or 2")
      end select

      ! -- enforce 3j rules
      if( .false. .eqv. istriangle(nlo, nup, lambda) ) cycle multipoles

      ! -- do the summation over angular momenta and multipole components
      !    or use the CDMS Einstein A coeffs for the dipole case
      select case(lambda)
      case(1)
        if((G%USE_CDMS_EINSTA .eqv. .false.) .OR. einsta .eq. 0) then
          ! -- calculate the am_summation ourselves if not using the CDMS
          !    or if it was not present in the CDMS data
          call multipole_am_summation_asym(am_summation, lambda, nlo, nup, multipole_moments, eigveclo, eigvecup)
        else
          ! -- use pre-existing Einstein A coeff from the CDMS read
          nullify(multipole_moments)
          cycle multipoles
        endif
      case(2)
        call multipole_am_summation_asym(am_summation, lambda, nlo, nup, multipole_moments, eigveclo, eigvecup)
      end select

      nullify(multipole_moments)

      select case(lambda)
      case(1)
        einsta = einsta &
               + am_summation * 4._dp/3._dp * (omega * invc)**3 * (2*nlo+1)
               ! + am_summation * 4._dp * (omega * invc)**3 * (2*nlo+1)
      case(2)
        ! einsta = einsta + am_summation * 2._dp/15._dp * (omega * invc)**5 * (2*nlo+1)
      end select

    enddo multipoles

  end subroutine get_einsta_only

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_CB_xs_asym( energies           &
                                  , sigma              &
                                  , N                  &
                                  , Np                 &
                                  , E                  &
                                  , Ep                 &
                                  , eigvec             &
                                  , eigvecp            &
                                  , einsta             &
                                  , spherical_dipole_moments     &
                                  , spherical_quadrupole_moments &
                                  , analytic_total_cb  &
                                  , lmax               &
                                  )
                               ! , analytic_total_cb, lmax, atol, rtol)
    !! Calculate the excitation and de-excitation cross sections (xs) for an asymmetric top up.
    !! The sum over partial waves is either truncated to \(l_\text{max}\) or determined analytically.
    !! The summation over the angular momentum components, multipole terms, and partial waves are separable
    !! for each value of λ as summation = (sum over angular momentum and multipole moments) \(\times\) (sum over partial waves).
    !! The first sum does not involve the scattering electron, which means that it is independent
    !! of the target's charge. The second sum is the only one that changes whether the target
    !! is neutral.

    use rotex__kinds,      only: dp
    use rotex__system,     only: die
    use rotex__constants,  only: invc, pi
    use rotex__functions,  only: istriangle
    use rotex__wigner,     only: wigner3j

    implicit none (type, external)

    real(dp), intent(in) :: energies(:)
      !! Array of scattering energies to consider
    real(dp), intent(out), allocatable :: sigma(:)
      !! Cross sections calculated on a grid of scattering energies for \(Nτ \rightarrow N'τ'\)
      !! Returned with the same size as `energies`
    integer, intent(in) :: N
      !! the angular momentum quantum number \(N\)
    integer, intent(in) :: Np
      !! the angular momentum quantum number \(N'\)
    real(dp), intent(in) :: E
      !! The energy (hartrees) of the initial state
    real(dp), intent(in) :: Ep
      !! The energy (hartrees) of the final state
    complex(dp), intent(in) :: eigvec(:)
      !! Eigenvector of the initial state in the basis of symmetric top wavefunctions; the coefficients
      !! \(c^{N,\tau}_K)
    complex(dp), intent(in) :: eigvecp(:)
      !! Eigenvector of the final state in the basis of symmetric top wavefunctions; the coefficients
      !! \(c^{N',\tau'}_{K})
    real(dp), intent(inout) :: einsta
      !! The Einstein coefficient for the transition
    complex(dp), intent(in), target :: spherical_dipole_moments(3)
      !! The spherical dipole moments
    complex(dp), intent(in), target :: spherical_quadrupole_moments(5)
      !! The spherical quadrupole moments
    logical, intent(in) :: analytic_total_cb(:)
      !! Array of values telling us whether we want to use the analytic expression for lmax -> infintiy
    integer, intent(in) :: lmax
      !! The max value of the orbital angular momentum quantum number l to consider

    complex(dp), pointer :: multipole_moments(:)

    integer  :: nE
    integer  :: lambda
    real(dp) :: am_summation
    real(dp) :: omega
    real(dp), allocatable :: pw_summation(:)

    nE    = size(energies, 1)
    omega = Ep - E

    ! -- (re)allocate sigma arrays if necessary
    if(allocated(sigma)) then
      if(size(sigma, 1) .ne. nE) then
        deallocate(sigma)
        allocate(sigma(nE))
      endif
    else
      allocate(sigma(nE))
    endif

    allocate(pw_summation(nE))

    sigma   = 0

    ! -- loop over dipole, quadrupole terms
    multipoles: do lambda = 1, 1
    ! multipoles: do lambda = 1, 2

      am_summation = 0
      pw_summation = 0

      select case(lambda)
      case(1)
        if(G%DO_DIPOLE     .eqv. .false.) cycle multipoles
        multipole_moments => spherical_dipole_moments
      case(2)
        if(G%DO_QUADRUPOLE .eqv. .false.) cycle multipoles
        multipole_moments => spherical_quadrupole_moments
      case default
        call die("Somehow, lambda is neither 1 nor 2")
      end select

      ! -- enforce 3j rules to avoid extra computation
      if( .false. .eqv. istriangle(N, Np, lambda) ) cycle multipoles

      ! -- do the summation over angular momenta and multipole components
      !    or use the CDMS Einstein A coeffs for the dipole case
      select case(lambda)
      case(1)
        if((G%USE_CDMS_EINSTA .eqv. .false.) .OR. einsta .eq. 0) then
          ! -- calculate the am_summation ourselves
          select case(G%ROTOR_KIND)
          case("l")
            ! -- |μz|² * (N,N,λ;0,0,0)²
            am_summation = abs(multipole_moments(2))**2 * wigner3j(N,Np,lambda,0,0,0)**2
          case("s")
            ! -- the reduced angular momentum and multipole moment summation
            !    for symmetric tops
            call die("multipole am summation not defined yet for symmetric tops")
          case("a")
            ! -- the full angular momentum and multipole moment summation
            !    for asymmetric tops
            call multipole_am_summation_asym(am_summation, lambda, N, Np, multipole_moments, eigvec, eigvecp)
          end select
        else
          ! -- use pre-existing Einstein A coeff from the CDMS read
          am_summation = einsta * 3._dp/4._dp * (omega * invc)**(-3) / (2*N+1)
        endif
      case(2)
        call die("quadrupoles not implemented yet")
        ! call multipole_am_summation_asym(am_summation, lambda, N, Np, multipole_moments, eigvec, eigvecp)
      end select

      nullify(multipole_moments)

      if(abs(am_summation) .lt. 1e-16) cycle multipoles

      if(am_summation .lt. 0) call die("The angular momentum summation is negative, which is not physical !")

      ! -- if the summation is basically 0 for this multipole term, check the next one
      ! if(am_summation .lt. 1e-16_dp) cycle multipoles

      ! -- do the summation over partial waves; independent of rotor kind
      if(analytic_total_cb(lambda) .eqv. .true.) then
        call get_infinite_pwsum(pw_summation, energies, E, Ep, lambda)
      else
        call get_truncated_pwsum(pw_summation, energies, E, Ep, lambda, lmax)
      endif

      if(any(pw_summation .lt. 0)) &
        call die("The partial wave summation is negative for at least one energy, which is not physical !")

      !   σ =         4π/(2λ+1)              * (target)     * (electron)
      sigma = sigma + 4.0_dp*pi/(2*lambda+1) * am_summation * pw_summation

      select case(lambda)
      case(1)
        if((G%USE_CDMS_EINSTA .eqv. .false.) .OR. einsta .eq. 0) then
          einsta = einsta &
                 + am_summation * 4._dp/3._dp * (omega * invc)**3 * (2*N+1)
        endif
      case(2)
        ! einsta = einsta + am_summation * 2._dp/15._dp * (omega * invc)**5 * (2*N+1)
      end select

    enddo multipoles

    nancheck: block
      use rotex__system, only: stderr
      use ieee_arithmetic, only: ieee_is_nan
      integer :: i
      integer, allocatable :: indices(:)
      if(all(ieee_is_nan(sigma) .eqv. .false.)) exit nancheck
      indices = [(i,i=1,size(sigma,1))]
      indices = pack(indices, ieee_is_nan(sigma) .eqv. .true.)
      write(stderr, '(A)', advance = "no") "NaNs detected at the following indices:"
      write(stderr, *) indices
      write(stderr, '(A, I0)') "Lmax = ", lmax
      call die("NaN(s) detected in cross sections ! Maybe lmax is too big for wigner&
      & or Ei is too small for the hypergeometric functions (or their normalizations)")
    end block nancheck

    ! -- σ: N Ka Kc -> N' Ka' Kc' (2N'+1 factor for excitation cross sections)
    sigma = sigma * (2*Np + 1)

    ! -- A: N Ka Kc <- N' Ka' Kc' (divide out the initial 2N'+1)
    ! einsta  = einsta * (2*N+1)

  end subroutine get_CB_xs_asym

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine multipole_am_summation_asym(summation, lambda, N, Np, multipole_moments, eigvec, eigvecp)
    !! Carry out the summation of the angular momentum projections
    !! and multipole components
    use rotex__kinds,  only: dp
    use rotex__system, only: die
    use rotex__wigner, only: wigner3j
    use rotex__functions, only: neg

    implicit none (type, external)

    real(dp),    intent(out) :: summation
      !! The result of the summation
    integer,     intent(in)  :: lambda
      !! The multipole expansion term
      !! 1: dipole
      !! 2: quadrupole
    integer,     intent(in)  :: N
      !! the quantum number N
    integer,     intent(in)  :: Np
      !! the quantum number N'
    complex(dp), intent(in)  :: multipole_moments(:)
      !! array of spherical multipole moments
    complex(dp),    intent(in)  :: eigvec(:)
      !! The eigenvector for the initial state Nτ
    complex(dp),    intent(in)  :: eigvecp(:)
      !! The eigenvector for the final state N'τ'

    integer :: ilambda1, ilambda2
    integer :: mu1, mu2
    integer :: K1, K2
    integer :: K1p, K2p
    integer :: i_K1, i_K2
    integer :: i_K1p, i_K2p

    complex(dp) :: tmpz
    complex(dp) :: zummation

    zummation = 0

    ! -- sum over μ, K, etc
    do mu1 = -lambda, lambda
      do K1 = -N, N
        ! do K1p = -Np, Np

        ! -- the only nonzero 3j term, cycle if it's not physical
        K1p = K1 + mu1
        if(abs(K1p) .gt. Np) cycle

        ilambda1 = mu1 + lambda + 1
        i_K1  = (K1  + N ) + 1
        i_K1p = (K1p + Np) + 1

        do mu2 = -lambda, lambda
          do K2 = -N, N
            ! do K2p = -Np, Np

            ! -- the only nonzero 3j term, cycle if it's not physical
            K2p = K2 + mu2
            if(abs(K2p) .gt. Np) cycle

            ilambda2 = mu2 + lambda + 1
            i_K2  = (K2  + N ) + 1
            i_K2p = (K2p + Np) + 1

            tmpz =       multipole_moments(ilambda1)  &
                 * conjg(multipole_moments(ilambda2)) &
                 ! * eigvecp(i_K1p)        * conjg(eigvec(i_K1))  &
                 * neg(K1+K2) &
                 * eigvec(i_K1)   * conjg(eigvecp(i_K1p)) &
                 * eigvecp(i_K2p) * conjg(eigvec(i_K2))       &
                 * wigner3j(2*N, 2*Np, 2*lambda, -2*K1,  2*K1p, -2*mu1)   &
                 * wigner3j(2*N, 2*Np, 2*lambda, -2*K2,  2*K2p, -2*mu2)
            ! tmpz =       multipole_moments(ilambda1)  &
            !      * conjg(multipole_moments(ilambda2)) &
            !      * eigvecp(i_K1p)        * conjg(eigvec(i_K1))  &
            !      * conjg(eigvecp(i_K2p)) * eigvec(i_K2)       &
            !      * wigner3j(2*N, 2*Np, 2*lambda, -2*K1,  2*K1p, -2*mu1)   &
            !      * wigner3j(2*N, 2*Np, 2*lambda, -2*K2,  2*K2p, -2*mu2)
                 ! * wigner3j(2*N, 2*Np, 2*lambda,  2*K1, -2*K1p,  2*mu1)   &
                 ! * wigner3j(2*N, 2*Np, 2*lambda,  2*K2, -2*K2p,  2*mu2)
            zummation = zummation + tmpz

          ! enddo ; enddo
          enddo
        enddo

      enddo
      ! enddo ; enddo
    enddo

    if(abs(zummation%im) .gt. 1e-16) call die("Large imaginary cross section detected. Check your multipole moments.")

    summation = zummation % re

  end subroutine multipole_am_summation_asym

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine get_infinite_pwsum(summation, energies, E, Ep, lambda)
    !! Return the analytic expression for the partial wave sum in the limit \(l\to\infty\)
    !! for the Coulomb or neutral case
    use rotex__kinds,          only: dp
    use rotex__system,         only: die, stderr
    use ieee_arithmetic,       only: ieee_is_nan
    use rotex__constants,      only: im
    use rotex__functions,      only: expm1
    use rotex__hypergeometric, only: f21

    implicit none (type, external)

    real(dp), intent(out) :: summation(:)
      !! the result of the summation, for each energy
    real(dp), intent(in) :: energies(:)
      !! the scattering eneriges
    real(dp), intent(in) :: E
      !! the energy of the initial state Nτ
    real(dp), intent(in) :: Ep
      !! the energy of the final state N'τ'
    integer, intent(in) :: lambda
      !! the multipole term
      !! 1: dipole
      !! 2: quadrupole

    logical :: flag
    integer :: ie, nE
    integer :: iemin
    real(dp) :: k, kp
    real(dp) :: dE

    interface
      function func(k, kp) result(res)
        import dp
        real(dp), intent(in) :: k, kp
        real(dp) :: res
      end function func
    end interface

    procedure(func), pointer :: pw_limit => null()

    if(lambda .ne. 1) call die("analytic partial wave sum not available for lambda =/= 1")

    if(G%TARGCHARGE .eq. 0) then
      pw_limit => born_inf_dipole
    else
      pw_limit => coul_inf_dipole
    endif

    nE = size(energies, 1)
    dE = Ep - E

    flag = .false.
    summation = 0
    iemin = findloc(energies .gt. dE, .true., 1)

    !$omp parallel do default(none) &
    !$omp& shared(energies, summation, E, Ep, dE, lambda, G, iemin, nE, flag, pw_limit) &
    !$omp& private(ie, k, kp)
    do ie = iemin, nE
    ! do concurrent(ie = 1 : nE)

      k    = sqrt(2*energies(ie))
      kp   = sqrt(2*(energies(ie) - dE))

      summation(ie) = pw_limit(k, kp)

      if(summation(ie) .ge. 0) cycle
      if(flag .eqv. .true.) cycle
      flag = .true.

      block
        use rotex__constants, only: au2ev
        write(stderr, *)
        write(stderr, *) E*au2ev, Ep*au2ev
        write(stderr, *) energies(ie)*au2ev, (energies(ie) - dE)*au2ev
        write(stderr, *) ie, energies(ie), -G%TARGCHARGE/k, -G%TARGCHARGE/kp
        write(stderr, '(A)') "WARN: negative cross sections when evaluating the infinite partial wave sum."
      end block

    enddo
    !$omp end parallel do

    nullify(pw_limit)

  end subroutine get_infinite_pwsum

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  function coul_inf_dipole(k, kp) result(res)
    use rotex__kinds, only: dp
    use rotex__functions,      only: expm1
    use rotex__constants,      only: pi, im
    use rotex__hypergeometric, only: f21
    implicit none (type, external)
    real(dp), intent(in) :: k, kp
    real(dp) :: res
    real(dp), parameter :: logpi = log(pi)
    complex(dp), parameter :: onez = cmplx(1, kind=dp)
    real(dp) :: x, eta, etap, twopi_eta, twopi_etap, ea, eap, sgn
    real(dp) :: log_prefactor
    eta  = real(-G%TARGCHARGE, kind=dp)/k
    etap = real(-G%TARGCHARGE, kind=dp)/kp
    ! -- seems more stable than with η,η'
    x = -4._dp * k*kp / ( (k-kp)*(k-kp) )
    ! -- for the exponential terms
    twopi_eta  = 2*pi*eta
    twopi_etap = 2*pi*etap
    ea = expm1(twopi_eta)
    eap = expm1(twopi_etap)
    ! -- keep track of this sign because we do log(abs(ea))
    sgn = sign(1._dp, ea)*sign(1._dp, eap)
    ! -- prefactor in log space (combine negatives, take log of |x|)
    log_prefactor = log(4._dp) + 2*logpi - 3*log(k) - log(kp) + log(abs(x)) &
                  + twopi_eta - log(abs(ea)) - log(abs(eap))
    res = sgn * exp(log_prefactor) * real(                             &
            f21(im*eta, im*etap, onez, x)                        &
          * f21(-im*eta + 1, -im*etap + 1, 2*onez, x), kind = dp &
          )
  end function coul_inf_dipole

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  pure function born_inf_dipole(k, kp) result(res)
    !! The analytic integrated born cross section for the dipole term
    use rotex__kinds,     only: dp
    use rotex__constants, only: pi
    use rotex__system,    only: die
    implicit none (type, external)
    real(dp), intent(in) :: k, kp
    real(dp) :: res
    real(dp), parameter :: twopi = 2._dp * pi
    if(G%TARGCHARGE .ne. 0) call die("Z≠0 detected in born_inf_dipole")
    res = twopi / (k*kp) * (log(k+kp) - log(abs(k-kp)))
  end function born_inf_dipole

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine get_truncated_pwsum(summation, energies, E, Ep, lambda, lmax)
    !! Calculate the truncated partial wave sum
    !! \(\sum\limits_{ll'}^{\l_\text{max}} (2l+1)(2l'+1) (l,l,\lambda;0,0,0)^2 |M^\lambda_{ll'}|^2\), given in
    !!   "Electromagnetic Excitation: Theory of Coulomb Excitation with Heavy Ions " by Kurt Alder and Aage Winther, Chapter IX,
    !!    section 2, page 244, equation 14.
    !! The integral given in that text is \(I^\lambda_{ll'}=-4M^\lambda_{ll'}\) in atomic units, based on equations 7 and 11 of the
    !! same section. For values of λ larger than 1, a recursion formula is used for the integral \(M^\lambda_{ll'}\), given in
    !!   "Study of Nuclear Structure by Electromagnetic Excitation with Accelerated Ions" by K. Alder, A. Bohr, T. Huus,
    !!    B. Mottelson, and A. Winther, equation II B.70–71 on page 453 for the λ=1 case.
    use rotex__kinds,     only: dp
    use rotex__system,    only: stderr
    ! use WignerSymbol, only: wigner3j
    use rotex__wigner,    only: wigner3j
    use rotex__globals, only: CB_MINT_IMAG_THRESH
    use rotex__functions, only: istriangle
    ! use rotex__wigner, only: wigner3j => threej

    implicit none (type, external)

    real(dp), intent(out) :: summation(:)
      !! the summation to carry out for each energy
    real(dp), intent(in)  :: energies(:)
      !! scattering energies (au)
    real(dp), intent(in)  :: E
      !! energy of the initial (lower) state
    real(dp), intent(in)  :: Ep
      !! energy of the final (higher) state
    integer,  intent(in)  :: lambda
      !! the multipole>
      !! 1: dipole
      !! 2: quadrupole
    integer,  intent(in)  :: lmax
      !! the max partial wave to consider

    logical :: flag
    integer :: l
    integer :: iemin
    integer :: ie, nE
    real(dp) :: k, kp
    real(dp) :: dE

    ! real(dp), allocatable :: weights(:)
    ! real(dp), allocatable :: M1(:)
    ! real(dp), allocatable :: M2(:)
    real(dp) :: weights(0:lmax-lambda)
    complex(dp) :: M1(0:lmax-lambda)
    complex(dp) :: M2(0:lmax-lambda)

    abstract interface
      function Mint(lambda, ltarg, ki, kf, swapl) result(res)
        import dp
        implicit none (type, external)
        integer,  intent(in) :: lambda, ltarg
        real(dp), intent(in) :: ki, kf
        logical,  intent(in), optional :: swapl
        complex(dp) :: res(0:ltarg)
      end function Mint
    end interface

    procedure(Mint), pointer :: Mints => null()

    ! -- point to the respective integral to calculate
    if(G%TARGCHARGE .eq. 0) then
      Mints => Mborn_array
    else
      Mints => Mcoul_array
    endif

    nE = size(energies, 1)

    dE = Ep - E
    flag = .false.
    summation = 0

    ! -- the starting energy index to ensure positive energies
    iemin = findloc(energies .gt. dE, .true., 1)

    ! -- the sum over ll'
    !$omp parallel do default(none) &
    !$omp& shared(energies, summation, E, Ep, dE, lambda, G, iemin, nE, flag, lmax, Mints) &
    !$omp& private(ie, k, kp, weights, l, M1, M2)
    nrg: do ie = iemin, nE

      k  = sqrt(2*energies(ie))
      kp = sqrt(2*(energies(ie) - dE))

      ! -- build the 3j symbols (they're symmetric w.r.t l and lp for λ=1) and integrals that we'll sum over. We only
      !    go to l=lmax-λ because we will get the values M_c{l+λ, l} and the l+λ value is what will be at most lmax.
      !    Because the dipole recursion is second order and only mixed sequential values of l with the same λ, we can
      !    calculate each sequence separately,
      !      - M_{l+λ,l}(k,k')
      !      - M_{l+λ,l}(k',k)
      !    and then just take the sum over these arrays with their 3j weights
      weights = [( (2*l+1) * (2*(l+lambda)+1) * wigner3j(2*l, 2*(lambda+l), 2*lambda, 0, 0, 0)**2, l=0, lmax-lambda )]
      M1      = Mints(lambda, lmax-lambda, k, kp, swapl = .false.) ! l'=l+λ
      M2      = Mints(lambda, lmax-lambda, k, kp, swapl = .true.)  ! l'=l-λ

      if(any(abs(M1%im) .gt. CB_MINT_IMAG_THRESH) .OR. any(abs(M2%im) .gt. CB_MINT_IMAG_THRESH)) then
        write(stderr, *)
        write(stderr, '("WARN: Skipping energy ", I0, ": the Coulomb integrals are nonreal")') ie
        write(stderr, *)
        cycle nrg
      endif

      summation(ie) = 4.0_dp*kp/k * sum( weights * (abs(M1)**2 + abs(M2)**2) )
      ! summation(ie) = 16.0_dp*pi*kp/k * sum( weights * (abs(M1)**2 + abs(M2)**2) )

      if(summation(ie) .ge. 0)  cycle
      if(flag .eqv. .true.) cycle
      flag = .true.

      write(stderr, *)
      write(stderr, '(A)') "WARN: negative cross sections when evaluating the truncated partial wave sum."

    enddo nrg
   !$omp end parallel do

   nullify(Mints)

  end subroutine get_truncated_pwsum

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function Mcoul(l, ki, kf) result(res)
    !! Calculates the integral \(M^\lambda_{l+λ,l}\) via the expression given in
    !!   "Electromagnetic Excitation: Theory of Coulomb Excitation with Heavy Ions " by Kurt Alder and Aage Winther, Chapter IX,
    !!    section 2, page 244, equation 14.
    !! for λ = 1, where \(-4M^\lambda_{l+λ,l} = I^λ_{l+λ, l}\). There is an expression for λ = 2 in
    !!   "Study of Nuclear Structure by Electromagnetic Excitation with Accelerated Ions" by K. Alder, A. Bohr, T. Huus,
    !!   B. Mottelson, and A. Winther, but only the dipole is used (at least for now).
    !! This is the integral of Colomb wavefunctions that arises in Coulomb-scattering off of
    !! a target whose long-rage potential is expended into multipoles (although this is only the λ=1 case for now)

    use rotex__kinds,          only: dp
    use rotex__constants,      only: im, pi
    use rotex__polygamma,      only: lgamma => log_gamma
    use rotex__functions,      only: factorial
    use rotex__hypergeometric, only: f21

    implicit none (type, external)

    integer, intent(in) :: l
    real(dp), intent(in) :: ki, kf
    complex(dp) :: res

    real(dp) :: etai, etaf, deta, x0
    complex(dp) :: c1, c2

    etai = -G%TARGCHARGE/ki
    etaf = -G%TARGCHARGE/kf
    deta = etaf - etai
    x0   = -4*etaf*etai/deta**2

    c1 = cmplx(2*l+2, kind = dp)
    c2 = cmplx(2*l+4, kind = dp)

    ! -- the original formula assumes k > k' (ki > kf). This is the symmetrized version so that we can
    !    freely call this procedure while switching ki and kf without worry. The only changes are to the
    !    `etaf - etai` and `deta` terms. Of course, ki and kf are both assumed to be positive given the context
    !    of electronic excitation
    res = 2 * abs((etaf-etai)/(etaf+etai))**(im*(etaf+etai))                                               &
        * exp(pi*abs(deta)/2) * (-x0)**(l+1)/factorial(2*l+2)                                              &
        * abs( exp( lgamma(l+1+im*etai) + lgamma(l+2+im*etaf) ) )                                          &
        * ( etai * F21(l+1-im*etai, l+1-im*etaf, c1, x0)                                                   &
          - etaf * abs(l+1+im*etai)**2 / ((2*l+2)*(2*l+3)) * (-x0) * F21(l+2-im*etai, l+2-im*etaf, c2, x0) &
        )

    ! -- get the integral M, not I
    res = res/(-4.0_dp)

  end function Mcoul

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  function Mborn_array(lambda, ltarg, ki, kf, swapl) result(res)
    !! Calculate the integral \(M^\eta_{ll_0}(k_vk_0)\) (13c) in Feldt+Morrison2008
    !!   doi: 10.1103/PhysRevA.77.012726
    !! in the first Born approximation for scattering of a target that whose potential is
    !! expanded into multipoles. This is the integral of spherical Bessel functions resulting
    !! from the asymptotic wavefunctions in the case of a sub-Coulomb interaction potential at
    !! large distances.

    use rotex__kinds,     only: dp
    use rotex__system,    only: die
    use rotex__constants, only: pi
    use rotex__polygamma, only: lgamma => log_gamma
    use rotex__hypergeometric, only: f21

    implicit none (type, external)

    integer, intent(in) :: lambda
      !! The dipole (1), quadrupole (2), ..
    integer, intent(in) :: ltarg
      !! The largest value lmax of ((l+λ,l), l=0, lmax)
    real(dp), intent(in) :: ki, kf
    logical, intent(in), optional :: swapl
    complex(dp) :: res(0:ltarg)

    logical :: swapl_
    integer :: eta
      !! η = λ+1 means something very different here than in the Coulomb case
    integer :: li, lf, li_, lf_
    real(dp) :: a, b, c, x

    swapl_ = .false. ; if(present(swapl)) swapl_ = swapl

    if(G%TARGCHARGE .ne. 0) call die("Called Mborn with a nonzero target charge")

    res = (0.0_dp, 0.0_dp)

    eta = lambda + 1

    ! -- For the dipole case, we only want li = lf ± λ
    do li=0,ltarg

      lf = li + lambda

      ! -- because we want M_{ll'}(k,k') and M_{l'l}(k,k'). We cannot just
      !    swap ki and kf for this integral
      li_ = merge(lf, li, swapl)
      lf_ = merge(li, lf, swapl)

      a = real(li_+lf_-eta+3,kind=dp)/2.0_dp
      b = real(lf_-li_-eta+2,kind=dp)/2.0_dp
      c = real(lf_,kind=dp) + 3.0_dp/2.0_dp
      x = (kf/ki)**2

      res(li)%re = pi/(2.0_dp**eta) * kf**lf_ / ki**(lf_-eta+3) &
                 * gamma(a)/gamma(1._dp-b)/gamma(c)           &
                 * f21(a, b, c, x)

    enddo

  end function Mborn_array

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  function Mcoul_array(lambda, ltarg, ki, kf, swapl) result(res)
    !! Calculate all integrals \(M^\lambda_{λ,0}\) to (M^\lambda_{l_text[targ]+λ,l_\text{targ}}) via recursion formula
    !!   "Study of Nuclear Structure by Electromagnetic Excitation with Accelerated Ions" by K. Alder, A. Bohr, T. Huus,
    !!   B. Mottelson, and A. Winther, equation II B.70–71 on page 453 for the λ=1 case.

    use rotex__kinds,  only: dp, qp
    use rotex__system, only: die, stderr

    implicit none (type, external)

    integer, intent(in) :: lambda
    integer, intent(in) :: ltarg
      !! The truncating value l of (l+λ,l) for the recursion.
    real(dp), intent(in) :: ki, kf
    logical, intent(in), optional :: swapl
    complex(dp) :: res(0:ltarg)
    complex(dp) :: M0, M1
    complex(qp) :: resqp(0:ltarg)

    logical :: swapl_
    integer :: l
    real(dp) :: kf_, ki_
    real(qp) :: etai, etaf

    ! -- because we want to calcualte M_{ll'}(k,k') and M_{l'l}(k,k').
    !    We can just swap ki and kf for this integral.
    ki_ = merge(kf, ki, swapl)
    kf_ = merge(ki, kf, swapl)

    etai = real(-G%TARGCHARGE/ki_, kind = qp)
    etaf = real(-G%TARGCHARGE/kf_, kind = qp)

    if(ltarg .lt. 2) then
      if(swapl) then
        res = [( Mcoul(l, ki_, kf_), l=0, ltarg )]
      else
        res = [( Mcoul(l, kf_, ki_), l=0, ltarg )]
      endif
      block
        use ieee_arithmetic, only: isnan => ieee_is_nan
        if(any(isnan(abs(res)))) then
          write(stderr, *)
          write(stderr, '(6(A15))') "λ", "l", "kf", "ki", "ηf", "ηi"
          write(stderr, '(2I15, 4e15.6)') lambda, ltarg, kf_, ki_, etaf, etai
          write(stderr, *) res
          call die("NaNs detected")
        endif
      end block
      return
    endif

    ! -- starting values
    M0 = Mcoul(0, kf_, ki_)
    M1 = Mcoul(1, kf_, ki_)
    resqp(0) = cmplx(M0%re, M0%im, kind = qp)
    resqp(1) = cmplx(M1%re, M1%im, kind = qp)

    select case(lambda)
    case(1)

      do l=2, ltarg
        resqp(l) = -(y1(l, etai, etaf)*resqp(l-2) + y2(l, etai, etaf)*resqp(l-1)) / y3(l, etai, etaf)
      enddo

    case default

      call die("λ =/= 1 unsupported")

    end select

    res = cmplx(resqp % re, resqp % im, kind = dp)

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  contains
  ! ------------------------------------------------------------------------------------------------------------------------------ !

    ! -- The recursion terms
    impure elemental function y1(l, etai, etaf) result (res)
      implicit none (type, external)
      integer, intent(in)  :: l
      real(qp), intent(in) :: etai, etaf
      real(qp) :: res
      complex(qp), parameter :: im = cmplx(0,1,kind=qp)
      res = 2*etai*etaf * abs(l-1+im*etaf)*abs(l+im*etai)
    end function y1
    impure elemental function y2(l, etai, etaf) result (res)
      implicit none (type, external)
      integer, intent(in)  :: l
      real(qp), intent(in) :: etai, etaf
      real(qp) :: res
      res = -4*etai**2*etaf**2 - l*(2*l+1)*etai**2 - l*(2*l-1)*etaf**2
    end function y2
    impure elemental function y3(l, etai, etaf) result (res)
      implicit none (type, external)
      integer, intent(in)  :: l
      real(qp), intent(in) :: etai, etaf
      real(qp) :: res
      complex(qp), parameter :: im = cmplx(0,1,kind=qp)
      res = 2*etai*etaf * abs(l+im*etaf) * abs(l+1+im*etai)
    end function y3

  end function Mcoul_array

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine xtrapolate_cb_xs(Ei_xtrap, Ethresh, nE_xtrap, Eel, xs_pcb, xs_tcb)
    !! Extrapolate an excitation/de-excitation cross section to its excitation threshold.
    !! Use a power law scaling by fitting the first 2 points.
    !! Re-allocates Eel and xs to contain the extrapolated values

    use rotex__kinds,     only: dp
    use rotex__system,    only: die
    use rotex__functions, only: logrange

    implicit none (type, external)

    real(dp), intent(in) :: Ei_xtrap
      !! Extrapolate down Ethresh + Ei_xtrap
    real(dp), intent(in) :: Ethresh
      !! The excitation threshold
    integer, intent(in) :: nE_xtrap
      !! Number of extrapolation energies
    real(dp), intent(inout), allocatable :: Eel(:)
      !! On input, the electron energy grid.
      !! On output, the electron energy grid with extrapolated energies prepended
    real(dp), intent(inout), allocatable :: xs_pcb(:)
      !! On input, the partial Coulomb-Born cross sections.
      !! On output, the partial Coulomb-Born cross sections with extrapolated cross sections prepended
    real(dp), intent(inout), allocatable :: xs_tcb(:)
      !! On input, the total Coulomb-Born cross sections.
      !! On output, the total Coulomb-Born cross sections with extrapolated cross sections prepended

    real(dp) :: Apcb, Atcb, ppcb, ptcb
    real(dp), allocatable :: Eel_pre(:), xs_pre_pcb(:), xs_pre_tcb(:)

    if(nE_xtrap .le. 1) call die("NE_XTRAP must be > 1")

    ! -- extrapolated energies
    Eel_pre = logrange(Ethresh + Ei_xtrap, Eel(1), nE_xtrap, inclast = .false.)

    ! -- get power law fit parmeters for σPCB and σTCB
    call simple_powerlaw_fit(Apcb, ppcb, Eel(1:2), xs_pcb(1:2))
    call simple_powerlaw_fit(Atcb, ptcb, Eel(1:2), xs_tcb(1:2))

    ! -- extrapolate cross sections
    xs_pre_pcb = Apcb*Eel_pre**(-ppcb)
    xs_pre_tcb = Atcb*Eel_pre**(-ptcb)

    ! -- new energy grid
    Eel     = [Eel_pre, Eel]

    ! -- prepend extrapolated cross sections
    xs_pcb  = [xs_pre_pcb, xs_pcb]
    xs_tcb  = [xs_pre_tcb, xs_tcb]

  end subroutine xtrapolate_cb_xs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure subroutine simple_powerlaw_fit(A, p, E, xs)
    !! Just fit the first two points cross sections to a power law.
    !! Preserves continuity in the derivative. Also make sure the
    !! cross section does not diverge faster than 1/E as E->0⁺
    use rotex__kinds, only: dp
    implicit none (type, external)
    real(dp), intent(out) :: A, p
    real(dp), intent(in) :: E(2), xs(2)
    p = -log(xs(2)/xs(1)) / log(E(2)/E(1))
    p = min(p, 1.0_dp)
    A = xs(1) * E(1)**p
  end subroutine simple_powerlaw_fit

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine powerlaw_fit(A, p, Eel, xs, nfit)
    !! Power-law extrapolatioln to threshold for σCB using
    !! the first nfit available points, assuming:
    !!   σ(E) = A * (E-Ethresh)^(-p)
    !! Returns A and p

    use rotex__kinds,  only: dp
    use rotex__system, only: stderr, die
    use rotex__arrays, only: size_check

    implicit none (type, external)

    real(dp), intent(out) :: A, p
      !! Power law fit parameters
      ! real(dp), intent(in) :: Ethresh
    real(dp), intent(in)  :: Eel(:)
      !! Electron energy
    real(dp), intent(in)  :: xs(:)
      !! Cross sections
    integer,  intent(in)  :: nfit
      !! Number of points to sample/fit

    integer :: i, n
    real(dp) :: sx, sy, sxx, sxy, denom, slope, intercept, logx, logy
    n = size(Eel, 1)
    call size_check(xs, n, "XS")

    if(n .lt. nfit) then
      write(stderr, '("N: ", I0)') n
      write(stderr, '("NFIT: ", I0)') nfit
      call die("Number of fit points is larger than the number of available points !")
    endif

    ! -- least squares minimization of the linearized equation
    sx = 0  ! Σ xi
    sy = 0  ! Σ yi
    sxx = 0 ! Σ xi²
    sxy = 0 ! Σ xi*yi
    do i=1, n
      ! logx = log(Eel(i)-Ethresh)
      logx = log(Eel(i))
      logy = log(xs(i))
      sx = sx + logx
      sy = sy + logy
      sxx = sxx + logx*logx
      sxy = sxy + logx*logy
    enddo

    !                           y = b     +   m  x
    ! -- σ = A (ΔE)*(-p) => ln(σ) = ln(A) + (-p) ln(ΔE)
    !    intercept -------------------^     ^
    !    slope -----------------------------^
    denom     = n*sxx - sx*sx
    slope     = (n*sxy - sx*sy) / denom
    intercept = (sy - slope*sx) / real(n, kind=dp)
    p = -slope
    A = exp(intercept)

  end subroutine powerlaw_fit

! ================================================================================================================================ !
end module rotex__CBXS
! ================================================================================================================================ !
