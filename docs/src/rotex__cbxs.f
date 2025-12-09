! ================================================================================================================================ !
module rotex__CBXS
  !! Routines to calculate cross sections in the Coulomb-Born approximation

  implicit none

  private

  public :: get_einsta_only
  public :: get_CB_xs_asym
  public :: M
  public :: xtrapolate_cb_xs

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_einsta_only( einsta, nlo, nup, elo, eup, eigveclo, eigvecup &
                                   , use_CDMS, do_dipole, do_quadrupole, dipole_moments, quadrupole_moments )
    !! Calculate only the Einstein A coefficeints for a transition

    use rotex__types,      only: dp
    use rotex__system,     only: die
    use rotex__constants,  only: invc
    use rotex__functions,  only: istriangle

    implicit none

    integer, intent(in) :: nlo
      !! the angular momentum quantum number \(N\)
    integer, intent(in) :: nup
      !! the angular momentum quantum number \(N'\)
    real(dp), intent(in) :: elo
      !! The energy (hartrees) of the initial state
    real(dp), intent(in) :: eup
      !! The energy (hartrees) of the final state
    real(dp), intent(in) :: eigveclo(:)
      !! Eigenvector of the initial state in the basis of symmetric top wavefunctions; the coefficients
      !! \(c^{N,\tau}_K)
    real(dp), intent(in) :: eigvecup(:)
      !! Eigenvector of the final state in the basis of symmetric top wavefunctions; the coefficients
      !! \(c^{N',\tau'}_{K})
    real(dp), intent(inout) :: einsta
      !! The Einstein coefficient for the transition
    logical, intent(in) :: use_CDMS
      !! Whether to calculate the Einstein A coefficients ourselves (.true.) or to use values obtained
      !! from the CDMS catalogue (.false.)
    logical, intent(in) :: do_dipole, do_quadrupole
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
        if(do_dipole .eqv. .false.) cycle multipoles
        multipole_moments => dipole_moments
      case(2)
        if(do_quadrupole .eqv. .false.) cycle multipoles
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
        if((use_CDMS .eqv. .false.) .OR. einsta .eq. 0) then
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
      case(2)
        einsta = einsta + am_summation * 2._dp/15._dp * (omega * invc)**5 * (2*nlo+1)
      end select

    enddo multipoles

  end subroutine get_einsta_only

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_CB_xs_asym( energies, sigma, Z &
                               , N, Np, E, Ep, eigvec, eigvecp &
                               , einsta, use_CDMS &
                               , do_dipole, do_quadrupole &
                               , dipole_moments, quadrupole_moments &
                               , analytic_total_cb, lmax)
                               ! , analytic_total_cb, lmax, atol, rtol)
    !! Calculate the excitation and de-excitation cross sections (xs) for an asymmetric top up.
    !! The sum over partial waves is either truncated to \(l_\text{max}\) or determined analytically.
    !! The summation over the angular momentum components, multipole terms, and partial waves are separable
    !! for each value of λ as summation = (sum over angular momentum and multipole moments) \(\times\) (sum overpartial waves).

    use rotex__types,      only: dp
    use rotex__system,     only: die
    use rotex__constants,  only: invc
    use rotex__functions,  only: istriangle

    implicit none

    real(dp), intent(in) :: energies(:)
      !! Array of scattering energies to consider
    real(dp), intent(out), allocatable :: sigma(:)
      !! Cross sections calculated on a grid of scattering energies for \(Nτ \rightarrow N'τ'\)
      !! Returned with the same size as `energies`
    integer, intent(in) :: Z
      !! Target charge
    integer, intent(in) :: N
      !! the angular momentum quantum number \(N\)
    integer, intent(in) :: Np
      !! the angular momentum quantum number \(N'\)
    real(dp), intent(in) :: E
      !! The energy (hartrees) of the initial state
    real(dp), intent(in) :: Ep
      !! The energy (hartrees) of the final state
    real(dp), intent(in) :: eigvec(:)
      !! Eigenvector of the initial state in the basis of symmetric top wavefunctions; the coefficients
      !! \(c^{N,\tau}_K)
    real(dp), intent(in) :: eigvecp(:)
      !! Eigenvector of the final state in the basis of symmetric top wavefunctions; the coefficients
      !! \(c^{N',\tau'}_{K})
    real(dp), intent(inout) :: einsta
      !! The Einstein coefficient for the transition
    logical, intent(in) :: use_CDMS
      !! Whether to calculate the Einstein A coefficients ourselves (.true.) or to use values obtained
      !! from the CDMS catalogue (.false.)
    logical, intent(in) :: do_dipole, do_quadrupole
    complex(dp), intent(in), target :: dipole_moments(3)
      !! The spherical dipole moments
    complex(dp), intent(in), target :: quadrupole_moments(5)
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

    if(Z .eq. 0) call die("Coulomb-Born method should not be called on neutral targets !")

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
        if(do_dipole     .eqv. .false.) cycle multipoles
        multipole_moments => dipole_moments
      case(2)
        if(do_quadrupole .eqv. .false.) cycle multipoles
        multipole_moments => quadrupole_moments
      case default
        call die("Somehow, lambda is neither 1 nor 2")
      end select

      ! -- enforce 3j rules to avoid extra computation
      if( .false. .eqv. istriangle(N, Np, lambda) ) cycle multipoles

      ! -- do the summation over angular momenta and multipole components
      !    or use the CDMS Einstein A coeffs for the dipole case
      select case(lambda)
      case(1)
        if((use_CDMS .eqv. .false.) .OR. einsta .eq. 0) then
          ! -- calculate the am_summation ourselves
          call multipole_am_summation_asym(am_summation, lambda, N, Np, multipole_moments, eigvec, eigvecp)
        else
          ! -- use pre-existing Einstein A coeff from the CDMS read
          am_summation = einsta * 3._dp/4._dp * (omega * invc)**(-3) / (2*N+1)
        endif
      case(2)
        call multipole_am_summation_asym(am_summation, lambda, N, Np, multipole_moments, eigvec, eigvecp)
      end select

      nullify(multipole_moments)

      if(abs(am_summation) .lt. 1e-16) cycle multipoles

      if(am_summation .lt. 0) call die("The angular momentum summation is negative, which is not physical !")

      ! -- if the summation is basically 0 for this multipole term, check the next one
      ! if(am_summation .lt. 1e-16_dp) cycle multipoles

      ! -- do the summation over partial waves
      if(analytic_total_cb(lambda) .eqv. .true.) then
        call get_infinite_pwsum(pw_summation, energies, E, Ep, lambda, Z)
      else
        call get_truncated_pwsum(pw_summation, energies, E, Ep, lambda, lmax, Z)
      endif

      if(any(pw_summation .lt. 0)) &
        call die("The partial wave summation is negative for at least one energy, which is not physical !")

      sigma  = sigma + pw_summation * am_summation / (2*lambda+1)

      select case(lambda)
      case(1)
        if((use_CDMS .eqv. .false.) .OR. einsta .eq. 0) then
          einsta = einsta &
                 + am_summation * 4._dp/3._dp * (omega * invc)**3 * (2*N+1)
        endif
      case(2)
        einsta = einsta + am_summation * 2._dp/15._dp * (omega * invc)**5 * (2*N+1)
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
    use rotex__types,  only: dp
    use rotex__system, only: die
    use rotex__wigner, only: wigner3j

    implicit none

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
    real(dp),    intent(in)  :: eigvec(:)
      !! The eigenvector for the initial state Nτ
    real(dp),    intent(in)  :: eigvecp(:)
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
                 * eigvecp(i_K1p) * eigvec(i_K1)                          &
                 * eigvec(i_K2)   * eigvecp(i_K2p)                        &
                 * wigner3j(2*N, 2*Np, 2*lambda, -2*K1,  2*K1p, -2*mu1)   &
                 * wigner3j(2*N, 2*Np, 2*lambda, -2*K2,  2*K2p, -2*mu2)
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
  subroutine get_infinite_pwsum(summation, energies, E, Ep, lambda, Z)
    !! Return the analytic expression for the partial wave sum in the limit \(l\to\infty\)
    use rotex__types,          only: dp
    use rotex__system,         only: die, stderr
    use ieee_arithmetic,       only: ieee_is_nan
    use rotex__constants,      only: pi, im
    use rotex__functions,      only: expm1
    use rotex__hypergeometric, only: f21

    implicit none

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
    integer, intent(in) :: Z
      !! Target charge

    logical :: flag
    integer :: ie, nE
    integer :: iemin
    real(dp) :: k, kp
    real(dp) :: eta, etap
    real(dp) :: dE
    real(dp) :: x
    complex(dp), parameter :: onez = (1, 0)

    if(lambda .ne. 1) call die("analytic partial wave sum not available for lambda =/= 1")

    nE = size(energies, 1)
    dE = Ep - E

    flag = .false.
    summation = 0
    iemin = findloc(energies .gt. dE, .true., 1)

    !$omp parallel do default(none) &
    !$omp& shared(energies, summation, E, Ep, dE, lambda, Z, iemin, nE, flag) &
    !$omp& private(ie, k, kp, eta, etap, x)
    do ie = iemin, nE
    ! do concurrent(ie = 1 : nE)

      k    = sqrt(2*energies(ie))
      kp   = sqrt(2*(energies(ie) - dE))
      eta  = -Z/k
      etap = -Z/kp
      x = -4*eta*etap/((etap-eta)**2)

      summation(ie) = -16*pi*pi*pi/(k*k*k*kp) * x &
                    * exp(2*pi*eta) / ( expm1(2*pi*eta) * expm1(2*pi*etap) ) &
                    * real( &
                        f21(im*eta, im*etap, onez, x) &
                      * f21(-im*eta + 1, -im*etap + 1, 2*onez, x), kind = dp &
                    )
                    ! * real(F21(im*eta, im*etap, onez, zx) * F21(-im*eta + 1, -im*etap + 1, 2*onez, zx), kind = dp)

      if(summation(ie) .ge. 0) cycle
      if(flag .eqv. .true.) cycle
      flag = .true.

      block
        use rotex__constants, only: au2ev
        write(stderr, *)
        write(stderr, *) E*au2ev, Ep*au2ev
        write(stderr, *) energies(ie)*au2ev, (energies(ie) - dE)*au2ev
        write(stderr, *) ie, energies(ie), eta, etap, x, 1/(1-x)
        write(stderr, '(A)') "WARN: negative cross sections when evaluating the infinite partial wave sum."
      end block

    enddo
    !$omp end parallel do

  end subroutine get_infinite_pwsum

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine get_truncated_pwsum(summation, energies, E, Ep, lambda, lmax, Z)
    !! Calculate the truncated partial wave sum
    !! \(\sum\limits_{ll'}^{\l_\text{max}} (2l+1)(2l'+1) (l,l,\lambda;0,0,0)^2 |M^\lambda_{ll'}|^2\), given in
    !!   "Electromagnetic Excitation: Theory of Coulomb Excitation with Heavy Ions " by Kurt Alder and Aage Winther, Chapter IX,
    !!    section 2, page 244, equation 14.
    !! The integral given in that text is \(I^\lambda_{ll'}=-4M^\lambda_{ll'}\) in atomic units, based on equations 7 and 11 of the
    !! same section. For values of l larger than 1, a recursion formula is used for the integral \(M^\lambda_{ll'}\), given in
    !!   "Study of Nuclear Structure by Electromagnetic Excitation with Accelerated Ions" by K. Alder, A. Bohr, T. Huus,
    !!    B. Mottelson, and A. Winther, equation II B.70–71 on page 453 for the λ=1 case.
    use rotex__types,     only: dp
    use rotex__system,    only: stderr
    ! use WignerSymbol, only: wigner3j
    use rotex__wigner,    only: wigner3j
    use rotex__constants, only: pi, CB_MINT_IMAG_THRESH
    use rotex__functions, only: istriangle
    ! use rotex__wigner, only: wigner3j => threej

    implicit none

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
    integer, intent(in) :: Z
      !! Target charge

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

    nE = size(energies, 1)

    dE = Ep - E
    flag = .false.
    summation = 0

    ! -- the starting energy index to ensure positive energies
    iemin = findloc(energies .gt. dE, .true., 1)

    ! -- the sum over ll'
    !$omp parallel do default(none) &
    !$omp& shared(energies, summation, E, Ep, dE, lambda, Z, iemin, nE, flag, lmax) &
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
      M1      = Mrecur(lambda, lmax-lambda, kp, k , Z)
      M2      = Mrecur(lambda, lmax-lambda, k , kp, Z)

      if(any(abs(M1%im) .gt. CB_MINT_IMAG_THRESH) .OR. any(abs(M2%im) .gt. CB_MINT_IMAG_THRESH)) then
        write(stderr, *)
        write(stderr, '("WARN: Skipping energy ", I0, ": the Coulomb integrals are nonreal")') ie
        write(stderr, *)
        cycle nrg
      endif

      summation(ie) = 16*pi*kp/k * sum( weights * (abs(M1)**2 + abs(M2)**2) )

      if(summation(ie) .ge. 0)  cycle
      if(flag .eqv. .true.) cycle
      flag = .true.

      write(stderr, *)
      write(stderr, '(A)') "WARN: negative cross sections when evaluating the truncated partial wave sum."

    enddo nrg
   !$omp end parallel do

  end subroutine get_truncated_pwsum

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  impure elemental function M(l, ki, kf, Z) result(res)
    !! Calculates the integral \(M^\lambda_{l+λ,l}\) via the expression given in
    !!   "Electromagnetic Excitation: Theory of Coulomb Excitation with Heavy Ions " by Kurt Alder and Aage Winther, Chapter IX,
    !!    section 2, page 244, equation 14.
    !! for λ = 1, where \(-4M^\lambda_{l+λ,l} = I^λ_{l+λ, l}\). There is an expression for λ = 2 in
    !!   "Study of Nuclear Structure by Electromagnetic Excitation with Accelerated Ions" by K. Alder, A. Bohr, T. Huus,
    !!   B. Mottelson, and A. Winther, but only the dipole is used (at least for now).

    use rotex__types,          only: dp
    use rotex__constants,      only: im, pi
    use rotex__polygamma,      only: lgamma => log_gamma
    use rotex__functions,      only: factorial
    use rotex__hypergeometric, only: f21

    implicit none

    integer, intent(in) :: l
    real(dp), intent(in) :: ki, kf
    integer, intent(in) :: Z
    complex(dp) :: res

    real(dp) :: etai, etaf, deta, x0
    complex(dp) :: c1, c2

    etai = -Z/ki
    etaf = -Z/kf
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

  end function M

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  function Mrecur(lambda, ltarg, ki, kf, Z) result(res)
    !! Calculate all integrals \(M^\lambda_{λ,0}\) to (M^\lambda_{l_text[targ]+λ,l_\text{targ}}) via recursion formula
    !!   "Study of Nuclear Structure by Electromagnetic Excitation with Accelerated Ions" by K. Alder, A. Bohr, T. Huus,
    !!   B. Mottelson, and A. Winther, equation II B.70–71 on page 453 for the λ=1 case.

    use rotex__types,  only: dp, qp
    use rotex__system, only: die, stderr

    implicit none

    integer, intent(in) :: lambda
    integer, intent(in) :: ltarg
      !! The truncating value of l+λ, l for the recursion.
    real(dp), intent(in) :: ki, kf
    integer, intent(in) :: Z
      !! Target charge
    complex(dp) :: res(0:ltarg)
    complex(dp) :: M0, M1
    complex(qp) :: resqp(0:ltarg)

    integer  :: l
    real(qp) :: etai, etaf

    etai = real(-Z/ki, kind = qp)
    etaf = real(-Z/kf, kind = qp)

    if(ltarg .lt. 2) then
      res = [( M(l, kf, ki, Z), l=0, ltarg )]
      block
        use ieee_arithmetic, only: isnan => ieee_is_nan
        if(any(isnan(abs(res)))) then
          write(stderr, *)
          write(stderr, '(6(A15))') "λ", "l", "kf", "ki", "ηf", "ηi"
          write(stderr, '(2I15, 4e15.6)') lambda, ltarg, kf, ki, etaf, etai
          write(stderr, *) res
          call die("NaNs detected")
        endif
      end block
      return
    endif

    ! -- starting values
    M0 = M(0, kf, ki, Z)
    M1 = M(1, kf, ki, Z)
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
      use rotex__types,     only: dp
      use rotex__constants, only: im
      implicit none
      integer, intent(in)  :: l
      real(qp), intent(in) :: etai, etaf
      real(qp) :: res
      res = 2*etai*etaf * abs(l-1+im*etaf)*abs(l+im*etai)
    end function y1
    impure elemental function y2(l, etai, etaf) result (res)
      use rotex__types,     only: dp
      use rotex__constants, only: im
      implicit none
      integer, intent(in)  :: l
      real(qp), intent(in) :: etai, etaf
      real(qp) :: res
      res = -4*etai**2*etaf**2 - l*(2*l+1)*etai**2 - l*(2*l-1)*etaf**2
    end function y2
    impure elemental function y3(l, etai, etaf) result (res)
      use rotex__types,     only: dp
      use rotex__constants, only: im
      implicit none
      integer, intent(in)  :: l
      real(qp), intent(in) :: etai, etaf
      real(qp) :: res
      res = 2*etai*etaf * abs(l+im*etaf) * abs(l+1+im*etai)
    end function y3

  end function Mrecur

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure module subroutine xtrapolate_cb_xs(Ei_xtrap, Ethresh, nE_xtrap, Eel, xs_pcb, xs_tcb)
    !! Extrapolate an excitation/de-excitation cross section to its excitation threshold.
    !! Re-allocates Eel and xs to contain the extrapolated values

    use rotex__kinds,     only: dp
    use rotex__system,    only: die
    use rotex__functions, only: logrange

    implicit none

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

    real(dp), allocatable :: Eel_pre(:), xs_pre_pcb(:), xs_pre_tcb(:)

    if(nE_xtrap .le. 1) call die("NE_XTRAP must be > 1")

    ! -- extrapolated energies
    Eel_pre = logrange(Ethresh + Ei_xtrap, Eel(1), nE_xtrap, inclast = .false.)

    ! -- extrapolated cross sections σ1*E1 = σ2*E2
    xs_pre_pcb  = xs_pcb(1) * Eel(1) / Eel_pre
    xs_pre_tcb  = xs_tcb(1) * Eel(1) / Eel_pre

    ! -- new energy grid
    Eel     = [Eel_pre, Eel]

    ! -- prepend extrapolated cross sections
    xs_pcb  = [xs_pre_pcb, xs_pcb]
    xs_tcb  = [xs_pre_tcb, xs_tcb]

  end subroutine xtrapolate_cb_xs

! ================================================================================================================================ !
end module rotex__CBXS
! ================================================================================================================================ !
