! ================================================================================================================================ !
module rotex__MQDTXS
  !! Routines to calculate cross sections with MQDT + S-matrix

  implicit none

  private

  public :: get_smat_probs

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_smat_probs( &
        total_energy_grid           &
      , prob                        &
      , transitions                 &
      , nmin                        &
      , nmax                        &
      , smat_j                      &
      , jmin                        &
      , jmax                        &
      , channels_j                  &
      , channels_tot                &
      , spin_isomer_kind            &
      , symaxis                     &
    )
    !! Given a rotationally resolved S-matrix, calculate rotational (de-)excitation
    !! cross section probabilities for the supplied transitions.

    use rotex__types,      only: dp, rvector_type, asymtop_rot_channel_l_vector_type, asymtop_rot_channel_l_type &
                               , asymtop_rot_channel_type, cmatrix_type, n_states_type, asymtop_rot_transition_type &
                               , operator(.ne.), operator(.eq.), operator(.isin.) &
                               , trim_channel_l, get_channel_index
    use rotex__arrays,     only: append, size_check, realloc
    use rotex__symmetry,   only: is_spin_forbidden
    use rotex__system,     only: die, stdout, stderr
    use rotex__constants,  only: pi
    use rotex__characters, only: i2c => int2char
#ifdef USE_FORBEAR
    use rotex__progress,   only: progressbar_type
#endif
#ifdef USE_OPENMP
    use omp_lib, only: omp_get_thread_num
#endif

    implicit none

    real(dp), intent(in) :: total_energy_grid(:)
      !! The total energy grid on which the S-matrix will be evaluated
    type(rvector_type), intent(out), allocatable :: prob(:)
      !! Probability at each pair of channels (n,N,Ka,Kc) ←→ (n',N',Ka',Kc')
    type(asymtop_rot_transition_type), intent(inout), allocatable :: transitions(:)
      !! Array of transitions that will be considered for (de-)excitation
    integer, intent(in) :: nmin, nmax
      !! Min/max values of N to consider for excitation calculations
    type(cmatrix_type), intent(in) :: smat_j(jmin:jmax)
      !! Array of S-matrix sub-blocks for each J
    integer, intent(in) :: jmin, jmax
      !! Min/max values of J to consider for the S-matrix subblocks
    type(asymtop_rot_channel_l_vector_type), intent(in) :: channels_j(jmin:jmax)
      !! Contains the array of channels for each J
    type(asymtop_rot_channel_l_type), intent(in) :: channels_tot(:)
      !! Contains the array of channels ∀ J
    integer, intent(in) :: spin_isomer_kind
      !! What kinda spin symmetry we need to respect
    character(1), intent(in) :: symaxis
      !! The symmetry axis of the target

#ifdef USE_FORBEAR
    integer, parameter :: NDOTS2PRINT = 5
      !! Progress bar update percentage
    integer  :: iprogress, iprogress_last
    real(dp) :: rprogress, rprogress_inc
    type(progressbar_type) progressbar
#endif

    logical, allocatable :: keep_transition_mask(:)
    integer  :: ithread
    integer  :: j, ni, nf
    integer  :: neleclo, nlo, kalo, kclo
    integer  :: nelecup, nup, kaup, kcup
    integer  :: ichan, fchan, itrans
    integer  :: ie, ne, nchans_j, nchans_tot, nopen, nclosed
    real(dp) :: etot, elo, eup
    real(dp) :: sum_ll
    real(dp) :: prob_term
    real(dp),    allocatable :: beta(:)
    integer, allocatable :: indices_lo(:), indices_up(:)
    complex(dp), allocatable :: s(:,:)
    complex(dp), allocatable :: q(:)
    complex(dp), allocatable :: sphys(:,:)
    character(22) :: prefix_string
    type(asymtop_rot_channel_type) :: lo, up
    type(asymtop_rot_transition_type) :: transition
    type(asymtop_rot_channel_l_type), allocatable :: channels_this_j(:)

    nchans_tot = size(channels_tot, 1)

    ! -- build combinations of states (without l), de-excitations will be handled by symmetry
    do ichan = 1, nchans_tot
      Ni = channels_tot(ichan) % N
      if(Ni .lt.  Nmin) cycle
      if(Ni .gt.  Nmax) cycle
      lo = trim_channel_l(channels_tot(ichan))
      do fchan = ichan+1, nchans_tot
        Nf = channels_tot(fchan) % N
        if(Nf .lt. Nmin) cycle
        if(Nf .gt. Nmax) cycle
        up = trim_channel_l(channels_tot(fchan))
        ! -- skip elastic pairs
        if(lo .eq. up) cycle
        ! -- skip de-excitations for now. These shoud not show up here anyway; they'll be handled symmetrically
        !    when excitations are considered
        if(lo % E .ge. up % E) cycle
        !  -- respect ortho/para symmetry if applicable (returns true if theres nothing to respect)
        if(is_spin_forbidden(lo, up, spin_isomer_kind, symaxis)) cycle
        transition = asymtop_rot_transition_type(lo = lo,  up = up)
        ! -- only append transitions uniquely
        if(allocated(transitions) .eqv. .false.) then
          call append(transitions, transition)
          cycle
        endif
        if(transition .isin. transitions) cycle
        call append(transitions, transition)
      enddo
    enddo

    ! -- keep track of which transitions to keep based on the corresponding cross sections
    allocate(keep_transition_mask(size(transitions, 1)), source = .true.)

    call size_check(channels_J, Jmax - Jmin + 1, "CHANNELS_J")

    ne = size(total_energy_grid, 1)

    ! -- allocate and initialize the probabilities arrays for each transition
    allocate(prob(size(transitions, 1)))
    do itrans = 1, size(transitions, 1)
      allocate(prob(itrans) % vec(ne))
      prob(itrans) % vec(:) = 0
    enddo

    ! -- update user
    write(*,*)
    write(stdout, '(A)') "Performing closed-channel elimination on the S-matrix and accumulating cross section probabilities"

    ! -- Σ_J
    write(stdout, '(2X, 2(A5, " /"), A5)') "Jmin", "J", "Jmax"
    do J=Jmin,Jmax

      write(prefix_string, '(2X, 2(I5," /"),I5,X)') Jmin, J, Jmax
#ifdef USE_FORBEAR
      call progressbar % initialize( &
          filled_char_string = "|" &
        , empty_char_string = " " &
        , bracket_left_string = "[" &
        , prefix_string = prefix_string &
        , suffix_string = "] " &
        , add_progress_percent = .true. &
      )
      call progressbar % start
      call progressbar % update(current = 0.0_dp)
#else
      write(stdout, '(A)') prefix_string
#endif

      ! -- get the channels, S-matrix, and f/g normalization coefficients for this J
      channels_this_J = channels_J(J) % channels

      S = smat_J(J) % mtrx
      nchans_J = size(channels_this_J, 1)

      call size_check(S, [nchans_J, nchans_J], "S")

#ifdef USE_FORBEAR
      rprogress    = 0.0_dp
      iprogress_last = 0
      rprogress_inc = 1.0_dp / real(ne, kind = dp)
#endif

      ! -- loop over the total enrgy grid
      !$omp parallel default(none) &
      !$omp& shared(ne, channels_this_J, S, prob, transitions, J, nchans_J, total_energy_grid&
#ifdef USE_FORBEAR
      !$omp&   , progressbar, rprogress, rprogress_inc, iprogress, iprogress_last)&
#else
      !$omp&   )&
#endif
      !$omp& private(ie, Etot, Sphys, beta, nopen, nclosed, lo, up, indices_lo, indices_up, itrans&
      !$omp&   , q,  neleclo, nelecup, Nlo, Nup, Kalo, Kaup, Kclo, Kcup, Elo, Eup, prob_term &
      !$omp&   , sum_ll, ithread)
#ifdef USE_OPENMP
      ithread = omp_get_thread_num()
#else
      ithread = 0
#endif
      !$omp do schedule(static)
      do ie=1,ne
      ! do concurrent (ie=1:ne)

        call realloc(q, nchans_J)

        Etot = total_energy_grid(ie)
        q    = fg_norm_coeff_q(channels_this_J, Etot)

        nclosed = count(channels_this_J % E .gt. Etot)
        nopen   = count(channels_this_J % E .le. Etot)
        if(nopen + nclosed .ne. nchans_J) call die("Number of opened and closed channels does not add to the number of channels !")

        call realloc(beta, nclosed)
        call realloc(Sphys, nopen, nopen)
        beta = 0
        Sphys = 0

        call CCEP(S, channels_this_J, q, Etot, Sphys, beta, nopen, nclosed)

#ifdef USE_FORBEAR
        !$omp atomic
        rprogress = rprogress + rprogress_inc

        ! -- update progress
        update_progress: if(ithread .eq. 0) then
          iprogress = floor(rprogress * 100)
          if(iprogress .eq. iprogress_last) exit update_progress
          if(iprogress .eq. 100) exit update_progress
          iprogress_last = iprogress
          call progressbar % update(current = rprogress)
        endif update_progress
#endif

        ! -- loop over the pairs of states for excitation, accumulate probabilities in each transition
        !    for this J
        do itrans = 1, size(transitions, 1)
          lo = transitions(itrans) % lo
          up = transitions(itrans) % up
          ! -- make sure both channels in this transition are actually in this J block
          if((lo .isin. channels_this_J(1:nopen)) .eqv. .false.) cycle
          if((up .isin. channels_this_J(1:nopen)) .eqv. .false.) cycle
          neleclo = lo % nelec
          Nlo     = lo % N
          Kalo    = lo % Ka
          Kclo    = lo % Kc
          nelecup = up % nelec
          Nup     = up % N
          Kaup    = up % Ka
          Kcup    = up % Kc
          Elo     = lo % E
          Eup     = up % E
          if(Elo .ge. Etot) cycle
          if(Eup .gt. Etot) cycle
          ! -- get the indices of the S-matrix channels (with l) that match
          indices_lo = pack([(ichan,ichan=1,nopen)], lo .eq. channels_this_J(1:nopen))! .AND. channels_this_J % E .lt. Etot)
          indices_up = pack([(ichan,ichan=1,nopen)], up .eq. channels_this_J(1:nopen))! .AND. channels_this_J % E .lt. Etot)
          ! -- Σ_{ll'} |Sphys_{il,i'l'}|²
          prob_term = (2*J+1) * sum(abs(Sphys(indices_lo, indices_up))**2)

          prob(itrans) % vec(ie) = prob(itrans) % vec(ie) + prob_term

          if(prob_term .ge. 0) cycle

          ! -- error
          write(stderr, '(A)') "❌"
          write(stderr, '(2X, A, 3E15.6)') "Elo, Eup, Etot: ", Elo, Eup, Etot
          write(stderr, '(2X, A, 4I4)') "Lower state nelec, N, Ka, Kc: ", neleclo, Nlo, Kalo, Kclo
          write(stderr, '(2X, A, 4I4)') "Upper state nelec, N, Ka, Kc: ", nelecup, Nup, Kaup, Kcup
          write(stderr, '(2X, A20, E15.6)') "(2J+1) Σ_{ll'} |S_{ll'}|²", sum_ll
          call die("Negative probability from the S-matrix !")

        enddo

      enddo
      !$omp end do
      !$omp end parallel

#ifdef USE_FORBEAR
      call progressbar % update(current = 1.0_dp)
#endif

    enddo

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  contains
  ! ------------------------------------------------------------------------------------------------------------------------------ !

    ! ---------------------------------------------------------------------------------------------------------------------------- !
    pure function fg_norm_coeff_q(channels, Etot) result(res)
      !! Given an array of channels, determine the factor B for each channel at a particular total energy E
      use rotex__system,    only: die
      use rotex__constants, only: pi
      implicit none
      type(asymtop_rot_channel_l_type), intent(in) :: channels(:)
        !! Array of channels for which we want to get the factor B
      real(dp), intent(in) :: Etot
        !! The total energy
      complex(dp), allocatable :: res(:)
      integer :: i, l, n, iq
      real(dp) :: EE, Echan
      n = size(channels, 1)
      allocate(res(n))
      res = 1
      do i = 1, n
        iq = channels(i) % iq
        l  = channels(i) % l
        if(iq .eq. 4) cycle
        if(iq .ne. 0) call die("iq cannot be different from 4 and 0")
        Echan = channels(i) % E
        EE = Etot - Echan
        res(i) = sqrt(A_coulomb(2*EE, l))
        ! -- extra factor for closed channels only
        if(EE .gt. 0) cycle
        res(i) = res(i) / sqrt(1._dp - exp(-2*pi/sqrt(2*EE)))
      enddo
    end function fg_norm_coeff_q

    ! ---------------------------------------------------------------------------------------------------------------------------- !
    pure elemental function A_coulomb(e, l) result(res)
      !! Calcualte the factor A for the Coulomb functions (Seaton, 2002, Comp. Phys. Comm.)
      implicit none
      real(dp), intent(in) :: e
      integer,  intent(in) :: l
      real(dp) :: res
      integer :: n
      res = product( [( 1+n*n*e, n=0, l )] )
    end function A_coulomb

  end subroutine get_smat_probs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine CCEP(S, channels, q, Etot, Sphys, beta, nopen, nclosed)
    !! Given the S-matrix, its basis of channels, the total energy E, and the number of open/closed channels,
    !! carry out the MQDT Closed-Channel Elimination Procedure to obtain the nopen x nopen physical S matrix

    use rotex__types,      only: dp, asymtop_rot_channel_l_type
    use rotex__system,     only: die
    use rotex__linalg,     only: zgesv, right_divide
    use rotex__arrays,     only: size_check
    use rotex__constants,  only: im, pi
    use rotex__characters, only: i2c => int2char

    implicit none

    complex(dp), intent(in) :: S(:,:)
      !! The n x n S-matrix before channel elimination
    type(asymtop_rot_channel_l_type), intent(in) :: channels(:)
      !! The rotational channels
    complex(dp), intent(in) :: q(:)
      !! The f/g normalization factors for each channel
    real(dp), intent(in) :: Etot
      !! The total energy
    complex(dp), intent(inout) :: Sphys(:,:)
      !! The nopen x nopen S-matrix after channel elimination
    real(dp), intent(inout) :: beta(:)
      !! The quantum defects for this channel at the given total energy
    integer, intent(in) :: nopen
      !! The number of open channels
    integer, intent(in) :: nclosed
      !! The number of closed channels

    integer                  :: ichan, iopen, iclosed, nchan
    real(dp)                 :: Echan
    complex(dp)              :: numer, denom
    complex(dp), allocatable :: sinbeta(:), cosbeta(:)
    complex(dp), allocatable :: qmat(:,:), qmatinv(:,:)
    complex(dp), allocatable :: Soo(:,:), Soc(:,:), Sco(:,:), Scc(:,:)
    complex(dp), allocatable :: Sce(:,:)
    complex(dp), allocatable :: A(:,:), B(:,:)
    ! -- LAPACK specific
    integer :: info
    integer, allocatable :: ipiv(:)

    nchan = nopen + nclosed
    call size_check(S,        [nchan,nchan], "S")
    call size_check(q,        nchan,         "Q")
    call size_check(beta,     nclosed,       "BETA")
    call size_check(Sphys,    [nopen,nopen], "SPHYS")
    call size_check(channels, nchan,         "CHANNELS")

    ! -- build β
    do concurrent (iclosed=1:nclosed)
      ichan = iclosed + nopen
      Echan = channels(ichan) % E
      beta(iclosed) = pi / sqrt(2*(Echan - Etot))
    enddo

    ! -- no CCEP if all channels are open
    if(nclosed .eq. 0) then
      Sphys =  S
      return
    endif

    Soo = S(1:nopen,       1:nopen) ; Soc  = S(1:nopen,      nopen+1:nchan)
    Sco = S(nopen+1:nchan, 1:nopen) ; Scc  = S(nopen+1:nchan,nopen+1:nchan)

    call size_check(Soo, [nopen,   nopen],   "SOO")
    call size_check(Soc, [nopen,   nclosed], "SOC")
    call size_check(Sco, [nclosed, nopen],   "SCO")
    call size_check(Scc, [nclosed, nclosed], "SCC")

    sinbeta = sin(beta)
    cosbeta = cos(beta)

    do concurrent (iclosed=1:nclosed)
      ichan = iclosed + nopen
      numer = q(ichan)*q(ichan)*cosbeta(iclosed) - im*sinbeta(iclosed)
      denom = q(ichan)*q(ichan)*cosbeta(iclosed) + im*sinbeta(iclosed)
      Scc(iclosed, iclosed) = Scc(iclosed, iclosed) - numer / denom
    enddo

    allocate(ipiv(nclosed))
    call zgesv(nclosed, nopen,  Scc,  nclosed, ipiv, Sco, nclosed, info)

    if(info .ne. 0) call die("ZGESV exited with INFO = " // i2c(info))

    Sce = matmul(Soc, Sco)
    Sce = Soo - Sce

    deallocate(Scc, Sco, Soc, Soo)
    deallocate(ipiv)

    ! -- leave if there is not extra f/g function norm stuff to deal with
    if(all(q .eq. cmplx(1, kind = dp))) then
      Sphys = Sce
      return
    endif

    allocate(qmat(nopen, nopen))
    allocate(qmatinv(nopen, nopen))
    qmat    = 0
    qmatinv = 0

    do concurrent (iopen=1:nopen)
      ichan = iopen
      qmat(iopen, iopen)    = q(ichan)
      qmatinv(iopen, iopen) = 1/q(ichan)
    enddo

    A = (qmatinv - qmat) + matmul(qmatinv + qmat, Sce)
    B = (qmatinv + qmat) + matmul(qmatinv - qmat, Sce)

    ! -- S = AB⁻¹
    Sphys = right_divide(A, B)

  end subroutine CCEP

! ================================================================================================================================ !
end module rotex__MQDTXS
! ================================================================================================================================ !
