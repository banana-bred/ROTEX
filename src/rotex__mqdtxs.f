! ================================================================================================================================ !
module rotex__MQDTXS
  !! Routines to calculate cross sections with MQDT + S-matrix

  use rotex__globals, only: G, PACKMAT_TRIANGLE
  use rotex__kinds, only: dp
  use rotex__system, only: stdout, stderr, die

  implicit none (type, external)

  private

  public :: get_smat_probs
  public :: get_smat_probs_chunk

! ================================================================================================================================ !
contains
! ================================================================================================================================ !

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_smat_probs( &
        total_energy_grid           &
      , transition_probs            &
      , transitions                 &
      , smat_rot_flat               &
      , smat_eval_energies          &
      , J                           &
      , channels_this_J             &
    )
    !! Given a rotationally resolved S-matrix, calculate rotational (de-)excitation
    !! cross section probabilities for the supplied transitions.

    use rotex__types,      only: prob_vector_type, asymtop_rot_channel_l_vector_type, asymtop_rot_channel_l_type &
                               , asymtop_rot_channel_type, r3carr_type, n_states_type, asymtop_rot_transition_type
    use rotex__channel_ops, only: operator(.ne.), operator(.eq.), operator(.isin.), trim_channel_l, get_channel_index
    use rotex__arrays,     only: append, size_check, realloc, interp_array_at_energy
    use rotex__symmetry,   only: is_spin_forbidden
    use rotex__characters, only: i2c => int2char

    implicit none (type, external)

    real(dp), intent(in) :: total_energy_grid(:)
      !! The total energy grid on which the S-matrix will be evaluated
    type(prob_vector_type), intent(inout), allocatable :: transition_probs(:)
      !! Probability at each pair of channels (n,N,Ka,Kc) ←→ (n',N',Ka',Kc')
    type(asymtop_rot_transition_type), intent(in), allocatable :: transitions(:)
      !! Array of transitions that will be considered for (de-)excitation
    complex(dp), intent(in), contiguous :: smat_rot_flat(:,:)
      !! Flattened array of S-matrix sub-block for this J
    real(dp), intent(in) :: smat_eval_energies(:)
      !! Array of evaluation energies of the S-matrix
    integer, intent(in) :: J
      !! The current J
    type(asymtop_rot_channel_l_type), intent(in) :: channels_this_j(:)
      !! The array of channels for THIS J

    integer :: ne_mat

    ne_mat = size(smat_eval_energies, 1)

    if(G%ALLOW_EDFT_EGRID_OUT_OF_BOUNDS .eqv. .false.) then
      if(minval(total_energy_grid) .lt. smat_eval_energies(1)) then
        call die("TOTAL_ENERGY_GRID extends below the EDFT evaluation grid and strict EDFT interpolation is requested")
      endif
      if(maxval(total_energy_grid) .gt. smat_eval_energies(ne_mat)) then
        call die("TOTAL_ENERGY_GRID extends above the EDFT evaluation grid and strict EDFT interpolation is requested")
      endif
    endif

    call get_smat_probs_chunk(       &
        total_energy_grid            &
      , lbound(total_energy_grid, 1) &
      , ubound(total_energy_grid, 1) &
      , transition_probs             &
      , transitions                  &
      , smat_rot_flat                &
      , smat_eval_energies           &
      , J                            &
      , channels_this_J              &
      , .true., .true.               &
    )

  end subroutine get_smat_probs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  module subroutine get_smat_probs_chunk( &
        total_energy_grid           &
      , ietot0, ietot1              &
      , transition_probs            &
      , transitions                 &
      , smat_rot_flat_chunk         &
      , smat_eval_energies_chunk    &
      , J                           &
      , channels_this_J             &
      , is_first_chunk              &
      , is_last_chunk               &
    )
    !! Given a rotationally resolved S-matrix, calculate rotational (de-)excitation
    !! cross section probabilities for the supplied transitions.

    use rotex__types,      only: asymtop_rot_channel_l_type, prob_vector_type, asymtop_rot_channel_type &
                               , n_states_type, asymtop_rot_transition_type
    use rotex__channel_ops, only: operator(.ne.), operator(.eq.), operator(.isin.), trim_channel_l, get_channel_index
    use rotex__arrays,     only: append, size_check, realloc, interp_array_at_energy, unpackmat_ch
    use rotex__symmetry,   only: is_spin_forbidden
    use rotex__characters, only: i2c => int2char

    implicit none (type, external)

    real(dp),                          intent(in)             :: total_energy_grid(:)
      !! The total energy grid on which the S-matrix will be evaluated
    integer,                           intent(in)             :: ietot0, ietot1
      !! Total energy grid boundary points for the current chunk
    type(prob_vector_type),            intent(inout)          :: transition_probs(:)
      !! Probability at each pair of channels (n,N,Ka,Kc) ←→ (n',N',Ka',Kc')
    type(asymtop_rot_transition_type), intent(in)             :: transitions(:)
      !! Array of transitions that will be considered for (de-)excitation
    complex(dp),                       intent(in), contiguous :: smat_rot_flat_chunk(:,:)
      !! Flattened array of S-matrix sub-block for this J
    real(dp),                          intent(in)             :: smat_eval_energies_chunk(:)
      !! Array of evaluation energies of the S-matrix
    integer,                           intent(in)             :: J
      !! The current J
    type(asymtop_rot_channel_l_type),  intent(in)             :: channels_this_j(:)
      !! The array of channels for THIS J
    logical,                           intent(in)             :: is_first_chunk, is_last_chunk
      !! True if this is the first/last chunk to process; false otherwise

    logical :: allow_left_oob, allow_right_oob

    integer :: neleclo, nlo, kalo, kclo
    integer :: nelecup, nup, kaup, kcup
    integer :: nchans_J_flat, ne_chunk, ntrans
    integer :: ichan, itrans, ntrans_keep
    integer :: ie, ne_tot, nchans_j, nopen, nclosed, imap
    integer :: max_nup, max_nlo

    integer, allocatable :: idx_trans(:), idx_tmp(:)
    integer, allocatable :: idx_lo(:,:), idx_up(:,:)
    integer, allocatable :: n_idx_lo(:), n_idx_up(:)

    real(dp) :: etot, elo, eup
    real(dp) :: prob_term

    real(dp),    allocatable :: beta(:)

    complex(dp), allocatable :: s(:,:), s_flat(:)
    complex(dp), allocatable :: q(:)
    complex(dp), allocatable :: sphys(:,:)

    type(asymtop_rot_channel_type) :: lo, up

    ntrans        = size(transitions, 1)
    ne_tot        = size(total_energy_grid, 1)
    ne_chunk      = size(smat_rot_flat_chunk, 2)
    nchans_J      = size(channels_this_J, 1)
    nchans_J_flat = (nchans_J*(nchans_J+1))/2

    call size_check(smat_rot_flat_chunk,      [nchans_J_flat, ne_chunk], "SMAT_ROT_FLAT_CHUNK")
    call size_check(transition_probs,         [ntrans],                  "TRANSITION_PROBS")
    call size_check(smat_eval_energies_chunk, [ne_chunk],                "SMAT_EVAL_ENERGIES_CHUNK")

    ! -- get the transition indices that matter
    idx_trans = get_transition_map_this_J(transitions, channels_this_J)
    ntrans_keep = size(idx_trans, 1)
    if(ntrans_keep .eq. 0) then
      write(stderr, '("WARN: No transitions detected for J = ", I0, ". Exiting")') J
      return
    endif

    ! -- get the mapping of transitions that matter -> S^J
    max_nlo = 0
    max_nup = 0
    allocate(n_idx_lo(ntrans_keep))
    allocate(n_idx_up(ntrans_keep))
    do imap = 1, ntrans_keep
      itrans = idx_trans(imap)
      lo = transitions(itrans) % lo
      up = transitions(itrans) % up

      idx_tmp = pack([(ichan, ichan=1, nchans_J)], lo .eq. channels_this_J)
      n_idx_lo(imap) = size(idx_tmp, 1)
      max_nlo = max(max_nlo, n_idx_lo(imap))

      idx_tmp = pack([(ichan, ichan=1, nchans_J)], up .eq. channels_this_J)
      n_idx_up(imap) = size(idx_tmp, 1)
      max_nup = max(max_nup, n_idx_up(imap))
    enddo

    if(max_nlo .eq. 0 .OR. max_nup .eq. 0) then
      write(stderr, '("MAX_NLO: " ,I0)') max_nlo
      write(stderr, '("MAX_NUP: " ,I0)') max_nup
      call die("Neither max_nlo nor max_nup is allowed to be 0 at this point")
    endif

    allocate(idx_lo(max_nlo, ntrans_keep), source=0)
    allocate(idx_up(max_nup, ntrans_keep), source=0)

    do imap = 1, ntrans_keep
      itrans = idx_trans(imap)
      lo = transitions(itrans) % lo
      up = transitions(itrans) % up

      idx_tmp = pack([(ichan, ichan=1, nchans_J)], lo .eq. channels_this_J)
      idx_lo(1:n_idx_lo(imap), imap) = idx_tmp

      idx_tmp = pack([(ichan, ichan=1, nchans_J)], up .eq. channels_this_J)
      idx_up(1:n_idx_up(imap), imap) = idx_tmp
    enddo
    deallocate(idx_tmp)

    ! -- determine energy grid boundary conditions
    allow_left_oob  = is_first_chunk
    allow_right_oob = is_last_chunk

    ! -- loop over the total enrgy grid
    !$omp parallel default(none) &
    !$omp& shared(ne_tot, channels_this_J, transition_probs, transitions, J, nchans_J, G, total_energy_grid &
    !$omp&   , smat_eval_energies_chunk, smat_rot_flat_chunk, idx_lo, idx_up, n_idx_up, n_idx_lo, idx_trans &
    !$omp&   , allow_left_oob, allow_right_oob, ietot0, ietot1) &
    !$omp& private(ie, Etot, Sphys, beta, nopen, nclosed, lo, up, itrans&
    !$omp&   , q,  neleclo, nelecup, Nlo, Nup, Kalo, Kaup, Kclo, Kcup, Elo, Eup, prob_term &
    !$omp&   , S, S_flat)

    ! -- allocate S (for each thread) so that it can be used in the EDFT or EIFT case
    call realloc(S, nchans_J, nchans_J)
    call realloc(q, nchans_J)

    ! -- make nthreads copies of the energy independent S-matrix, or allocate nthreads
    !    vectors for each energy's flattened S-matrix
    if(G%EDFT .eqv. .false.) then
      call unpackmat_ch(smat_rot_flat_chunk(:,1), S, PACKMAT_TRIANGLE)
    else
      call realloc(s_flat, (nchans_J*(nchans_J+1))/2)
    endif

    !$omp do schedule(static)
    nrg: do ie=ietot0, ietot1

      Etot = total_energy_grid(ie)
      q    = fg_norm_coeff_q(channels_this_J, Etot)

      ! -- linear interpolation of S-matrix. OpenMP will use local copies of the matrices
      if(G%EDFT) then
        call interp_array_at_energy(Etot, smat_eval_energies_chunk, smat_rot_flat_chunk, S_flat, allow_left_oob, allow_right_oob)
        call unpackmat_ch(s_flat(:), s, PACKMAT_TRIANGLE)
      endif

      nclosed = count(channels_this_J % E .gt. Etot)
      nopen   = count(channels_this_J % E .le. Etot)

      if(nopen + nclosed .ne. nchans_J) call die("Number of opened and closed channels does not add to the number of channels !")

      call realloc(beta, nclosed)
      call realloc(Sphys, nopen, nopen)
      beta = 0
      Sphys = 0

      if(G%TARGCHARGE .eq. 0) then
        !TODO: CCEP from neutral BC ?
        Sphys = S(1:nopen,1:nopen)
      elseif(G%TARGCHARGE .gt. 0) then
        call CCEP(S, channels_this_J, q, Etot, Sphys, beta, nopen, nclosed)
      else
        call die("CCEP not implemented for negative ions")
      endif

      ! -- loop over the pairs of states for excitation, accumulate probabilities in each transition
      !    for this J
      call accumulate_probs( &
          transitions        &
        , idx_trans          &
        , idx_lo             &
        , n_idx_lo           &
        , idx_up             &
        , n_idx_up           &
        , ie                 &
        , Etot               &
        , J                  &
        , Sphys              &
        , transition_probs   &
      )

    enddo nrg
    !$omp end do
    !$omp end parallel

  end subroutine get_smat_probs_chunk

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  subroutine CCEP(S, channels, q, Etot, Sphys, beta, nopen, nclosed)
    !! Given the S-matrix, its basis of channels, the total energy E, and the number of open/closed channels,
    !! carry out the MQDT Closed-Channel Elimination Procedure to obtain the nopen x nopen physical S matrix

    use rotex__types,      only: asymtop_rot_channel_l_type
    use rotex__linalg,     only: zgesv, right_divide
    use rotex__arrays,     only: size_check
    use rotex__constants,  only: im, pi
    use rotex__characters, only: i2c => int2char

    implicit none (type, external)

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

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  pure function fg_norm_coeff_q(channels, Etot) result(res)
    !! Given an array of channels, determine the factor B for each channel at a particular total energy E
    use rotex__types,     only: asymtop_rot_channel_l_type
    use rotex__constants, only: pi
    implicit none (type, external)
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
    implicit none (type, external)
    real(dp), intent(in) :: e
    integer,  intent(in) :: l
    real(dp) :: res
    integer :: n
    res = product( [( 1+n*n*e, n=0, l )] )
  end function A_coulomb

  ! ---------------------------------------------------------------------------------------------------------------------------- !
  pure subroutine accumulate_probs( &
        transitions                 &
      , idx_trans                   &
      , idx_lo                      &
      , n_idx_lo                    &
      , idx_up                      &
      , n_idx_up                    &
      , ie                          &
      , Etot                        &
      , J                           &
      , Sphys                       &
      , transition_probs            &
    )
    !! Accumulate transition probabilities for the current energy and J

    use rotex__kinds, only: prob_rk
    use rotex__types, only: asymtop_rot_transition_type, prob_vector_type, asymtop_rot_channel_type

    implicit none (type, external)

    type(asymtop_rot_transition_type), intent(in)             :: transitions(:)
      !! Array of all transitions for the current spin
    integer,                           intent(in)             :: idx_trans(:)
      !! Array of indices that slice the transitions for the current J
    integer,                           intent(in)             :: idx_lo(:,:)
      !! Array map from transitions % lo -> S^J(ilo,:)
    integer,                           intent(in)             :: n_idx_lo(:)
      !! Number of indices for each transition
    integer,                           intent(in)             :: idx_up(:,:)
      !! Array map from transitions % up -> S^J(:,iup)
    integer,                           intent(in)             :: n_idx_up(:)
      !! Number of indices for each transition
    integer,                           intent(in)             :: ie
      !! The current energy index
    real(dp),                          intent(in)             :: Etot
      !! The current total energy
    integer,                           intent(in)             :: J
      !! The current J
    complex(dp),                       intent(in), contiguous :: Sphys(:,:)
      !! The current physical S-matrix
    type(prob_vector_type),            intent(inout)          :: transition_probs(:)
      !! Array of probabilities for each transition

    integer :: itrans, imap
    ! integer :: neleclo, Nlo, Kalo, Kclo
    ! integer :: nelecup, Nup, Kaup, Kcup
    real(dp) :: Elo, Eup
    real(dp) :: prob_term
    type(asymtop_rot_channel_type) :: lo, up

    do imap = 1, size(idx_trans, 1)
      itrans = idx_trans(imap)
      lo = transitions(itrans) % lo
      up = transitions(itrans) % up
      Elo = lo % E
      Eup = up % E

      ! -- skip energeticaly unavailable transitions
      if(Elo .ge. Etot) cycle
      if(Eup .gt. Etot) cycle

      ! -- Σ_{ll'} |Sphys_{il,i'l'}|²
      prob_term = (2*J+1) * sum(abs(            &
        Sphys( idx_lo(1:n_idx_lo(imap), imap)   &
             , idx_up(1:n_idx_up(imap), imap) ) &
      )**2)

      ! -- accumulate probabilities for this J
      transition_probs(itrans) % vec(ie) = transition_probs(itrans) % vec(ie) + real(prob_term, kind=prob_rk)

      if(prob_term .ge. 0) cycle

      ! -- error
      ! neleclo = lo % nelec
      ! Nlo     = lo % N
      ! Kalo    = lo % Ka
      ! Kclo    = lo % Kc
      ! nelecup = up % nelec
      ! Nup     = up % N
      ! Kaup    = up % Ka
      ! Kcup    = up % Kc
      ! write(stderr, '(A)') "❌"
      ! write(stderr, '(2X, A, 3E15.6)') "Elo, Eup, Etot: ", Elo, Eup, Etot
      ! write(stderr, '(2X, A, 4I4)') "Lower state nelec, N, Ka, Kc: ", neleclo, Nlo, Kalo, Kclo
      ! write(stderr, '(2X, A, 4I4)') "Upper state nelec, N, Ka, Kc: ", nelecup, Nup, Kaup, Kcup
      call die("Negative probability from the S-matrix !")

    enddo

  end subroutine accumulate_probs

  ! ------------------------------------------------------------------------------------------------------------------------------ !
  pure function get_transition_map_this_J(transitions, channels) result(idx)
    !! Determine which transitions are to be considered for the current J,
    !! given the current channels for this J
    use rotex__types,       only: asymtop_rot_transition_type, asymtop_rot_channel_l_type &
                                , asymtop_rot_channel_type
    use rotex__channel_ops, only: operator(.eq.), operator(.isin.)
    implicit none (type, external)
    type(asymtop_rot_transition_type), intent(in) :: transitions(:)
    type(asymtop_rot_channel_l_type),  intent(in) :: channels(:)
    logical, allocatable :: mask(:)
    integer, allocatable :: idx(:)
    integer :: itrans, ntrans
    type(asymtop_rot_channel_type) :: lo, up
    ntrans = size(transitions, 1)
    allocate(mask(ntrans), source=.false.)
    do itrans=1, ntrans
      lo = transitions(itrans) % lo
      up = transitions(itrans) % up
      ! -- skip elastic scattering
      if(lo .eq. up) cycle
      ! -- make sure both channels in this transition are actually in this J block
      if((lo .isin. channels) .eqv. .false.) cycle
      if((up .isin. channels) .eqv. .false.) cycle
      mask(itrans) = .true.
    enddo
    idx = pack([(itrans, itrans=1, ntrans)], mask)
  end function get_transition_map_this_J


! ================================================================================================================================ !
end module rotex__MQDTXS
! ================================================================================================================================ !
