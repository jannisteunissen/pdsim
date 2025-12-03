module m_photoi
  use m_random
  use m_config
  use m_units_constants
  use m_gas

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  real(dp)           :: pi_quench_fac
  real(dp)           :: pi_efficiency

  !> Minimum absorption length of photons
  real(dp)           :: pi_k1
  !> Maximum absorption length of photons
  real(dp)           :: pi_k2
  !> If > 0: length scale for non-ionizing absorption
  real(dp)           :: pi_k3
  logical, public    :: photoi_enabled = .true.

  ! Public methods
  public :: photoi_create_cfg
  public :: photoi_initialize
  public :: photoi_from_events
  public :: photoi_sample_photons

contains

  subroutine photoi_create_cfg(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "photoi%enabled", .true., &
         "Whether photoionization is used")
    call CFG_add(cfg, "photoi%model", "Zheleznyak", &
         "Photoionization model. Options are: " // &
         "Zheleznyak, Naidis_humid, Aints_humid, custom")
    call CFG_add(cfg, "photoi%absorp_inv_lengths", [0.0_dp, 0.0_dp], &
         "For custom model: inverse min/max absorption length. " // &
         "Will be scaled by gas pressure.")
    call CFG_add(cfg, "photoi%custom_k1_k2_k3", [0.0_dp, 0.0_dp, 0.0_dp], &
         "For custom model: coefficients k1, k2, k3. " // &
         "Will be scaled by gas pressure.")
    call CFG_add(cfg, "photoi%efficiency", 0.075_dp, &
         "Photoionization efficiency factor (scalar)")
  end subroutine photoi_create_cfg

  subroutine photoi_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg
    real(dp)                   :: temp_vec(3), p_O2, p_H2O
    character(len=40)          :: model
    real(dp), parameter        :: pq_air = 30.0_dp * UC_torr_to_bar
    real(dp), parameter        :: pq_H2O = 0.3_dp * UC_torr_to_bar


    call CFG_get(cfg, "photoi%enabled", photoi_enabled)

    p_O2 = GAS_get_fraction("O2") * GAS_pressure
    p_H2O = GAS_get_fraction("H2O") * GAS_pressure

    if (p_O2 <= epsilon(1.0_dp) .and. photoi_enabled) &
       error stop "There is no oxygen, you should disable photoionzation"

    call CFG_get(cfg, "photoi%model", model)

    select case (model)
    case ("Zheleznyak")
       pi_k1 = 3.5_dp / UC_torr_to_bar * p_O2
       pi_k2 = 200_dp / UC_torr_to_bar * p_O2
       pi_k3 = -1.0_dp
       pi_quench_fac = (pq_air) / (GAS_pressure + (pq_air))
    case ("Aints_humid")
       pi_k1 = 3.5_dp / UC_torr_to_bar * p_O2 + 0.13e2 / UC_torr_to_bar * p_H2O
       pi_k2 = 200_dp / UC_torr_to_bar * p_O2 + 0.57e2 / UC_torr_to_bar * p_H2O
       pi_k3 = -1.0_dp
       pi_quench_fac = 1/(1 + (GAS_pressure - p_H2O)/pq_air + p_H2O / pq_H2O)
    case ("Naidis_humid")
       pi_k1 = 3.5_dp / UC_torr_to_bar * p_O2
       pi_k2 = 200_dp / UC_torr_to_bar * p_O2
       pi_k3 = 0.26e2_dp / UC_torr_to_bar * p_H2O
       pi_quench_fac = (pq_air) / (GAS_pressure + (pq_air))
    case ("custom")
       call CFG_get(cfg, "photoi%custom_k1_k2_k3", temp_vec)
       if (minval(abs(temp_vec)) <= 0) &
            error stop "photoi%custom_k1_k2_k3 not set correctly"
       pi_k1 = temp_vec(1)
       pi_k2 = temp_vec(2)
       pi_k3 = temp_vec(3)
    case default
       error stop "Invalid photoi%model specified"
    end select

    call CFG_get(cfg, "photoi%efficiency", pi_efficiency)
    print *, pi_quench_fac
  end subroutine photoi_initialize

  subroutine photoi_from_events(n_events, events, rng, photo_pos, n_photons)
    use m_cross_sec
    use m_particle_core
    use m_units_constants

    integer, intent(in)          :: n_events !< Number of events
    type(PC_event_t), intent(in) :: events(n_events) !< list of events
    type(RNG_t), intent(inout)   :: rng
    real(dp), intent(inout)      :: photo_pos(:, :)
    integer, intent(out)         :: n_photons

    integer  :: i, n

    n_photons = 0
    do n = 1, n_events
       if (events(n)%ctype == CS_ionize_t) then
          call photoi_sample_photons(rng, events(n)%part%x, events(n)%part%w, &
               huge(1), i, photo_pos(:, n_photons+1:))
          n_photons = n_photons + i
       end if
    end do
  end subroutine photoi_from_events

  subroutine photoi_sample_photons(rng, x_source, num_ionizations, &
       max_photons, n_photons, absorption_locations)
    type(RNG_t), intent(inout) :: rng
    real(dp), intent(in)       :: x_source(3)
    real(dp), intent(in)       :: num_ionizations
    integer, intent(in)        :: max_photons !< Maximum number of photons
    integer, intent(out)       :: n_photons
    real(dp), intent(inout)    :: absorption_locations(:, :)

    integer  :: n, n_photons_max
    real(dp) :: num_mean, en_frac, r, lambda, p_no_ionization, u01

    num_mean  = pi_efficiency * pi_quench_fac * num_ionizations
    n_photons_max = min(rng%poisson(num_mean), max_photons)

    if (size(absorption_locations, 1) /= 3) &
         error stop "absorption_locations should be sized (3, max_photons)"

    n_photons = 0

    ! If there are too many photons, do not store all of them
    do n = 1, min(n_photons_max, size(absorption_locations, 2))
       en_frac = rng%unif_01()
       lambda = get_photoi_lambda(en_frac)

       ! Sample from exponential distribution
       u01 = rng%unif_01()
       r = -rng%log1p(-u01) / lambda

       ! Check if photon was absorbed without ionization
       if (pi_k3 > 0) then
          p_no_ionization = exp(-pi_k3 * r)
          if (rng%unif_01() > p_no_ionization) cycle
       end if

       n_photons = n_photons + 1
       absorption_locations(:, n_photons) = x_source + rng%sphere(r)
    end do

  end subroutine photoi_sample_photons

  ! Returns the inverse mean free path for a photon.
  real(dp) function get_photoi_lambda(en_frac)
    real(dp), intent(in) :: en_frac
    get_photoi_lambda = pi_k1 * (pi_k2/pi_k1)**en_frac
  end function get_photoi_lambda

end module m_photoi
