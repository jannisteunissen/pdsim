module m_photoi
  use m_random
  use m_config
  use m_units_constants
  use m_gas

  implicit none
  private

  integer, parameter    :: dp = kind(0.0d0)
  real(dp)              :: pi_quench_fac
  real(dp)              :: pi_efficiency
  real(dp)              :: pi_min_inv_abs_len, pi_max_inv_abs_len

  logical, public :: photoi_enabled = .true.

  ! Public methods
  public :: photoi_create_cfg
  public :: photoi_initialize
  public :: photoi_from_events

contains

  subroutine photoi_create_cfg(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "photoi%enabled", .true., &
         "Whether photoionization is used")
    call CFG_add(cfg, "photoi%absorp_inv_lengths", &
         [3.5D0 / UC_torr_to_bar, 200D0 / UC_torr_to_bar], &
         "The inverse min/max absorption length, will be scaled by pO2")
    call CFG_add(cfg, "photoi%efficiency", 0.075_dp, &
         "Photoionization efficiency factor (instead of table)")
  end subroutine photoi_create_cfg

  subroutine photoi_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg
    real(dp)                   :: frac_O2, temp_vec(2)

    call CFG_get(cfg, "photoi%enabled", photoi_enabled)

    if (photoi_enabled) then
       frac_O2 = GAS_get_fraction("O2")
       if (frac_O2 <= epsilon(1.0_dp)) then
          error stop "There is no oxygen, you should disable photoionzation"
       end if

       call CFG_get(cfg, "photoi%absorp_inv_lengths", temp_vec)
       pi_min_inv_abs_len = temp_vec(1) * frac_O2 * GAS_pressure
       pi_max_inv_abs_len = temp_vec(2) * frac_O2 * GAS_pressure

       call CFG_get(cfg, "photoi%efficiency", pi_efficiency)

       pi_quench_fac = (30.0D0 * UC_torr_to_bar) / &
            (GAS_pressure + (30.0D0 * UC_torr_to_bar))
    end if
  end subroutine photoi_initialize

  subroutine photoi_from_events(n_events, events, rng, photo_pos, photo_w, n_photons)
    use m_cross_sec
    use m_particle_core
    use m_units_constants

    integer, intent(in)          :: n_events !< Number of events
    type(PC_event_t), intent(in) :: events(n_events) !< list of events
    type(RNG_t), intent(inout)   :: rng
    real(dp), intent(inout)      :: photo_pos(:, :)
    real(dp), intent(inout)      :: photo_w(:)
    integer, intent(out)         :: n_photons

    integer  :: i, n

    n_photons = 0
    do n = 1, n_events
       if (events(n)%ctype == CS_ionize_t) then
          call photoi_sample_photons(rng, events(n)%part%x, events(n)%part%w, &
               i, photo_pos(:, n_photons+1:))
          n_photons = n_photons + i
       end if
    end do
  end subroutine photoi_from_events

  subroutine photoi_sample_photons(rng, x_source, num_ionizations, &
       n_photons, absorption_locations)
    type(RNG_t), intent(inout) :: rng
    real(dp), intent(in)       :: x_source(3)
    real(dp), intent(in)       :: num_ionizations
    integer, intent(out)       :: n_photons
    real(dp), intent(inout)    :: absorption_locations(:, :)

    integer  :: n
    real(dp) :: num_mean, en_frac, fly_len

    num_mean  = pi_efficiency * pi_quench_fac * num_ionizations
    n_photons = rng%poisson(num_mean)

    if (size(absorption_locations, 1) /= 3) &
         error stop "absorption_locations should be sized (3, max_photons)"
    if (n_photons > size(absorption_locations, 2)) &
         error stop "Too many photons were generated"

    do n = 1, n_photons
       en_frac = rng%unif_01()
       fly_len = -log(1.0_dp - rng%unif_01()) / get_photoi_lambda(en_frac)
       absorption_locations(:, n) = x_source + rng%sphere(fly_len)
    end do

  end subroutine photoi_sample_photons

  ! Returns the inverse mean free path for a photon.
  real(dp) function get_photoi_lambda(en_frac)
    real(dp), intent(in) :: en_frac
    get_photoi_lambda = pi_min_inv_abs_len * &
         (pi_max_inv_abs_len/pi_min_inv_abs_len)**en_frac
  end function get_photoi_lambda

end module m_photoi
