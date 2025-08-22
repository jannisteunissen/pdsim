!> Module for computing ionization integral
module m_integral
  use m_config
  use m_lookup_table
  use m_interp_unstructured
  use m_pdsim

  implicit none
  private

  ! Index for storing travel time
  integer, parameter :: i_travel_time = 1
  ! Index for storing integral of alpha - eta
  integer, parameter :: i_Kint = 2
  ! Index for storing integral of alpha + eta
  integer, parameter :: i_Lint = 3

  ! Public methods
  public :: integral_create_config
  public :: integral_compute

contains

  subroutine integral_create_config(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "integral%max_steps", 1000, &
         "Maximum number of steps for K integral")
    call CFG_add(cfg, "integral%rtol", 1e-5_dp, &
         "Relative tolerance per step for K integral")
    call CFG_add(cfg, "integral%atol", 1e-5_dp, &
         "Absolute tolerance per step for K integral")
    call CFG_add(cfg, "integral%min_dx", 1e-7_dp, &
         "Minimum step size for K integral")
    call CFG_add(cfg, "integral%max_dx", 1e-3_dp, &
         "Maximum step size for K integral")
    call CFG_add(cfg, "integral%boundary_distance", 1e-6_dp, &
         "For integration, keep this distance away from boundaries")
    call CFG_add(cfg, "integral%move_distance", 1e-9_dp, &
         "Move points this distance to avoid problems at interfaces")
  end subroutine integral_create_config

  subroutine integral_compute(cfg)
    type(cfg_t), intent(inout) :: cfg
    integer                    :: n
    integer                    :: max_steps, n_steps
    integer                    :: boundary_material
    real(dp)                   :: min_dx, max_dx, boundary_distance
    real(dp)                   :: move_distance
    real(dp)                   :: rtol, atol, r(3), w
    real(dp)                   :: field(pdsim_ndim), td(pdsim_ncols)
    real(dp)                   :: r_final(3), p_m1, sum_ioniz, t_avg
    real(dp), allocatable      :: y(:, :), y_field(:, :)
    real(dp), allocatable      :: moved_points(:, :)
    logical                    :: reverse
    integer, parameter         :: nvar = 3

    call CFG_get(cfg, "integral%max_steps", max_steps)
    call CFG_get(cfg, "integral%rtol", rtol)
    call CFG_get(cfg, "integral%atol", atol)
    call CFG_get(cfg, "integral%min_dx", min_dx)
    call CFG_get(cfg, "integral%max_dx", max_dx)
    call CFG_get(cfg, "integral%boundary_distance", boundary_distance)
    call CFG_get(cfg, "integral%move_distance", move_distance)

    call iu_add_point_data(pdsim_ug, "K_integral", i_k_integral)
    call iu_add_point_data(pdsim_ug, "alpha", i_alpha)
    call iu_add_point_data(pdsim_ug, "eta", i_eta)
    call iu_add_point_data(pdsim_ug, "alpha_eff", i_alpha_eff)
    call iu_add_point_data(pdsim_ug, "avalanche_p0", i_p0)
    call iu_add_point_data(pdsim_ug, "avalanche_w", i_w)
    call iu_add_point_data(pdsim_ug, "avalanche_time", i_avalanche_time)
    call iu_add_point_data(pdsim_ug, "avalanche_x1", i_x1)
    call iu_add_point_data(pdsim_ug, "avalanche_x2", i_x2)
    call iu_add_point_data(pdsim_ug, "avalanche_x3", i_x3)

    call iu_add_point_data(pdsim_ug, "avalanche_p_m1", i_p_m1)
    call iu_add_point_data(pdsim_ug, "K_star", i_kstar)

    call iu_add_point_data(pdsim_ug, "ion_time", i_ion_time)
    call iu_add_point_data(pdsim_ug, "ion_x1", i_ion_x1)
    call iu_add_point_data(pdsim_ug, "ion_x2", i_ion_x2)
    call iu_add_point_data(pdsim_ug, "ion_x3", i_ion_x3)
    call iu_add_point_data(pdsim_ug, "ion_gamma", i_ion_gamma)

    ! Move mesh points slightly away from domain boundaries
    call move_mesh_points_slightly(pdsim_ug, boundary_distance, &
         move_distance, moved_points)

    allocate(y(pdsim_ndim+nvar, max_steps))
    allocate(y_field(pdsim_ndim, max_steps))

    !$omp parallel do private(r, y, y_field, n_steps, field, td, &
    !$omp w, reverse, sum_ioniz, p_m1, r_final, t_avg, boundary_material)
    do n = 1, pdsim_ug%n_points
       r = moved_points(:, n)

       ! Trace electron avalanche
       reverse = .true.
       call iu_integrate_along_field(pdsim_ug, pdsim_ndim, electron_sub, &
            r(1:pdsim_ndim), pdsim_pdata_field(1:pdsim_ndim), &
            min_dx, max_dx, max_steps, rtol, atol, reverse, nvar, y, y_field, &
            n_steps, pdsim_axisymmetric, &
            pdsim_icdata_material, pdsim_gas_material_value)

       if (n_steps > max_steps) then
          print *, "Error for electron start position", r(1:pdsim_ndim)
          print *, "n_steps = ", n_steps, ", max_steps = ", max_steps
          error stop "Increase integral%max_steps"
       end if

       ! Compute integrals along path
       call compute_path_integrals(pdsim_ndim, nvar, n_steps, y(:, 1:n_steps), &
            y_field(:, 1:n_steps), w, sum_ioniz, p_m1, &
            r_final(1:pdsim_ndim), t_avg)

       pdsim_ug%point_data(n, i_k_integral) = maxval(y(pdsim_ndim+i_Kint, 1:n_steps))
       pdsim_ug%point_data(n, i_avalanche_time) = t_avg

       ! Store position
       r_final(pdsim_ndim+1:) = 0.0_dp
       pdsim_ug%point_data(n, [i_x1, i_x2, i_x3]) = r_final

       ! Store alpha and eta
       field = pdsim_ug%point_data(n, pdsim_pdata_field(1:pdsim_ndim))
       td = LT_get_mcol(pdsim_tdtbl, norm2(field))
       pdsim_ug%point_data(n, i_alpha) = td(pdsim_col_alpha)
       pdsim_ug%point_data(n, i_eta) = td(pdsim_col_eta)
       pdsim_ug%point_data(n, i_alpha_eff) = td(pdsim_col_alpha) - td(pdsim_col_eta)

       pdsim_ug%point_data(n, i_w) = w
       pdsim_ug%point_data(n, i_p0) = 1 - exp(y(pdsim_ndim+i_Kint, n_steps))/w
       pdsim_ug%point_data(n, i_p_m1) = p_m1
       pdsim_ug%point_data(n, i_kstar) = log(1 + sum_ioniz)

       ! Trace positive ions
       reverse = .false.
       call iu_integrate_along_field(pdsim_ug, pdsim_ndim, ion_sub, &
            r(1:pdsim_ndim), pdsim_pdata_field(1:pdsim_ndim), &
            min_dx, max_dx, max_steps, rtol, atol, reverse, nvar, y, y_field, &
            n_steps, pdsim_axisymmetric, &
            pdsim_icdata_material, pdsim_gas_material_value, boundary_material)

       if (n_steps > max_steps) then
          print *, "Error for ion start position", r(1:pdsim_ndim)
          print *, "n_steps = ", n_steps, ", max_steps = ", max_steps
          error stop "Increase integral%max_steps"
       end if

       pdsim_ug%point_data(n, i_ion_time) = y(pdsim_ndim+i_travel_time, n_steps)

       ! Store final position
       r_final(pdsim_ndim+1:) = 0.0_dp
       r_final(1:pdsim_ndim) = y(1:pdsim_ndim, n_steps)
       pdsim_ug%point_data(n, [i_ion_x1, i_ion_x2, i_ion_x3]) = r_final

       ! Store secondary emission coefficient
       if (boundary_material < 0) then
          ! Domain boundary
          pdsim_ug%point_data(n, i_ion_gamma) = pdsim_ion_gamma(1)
       else
          ! Set coefficient for this material
          pdsim_ug%point_data(n, i_ion_gamma) = &
               pdsim_ion_gamma(1+boundary_material)
       end if

    end do
    !$omp end parallel do

  end subroutine integral_compute

  !> Move mesh points slightly away from domain boundaries. Other points are
  !> also slightly moved, to avoid them being exactly on a material boundary
  subroutine move_mesh_points_slightly(ug, bdist, mdist, points)
    type(iu_grid_t), intent(in)          :: ug
    !> How far points are moved from domain boundaries
    real(dp), intent(in)                 :: bdist
    !> How far other points are moved
    real(dp), intent(in)                 :: mdist
    real(dp), allocatable, intent(inout) :: points(:, :)

    integer              :: n, k, i_point
    real(dp)             :: center(3), direction(3), dist
    logical, allocatable :: moved(:)

    allocate(points(3, ug%n_points))
    allocate(moved(ug%n_points))

    points(:, :) = ug%points
    moved(:) = .false.

    do n = 1, ug%n_cells
       ! For cells in the gas, move boundary points inwards
       if (ug%icell_data(n, pdsim_icdata_material) == &
            pdsim_gas_material_value) then
          center = iu_get_cell_center(ug, n)
          do k = 1, ug%n_points_per_cell
             i_point = ug%cells(k, n)

             ! Move points only once
             if (.not. moved(i_point)) then
                if (ug%point_is_at_boundary(i_point)) then
                   dist = bdist
                else
                   dist = mdist
                end if

                direction = center - points(:, i_point)
                points(:, i_point) = points(:, i_point) + dist * &
                     direction / norm2(direction)
                moved(i_point) = .true.
             end if
          end do
       end if
    end do
  end subroutine move_mesh_points_slightly

  !> Helper routine for integration along field lines
  subroutine electron_sub(ndim, r, field, nvar, integrand)
    integer, intent(in)   :: ndim
    real(dp), intent(in)  :: r(ndim)
    real(dp), intent(in)  :: field(ndim)
    integer, intent(in)   :: nvar
    real(dp), intent(out) :: integrand(nvar)
    real(dp)              :: field_norm, td(pdsim_ncols)

    field_norm = norm2(field)
    td = LT_get_mcol(pdsim_tdtbl, field_norm)

    ! 1/velocity
    integrand(i_travel_time) = 1/(field_norm * td(pdsim_col_mu))
    ! alpha - eta
    integrand(i_Kint) = td(pdsim_col_alpha) - td(pdsim_col_eta)
    ! alpha + eta
    integrand(i_Lint) = td(pdsim_col_alpha) + td(pdsim_col_eta)

    integrand(4:nvar) = 0.0_dp
  end subroutine electron_sub

  !> Helper routine for integration along field lines
  subroutine ion_sub(ndim, r, field, nvar, integrand)
    integer, intent(in)   :: ndim
    real(dp), intent(in)  :: r(ndim)
    real(dp), intent(in)  :: field(ndim)
    integer, intent(in)   :: nvar
    real(dp), intent(out) :: integrand(nvar)
    real(dp)              :: field_norm

    field_norm = norm2(field)
    ! 1/velocity
    integrand(i_travel_time) = 1/(field_norm * pdsim_ion_mobility)
    integrand(2:nvar) = 0.0_dp
  end subroutine ion_sub

  !> Compute several integrals along the electric field streamline
  subroutine compute_path_integrals(ndim, nvar, n_steps, y, y_field, w, &
       sum_ioniz, p_m1, r_avg, t_avg)
    integer, intent(in)   :: ndim
    integer, intent(in)   :: nvar
    integer, intent(in)   :: n_steps
    real(dp), intent(in)  :: y(ndim+nvar, n_steps)
    real(dp), intent(in)  :: y_field(ndim, n_steps)
    !> The w value at the final point, as defined in Kendall's 1948 paper:
    !> w = 1 + exp(-rho(x)) * integral(exp(rho(x')) * alpha(x') dx')
    !> for x' between 0 and x. Here rho(x) = -K(x).
    real(dp), intent(out) :: w
    !> The total number of ionizations given by integral(exp(-rho(x')) *
    !> alpha(x') dx')
    real(dp), intent(out) :: sum_ioniz
    !> The probability that there are no ionizations
    real(dp), intent(out) :: p_m1
    !> The 'average' coordinate at which ionization is produced
    real(dp), intent(out) :: r_avg(ndim)
    !> The 'average' time at which ionization is produced
    real(dp), intent(out) :: t_avg

    integer  :: n
    real(dp) :: alpha, eta, dx, field_norm, tmp, k0, k1
    real(dp) :: f1(n_steps), f2(n_steps), f3(n_steps)

    ! Threshold for non-zero ionization
    real(dp), parameter :: ionization_threshold = 1e-3_dp

    ! Determine integrands at every point
    do n = 1, n_steps
       field_norm = norm2(y_field(:, n))
       alpha = LT_get_col(pdsim_tdtbl, pdsim_col_alpha, field_norm)
       eta = LT_get_col(pdsim_tdtbl, pdsim_col_eta, field_norm)
       f1(n) = exp(-y(ndim+i_Kint, n)) * alpha
       f2(n) = exp(y(ndim+i_Kint, n)) * alpha
       f3(n) = eta
    end do

    ! Use composite trapezoidal rule to evaluate integrals
    w = 0.0_dp
    sum_ioniz = 0.0_dp
    p_m1 = 0.0_dp
    r_avg = 0.0_dp
    t_avg = 0.0_dp

    do n = 1, n_steps-1
       dx = norm2(y(1:ndim, n+1) - y(1:ndim, n))
       w = w + dx * 0.5_dp * (f1(n) + f1(n+1))

       tmp = dx * 0.5_dp * (f2(n) + f2(n+1))
       sum_ioniz = sum_ioniz + tmp

       ! Average coordinate where ionization is produced
       r_avg = r_avg + tmp * 0.5_dp * (y(1:ndim, n) + y(1:ndim, n+1))
       t_avg = t_avg + tmp * 0.5_dp * (y(ndim+i_travel_time, n) + &
            y(ndim+i_travel_time, n+1))

       ! Assume linear variation of term inside exponential, of the form k0 +
       ! k1 * (x - x0). Then integrate exponential term analytically,
       ! since it rapidly decays.
       k0 = y(ndim+i_Lint, n)
       k1 = (y(ndim+i_Lint, n+1) - y(ndim+i_Lint, n))/dx
       p_m1 = p_m1 + 0.5_dp * (f3(n) + f3(n+1)) * exp(-k0) * (1 - exp(-k1*dx))/k1
    end do

    ! 'Average' position and time where ionization was produced
    if (sum_ioniz > ionization_threshold) then
       r_avg = r_avg / sum_ioniz
       t_avg = t_avg / sum_ioniz
    else
       ! No ionization produced; take start location and time
       r_avg = y(1:ndim, 1)
       t_avg = y(ndim+i_travel_time, 1)
    end if

    w = 1 + exp(y(ndim+i_Kint, n_steps)) * w
    p_m1 = exp(-y(ndim+i_Lint, n_steps)) + p_m1
  end subroutine compute_path_integrals

end module m_integral
