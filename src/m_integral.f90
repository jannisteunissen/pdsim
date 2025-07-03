!> Module for computing ionization integral
module m_integral
  use m_config
  use m_lookup_table
  use m_interp_unstructured
  use m_pdsim

  implicit none
  private

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
         "Initially, keep this distance away from boundaries")
  end subroutine integral_create_config

  subroutine integral_compute(cfg)
    type(cfg_t), intent(inout) :: cfg
    integer                    :: n
    integer                    :: max_steps, n_steps
    real(dp)                   :: min_dx, max_dx, boundary_distance
    real(dp)                   :: rtol, atol, r(3), w
    real(dp)                   :: field(pdsim_ndim), td(pdsim_ncols)
    real(dp)                   :: r_final(3)
    real(dp), allocatable      :: y(:, :), y_field(:, :)
    real(dp), allocatable      :: moved_points(:, :)
    logical                    :: reverse
    integer, parameter         :: nvar = 2

    call CFG_get(cfg, "integral%max_steps", max_steps)
    call CFG_get(cfg, "integral%rtol", rtol)
    call CFG_get(cfg, "integral%atol", atol)
    call CFG_get(cfg, "integral%min_dx", min_dx)
    call CFG_get(cfg, "integral%max_dx", max_dx)
    call CFG_get(cfg, "integral%boundary_distance", boundary_distance)

    call iu_add_point_data(pdsim_ug, "K_integral", i_k_integral)
    call iu_add_point_data(pdsim_ug, "alpha", i_alpha)
    call iu_add_point_data(pdsim_ug, "eta", i_eta)
    call iu_add_point_data(pdsim_ug, "avalanche_p0", i_p0)
    call iu_add_point_data(pdsim_ug, "avalanche_w", i_w)
    call iu_add_point_data(pdsim_ug, "avalanche_time", i_avalanche_time)
    call iu_add_point_data(pdsim_ug, "avalanche_x1", i_x1)
    call iu_add_point_data(pdsim_ug, "avalanche_x2", i_x2)
    call iu_add_point_data(pdsim_ug, "avalanche_x3", i_x3)

    call iu_add_point_data(pdsim_ug, "ion_time", i_ion_time)
    call iu_add_point_data(pdsim_ug, "ion_x1", i_ion_x1)
    call iu_add_point_data(pdsim_ug, "ion_x2", i_ion_x2)
    call iu_add_point_data(pdsim_ug, "ion_x3", i_ion_x3)

    ! Move mesh points slightly away from domain boundaries
    call move_mesh_points_from_boundary(pdsim_ug, boundary_distance, &
         moved_points)

    allocate(y(pdsim_ndim+nvar, max_steps))
    allocate(y_field(pdsim_ndim, max_steps))

    !$omp parallel do private(r, y, y_field, n_steps, field, td, &
    !$omp w, reverse, r_final)
    do n = 1, pdsim_ug%n_points
       r = moved_points(:, n)

       ! Trace electron avalanche
       reverse = .true.
       call iu_integrate_along_field(pdsim_ug, pdsim_ndim, electron_sub, &
            r(1:pdsim_ndim), pdsim_pdata_field(1:pdsim_ndim), &
            min_dx, max_dx, max_steps, rtol, atol, reverse, nvar, y, y_field, &
            n_steps, pdsim_axisymmetric, &
            pdsim_cdata_material, pdsim_gas_material_value)

       if (n_steps > max_steps) then
          print *, "Error for electron start position", r(1:pdsim_ndim)
          print *, "n_steps = ", n_steps, ", max_steps = ", max_steps
          error stop "Increase integral%max_steps"
       end if

       pdsim_ug%point_data(n, i_k_integral) = y(pdsim_ndim+1, n_steps)
       pdsim_ug%point_data(n, i_avalanche_time) = y(pdsim_ndim+2, n_steps)

       ! Store final position
       r_final(pdsim_ndim+1:) = 0.0_dp
       r_final(1:pdsim_ndim) = y(1:pdsim_ndim, n_steps)
       pdsim_ug%point_data(n, [i_x1, i_x2, i_x3]) = r_final

       ! Store alpha and eta
       field = pdsim_ug%point_data(n, pdsim_pdata_field(1:pdsim_ndim))
       td = LT_get_mcol(pdsim_tdtbl, norm2(field))
       pdsim_ug%point_data(n, i_alpha) = td(pdsim_col_alpha)
       pdsim_ug%point_data(n, i_eta) = td(pdsim_col_eta)

       ! Compute w according to Kendall's 1948 paper
       call compute_kendall_w(pdsim_ndim, nvar, n_steps, y(:, 1:n_steps), &
            y_field(:, 1:n_steps), w)
       pdsim_ug%point_data(n, i_w) = w
       pdsim_ug%point_data(n, i_p0) = 1 - exp(y(pdsim_ndim+1, n_steps))/w

       ! Trace positive ions
       reverse = .false.
       call iu_integrate_along_field(pdsim_ug, pdsim_ndim, ion_sub, &
            r(1:pdsim_ndim), pdsim_pdata_field(1:pdsim_ndim), &
            min_dx, max_dx, max_steps, rtol, atol, reverse, nvar, y, y_field, &
            n_steps, pdsim_axisymmetric, &
            pdsim_cdata_material, pdsim_gas_material_value)

       if (n_steps > max_steps) then
          print *, "Error for ion start position", r(1:pdsim_ndim)
          print *, "n_steps = ", n_steps, ", max_steps = ", max_steps
          error stop "Increase integral%max_steps"
       end if

       pdsim_ug%point_data(n, i_ion_time) = y(pdsim_ndim+1, n_steps)

       ! Store final position
       r_final(pdsim_ndim+1:) = 0.0_dp
       r_final(1:pdsim_ndim) = y(1:pdsim_ndim, n_steps)
       pdsim_ug%point_data(n, [i_ion_x1, i_ion_x2, i_ion_x3]) = r_final
    end do
    !$omp end parallel do

  end subroutine integral_compute

  !> Move mesh points slightly away from domain boundaries
  subroutine move_mesh_points_from_boundary(ug, distance, points)
    type(iu_grid_t), intent(in)          :: ug
    real(dp), intent(in)                 :: distance
    real(dp), allocatable, intent(inout) :: points(:, :)

    integer              :: n, k, i_point
    real(dp)             :: center(3), direction(3)
    logical, allocatable :: moved(:)

    allocate(points(3, ug%n_points))
    allocate(moved(ug%n_points))

    points = ug%points
    moved = .false.

    do n = 1, ug%n_cells
       ! For cells in the gas, move boundary points inwards
       if (abs(ug%cell_data(n, pdsim_cdata_material) - &
            pdsim_gas_material_value) <= 0) then
          do k = 1, ug%n_points_per_cell
             i_point = ug%cells(k, n)

             ! Move points only once
             if (ug%point_is_at_boundary(i_point) .and. &
                  .not. moved(i_point)) then
                center = iu_get_cell_center(ug, n)
                direction = center - points(:, i_point)
                points(:, i_point) = points(:, i_point) + distance * &
                     direction / norm2(direction)
                moved(i_point) = .true.
             end if
          end do
       end if
    end do
  end subroutine move_mesh_points_from_boundary

  !> Computes alpha effective and 1/velocity for electrons
  subroutine electron_sub(ndim, r, field, nvar, integrand)
    integer, intent(in)   :: ndim
    real(dp), intent(in)  :: r(ndim)
    real(dp), intent(in)  :: field(ndim)
    integer, intent(in)   :: nvar
    real(dp), intent(out) :: integrand(nvar)
    real(dp)              :: field_norm, td(pdsim_ncols)

    field_norm = norm2(field)
    td = LT_get_mcol(pdsim_tdtbl, field_norm)
    integrand(1) = td(pdsim_col_alpha) - td(pdsim_col_eta)
    integrand(2) = 1/(field_norm * td(pdsim_col_mu))
    integrand(3:nvar) = 0.0_dp
  end subroutine electron_sub

  !> Computes 1/velocity for ions
  subroutine ion_sub(ndim, r, field, nvar, integrand)
    integer, intent(in)   :: ndim
    real(dp), intent(in)  :: r(ndim)
    real(dp), intent(in)  :: field(ndim)
    integer, intent(in)   :: nvar
    real(dp), intent(out) :: integrand(nvar)
    real(dp)              :: field_norm

    field_norm = norm2(field)
    integrand(1) = 1/(field_norm * pdsim_ion_mobility)
    integrand(2:nvar) = 0.0_dp
  end subroutine ion_sub

  !> Compute the w value at the final point, as defined in Kendall's 1948
  !> paper. This can be written as:
  !> w = 1 + exp(-rho(x)) * integral(exp(rho(x')) * alpha(x') dx')
  !> for x' between 0 and x. Here rho(x) = -K(x).
  subroutine compute_kendall_w(ndim, nvar, n_steps, y, y_field, w)
    integer, intent(in)   :: ndim
    integer, intent(in)   :: nvar
    integer, intent(in)   :: n_steps
    real(dp), intent(in)  :: y(ndim+nvar, n_steps)
    real(dp), intent(in)  :: y_field(ndim, n_steps)
    real(dp), intent(out) :: w

    integer               :: n
    real(dp)              :: alpha, dx, rho
    real(dp)              :: integrand(n_steps)

    ! Determine integrand at every point
    do n = 1, n_steps
       alpha = LT_get_col(pdsim_tdtbl, pdsim_col_alpha, norm2(y_field(:, n)))
       rho = -y(ndim+1, n)
       integrand(n) = exp(rho) * alpha
    end do

    ! Use composite trapezoidal rule to evaluate integral
    w = 0.0_dp
    do n = 1, n_steps-1
       dx = norm2(y(1:ndim, n+1) - y(1:ndim, n))
       w = w + dx * 0.5_dp * (integrand(n) + integrand(n+1))
    end do

    w = 1 + exp(y(ndim+1, n_steps)) * w
  end subroutine compute_kendall_w

end module m_integral
