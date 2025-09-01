!> Module with settings and routines for tabulated data
module m_table_data
  use m_lookup_table

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! The maximum number of rows per entry
  integer, parameter :: table_max_rows = 1500

  ! Maximum length of a line
  integer, parameter :: string_len = 200

  ! Interpolation methods
  integer, parameter, public :: table_interp_linear = 1
  integer, parameter, public :: table_interp_cubic_spline = 2

  ! Public methods
  public :: table_from_file
  public :: table_strings_numbers_from_file
  public :: table_set_column

contains

  !> Interpolate data and store in lookup table
  subroutine table_set_column(tbl, i_col, x, y, interpolation, max_err, add)
    use m_spline_interp
    use m_lookup_table
    type(LT_t), intent(inout)       :: tbl     !< Lookup table
    integer, intent(in)             :: i_col   !< Index of column
    real(dp), intent(in)            :: x(:)
    real(dp), intent(in)            :: y(:)
    integer, intent(in), optional   :: interpolation
    real(dp), intent(out), optional :: max_err !< Estimate of maximal error
    logical, intent(in), optional   :: add     !< Whether to add to existing data
    real(dp), allocatable           :: y_table(:)
    type(spline_t)                  :: spl
    integer                         :: input_interpolation

    if (size(x) /= size(y)) error stop "size(x) /= size(y)"

    input_interpolation = table_interp_linear
    if (present(interpolation)) input_interpolation = interpolation

    select case (input_interpolation)
    case (table_interp_linear)
       call LT_set_col(tbl, i_col, x, y, add)
    case (table_interp_cubic_spline)
       ! Perform cubic spline interpolation
       call spline_set_coeffs(x, y, size(x), spl)
       y_table = spline_evaluate(tbl%x, spl)

       if (minval(y) >= 0.0_dp) then
          ! If original data is non-negative, ensure interpolated data is also
          ! non-negative (important for e.g. rate coefficients)
          y_table = max(0.0_dp, y_table)
       end if

       call LT_set_col_data(tbl, i_col, y_table, add)
    case default
       error stop "invalid interpolation"
    end select

    if (present(max_err)) then
       ! Largest difference divided by maximum value of |y|
       max_err = maxval(abs(y - LT_get_col(tbl, i_col, x)))/maxval(abs(y))
    end if
  end subroutine table_set_column

  !> Routine to read in tabulated data from a file
  subroutine table_from_file(file_name, data_name, x_data, y_data)
    character(len=*), intent(in)       :: file_name, data_name
    real(dp), allocatable, intent(out) :: x_data(:), y_data(:)

    ! Temporary variables
    integer                   :: ioState, nL
    integer                   :: n_rows
    integer                   :: my_unit
    character(len=40)         :: line_fmt
    character(len=string_len) :: line
    real(dp)                  :: temp_table(2, table_max_rows)
    real(dp)                  :: factor

    nL = 0 ! Set the number of lines to 0

    ! Set the line format to read, only depends on string_len currently
    write(line_fmt, FMT = "(I6)") string_len
    line_fmt = "(A" // trim(adjustl(line_fmt)) // ")"

    ! Open 'file_name' (with error checking)
    open(newunit=my_unit, file = trim(file_name), action = "read", &
         err = 999, iostat = ioState, status="old")

    ! Table format

    !     table_name
    !     FACTOR: 1.0                   [optional: multiply with this factor]
    !     [other lines]
    !     ------------------            [at least 5 dashes]
    !     xxx       xxx                 [data in two column format]
    !     ...       ...
    !     xxx       xxx
    !     ------------------

    ! The outer DO loop, running until the end of the file is reached
    do
       ! Search for 'data_name' in the file
       do
          read(my_unit, FMT = line_fmt, ERR = 999, end = 888) line; nL = nL+1
          if (line == data_name) exit
       end do

       factor = 1.0_dp

       ! Now we can check whether there is a comment, while scanning lines until
       ! dashes are found, which indicate the start of the data
       do
          read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
          line = adjustl(line)
          if ( line(1:5) == "-----" ) then
             exit
          else if (line(1:7) == "FACTOR:") then
             read(line(8:), *) factor
          else if (line(1:8) == "COMMENT:") then
             continue
          else
             print *, "In file ", trim(file_name), " at line", nL
             print *, trim(line)
             error stop "Unknown statement in input file"
          end if
       end do

       ! Read the data into a temporary array
       n_rows = 0
       do
          read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
          line = adjustl(line)
          if ( line(1:5) == "-----" ) then
             exit  ! Dashes mark the end of the data
          else if (trim(line) == "" .or. line(1:1) == "#") then
             cycle ! Ignore whitespace or comments
          else if (n_rows < table_max_rows) then
             n_rows = n_rows + 1
             read(line, FMT = *, ERR = 999, end = 777) temp_table(:, n_rows)
          else
             print *, "CS_read_file error: too many rows in ", &
                  file_name, " at line ", nL
          end if
       end do

       ! Store the data in the actual table
       if (allocated(x_data)) deallocate(x_data)
       if (allocated(y_data)) deallocate(y_data)
       allocate(x_data(n_rows))
       allocate(y_data(n_rows))

       x_data = temp_table(1, 1:n_rows)
       y_data = factor * temp_table(2, 1:n_rows)

       exit                   ! Done
    end do

    close(my_unit)
    return

777 continue ! If the end of the file is reached after finding data
    print *, "table_from_file unexpectedly reached end of " // trim(file_name)
    print *, "searching '" // trim(data_name) // "'"
    call print_usage(file_name, data_name)
    error stop

888 continue ! If the end of the file is reached without finding data
    print *, "table_from_file: no data in " // trim(file_name)
    print *, "searching '" // trim(data_name) // "'"
    call print_usage(file_name, data_name)
    error stop

999 continue ! If there was an input error, the routine will end here
    print *, "table_from_file error at line", nL
    print *, "ioState = ", ioState, " in ", trim(file_name)
    print *, "searching '" // trim(data_name) // "'"
    call print_usage(file_name, data_name)
    error stop

  contains

    subroutine print_usage(file_name, data_name)
      character(len=*), intent(in) :: file_name
      character(len=*), intent(in) :: data_name

      print *, ""
      print *, "Expected a file '", trim(file_name), &
           "' with the following structure:"
      print *, trim(data_name)
      print *, "FACTOR: 1.0         [optional: multiply with this factor]"
      print *, "[other lines]"
      print *, "------------------  [at least 5 dashes]"
      print *, "xxx       xxx       [data in two column format]"
      print *, "...       ..."
      print *, "xxx       xxx"
      print *, "------------------"
      print *, ""
    end subroutine print_usage

  end subroutine table_from_file

  !> Read a list with strings and numbers from a file
  subroutine table_strings_numbers_from_file(filename, data_name, &
       strings, numbers)
    character(len=*), intent(in)                 :: filename
    character(len=*), intent(in)                 :: data_name
    character(len=*), allocatable, intent(inout) :: strings(:)
    real(dp), allocatable, intent(inout)         :: numbers(:)

    character(len=string_len) :: line
    integer                   :: my_unit, n_found, n
    character(len=string_len) :: tmp_strings(table_max_rows)
    real(dp)                  :: tmp_numbers(table_max_rows)

    n_found = 0

    open(newunit=my_unit, file=filename, action="read")

    do
       read(my_unit, "(A)", end=998) line
       line = adjustl(line)

       if (line == data_name) then
          ! Read next line starting with at least 5 dashes
          read(my_unit, "(A)") line
          if (line(1:5) /= "-----") &
               error stop "data name not followed by -----"
          exit
       end if
    end do

    ! Read ignored species, one per line
    do
       read(my_unit, "(A)", end=999) line
       line = adjustl(line)

       ! Ignore comments
       if (line(1:1) == "#") cycle

       ! Exit when we read a line of dashes
       if (line(1:5) == "-----") exit

       n_found = n_found + 1
       read(line, *) tmp_strings(n_found), tmp_numbers(n_found)
    end do

998 close(my_unit)

    allocate(strings(n_found))
    allocate(numbers(n_found))

    do n = 1, n_found
       strings(n) = trim(tmp_strings(n))
       numbers(n) = tmp_numbers(n)
    end do

    return

    ! Error messages
999 error stop "table_names_numbers: no closing dashes"
  end subroutine table_strings_numbers_from_file

end module m_table_data
