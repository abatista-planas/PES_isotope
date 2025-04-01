program example

    use mod_types, only : real64, int32
    use helper_functions,only: get_pes_coordinates
    implicit none

    integer(int32), parameter :: XDIM = 6
    integer(int32), parameter :: ISOTOPE_COORDINATE_SIZE = 6
    real(real64), parameter :: PII = DACOS(-1d0)
    character (len = *, kind = 1), parameter :: ISOTOPE_COORDINATE_FORMAT = "Autosurf"
    character (len = *, kind = 1), parameter :: PES_COORDINATE_FORMAT = "Autosurf"
    character (len = *, kind = 1), parameter :: PATH_TO_FILE = "./test/6DCase/input.dat"

    integer(int32), parameter :: dataset_size = 10
    integer(int32):: i
    real(real64), dimension(6) :: isotope_internal_coordinate
    real(real64), dimension(6,dataset_size) :: isotope_internal_coordinate_array
    real(real64), allocatable :: pes_internal_coordinate(:)


    isotope_internal_coordinate = [ 10d0,&
                                    20d0*PII/180d0,&
                                    30d0*PII/180d0,&
                                    40d0*PII/180d0,&
                                    50d0*PII/180d0,&
                                    60d0*PII/180d0]

    ! Return the coordinates of the calculated PES
    call get_pes_coordinates(pes_internal_coordinate, &
            isotope_internal_coordinate, &
            ISOTOPE_COORDINATE_SIZE, &
            XDIM, & ! System Dimension
            ISOTOPE_COORDINATE_FORMAT, &
            PES_COORDINATE_FORMAT, &
            PATH_TO_FILE) ! File must contain the cartesian coordinates and masses of atoms for each molecule



end program example
