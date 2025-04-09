## Introduction
The study of intermolecular forces is essential for predicting and
understanding the dynamics of molecular systems, which can be relevant to
numerous scientific disciplines, including atmospheric chemistry,
environmental chemistry, and astrochemistry. Recent discoveries of complex
organic and especially chiral molecules in the interstellar medium have drawn
a lot of interest, and theory and simulation are important partners to
experimental measurements in these efforts. For modeling purposes, the study of molecular 
interactions is essential and can be better understood by constructing a potential energy surface (PES) 
for the system of interest. The PES construction starts computing ab-initio energy points for a list of geometrical configurations (tipically the intermolecular distance R, 
and a set the angles that define the relative orientation of the monomers). Every point is very time and computational consuming

Isotopes are molecules that have same atoms but at least one atom different in the number of neutrons. Since the 
neutrons has no charge, if the bond length in both molecules remain constant, the electric field remain invariant. However, the change of mass could change the center of mass and/or the way the principal axis are choseen. To be able to utilize the already computed PES its necessary find a mapping from the the new isosope coordinates to the PES system.

## Features
- Works for any pair of molecules independenly of symmetry's groups.
- The user can enter the coordinate in the Autosurf, Spherical, and Bi-shperical  or define their own custom internal coodinate system by providing the conversion to cartesian coordinates
- The software can be run in parallel 
## Getting Started

### Prerequisites

ReadmeAI requires Python 3.9 or higher, and one of the following installation methods:

| Requirement                          | Details                          |
|--------------------------------------|----------------------------------|
| • Fortran Compiler                   | Core runtime                     |
| • Make (optional)                    | file package manager             |



### Usage


```fortran
program example

    use mod_types, only : real64, int32
    use helper_functions,only: get_pes_coordinates
    implicit none

    integer(int32), parameter :: XDIM = 6
    integer(int32), parameter :: ISOTOPE_COORDINATE_SIZE = 6
    character (len = *, kind = 1), parameter :: ISOTOPE_COORDINATE_FORMAT = "Autosurf"
    character (len = *, kind = 1), parameter :: PES_COORDINATE_FORMAT = "Autosurf"
    character (len = *, kind = 1), parameter :: PATH_TO_FILE = "./test/6DCase/input.dat"

    integer(int32), parameter :: dataset_size = 10
    integer(int32):: i
    real(real64), dimension(6) :: isotope_internal_coordinate
    real(real64), dimension(6,dataset_size) :: isotope_internal_coordinate_array
    real(real64), allocatable :: pes_internal_coordinate(:)

   
    isotope_internal_coordinate = [ 10d0,&            ! R
                                    0.349065850d0,&   ! beta 1
                                    0.523598776d0,&   ! beta 2
                                    0.698131701d0,&   ! alpha
                                    0.872664626d0,&   ! gamma 1 
                                    1.047197551d0]    ! gamma 2

    ! Return the coordinates of the calculated PES
    call get_pes_coordinates(pes_internal_coordinate, &
            isotope_internal_coordinate, &
            ISOTOPE_COORDINATE_SIZE, &
            XDIM, & ! System Dimension
            ISOTOPE_COORDINATE_FORMAT, &
            PES_COORDINATE_FORMAT, &
            PATH_TO_FILE) ! File must contain the cartesian coordinates and masses of atoms for each molecule



end program example
```



### Testing
Perform a test for several different systems symmetries and system dimmensions 
```sh
❯ make test
```



## Credits
Research Group:
- Adrian Batista-Planas
- Ernesto Quintas-Sanchez
- Richard Dawes (Advisor)

This work was partially supported by the Missouri University of Science and Technology’s Kummer Institute for Student Success and the United States Department of Energy (DOE), grant numbers DE-SC0019740 and DE-SC0025420.



