PROGRAM TESTING
  
  IMPLICIT NONE
  

  real*8:: internal(5),internal0(5)
  integer :: XDIM,i
  ! 'interFunc' defines the Int_to_Cart Function 
  ![Options: 0 --> User int_to_Cart , 1 --> Spherical_to_Cart, 2 --> ZYZ_to_Cart  ]
  integer :: interFunc 
  real*8::pii,start, finish
  
  
  pii = DACOS(-1d0)
  
  XDIM = 5
  interFunc = 0

  internal(1) = 5d0
  internal(2) = 20d0*pii/180d0
  internal(3) = 30d0*pii/180d0  
  internal(4) = 40d0*2d0*pii/180d0
  internal(5) = 50d0*2d0*pii/180d0


  Call Get_ISOTOP_COORDINATES(internal,internal0,XDIM,interFunc,'./input.dat')
 
  Write(*,*) "R : ",internal(1),internal0(1)
  Write(*,*) "theta1 :",internal(2),internal0(2),DABS(internal(2)-internal0(2))
  Write(*,*) "theta2 :",internal(3),internal0(3),DABS(internal(3)-internal0(3))
  Write(*,*) "phi1 :",internal(4),internal0(4),DABS(internal(4)-internal0(4))
  Write(*,*) "phi2 :",internal(5),internal0(5),DABS(internal(5)-internal0(5))
  
  

  call cpu_time(start)
  do i=1,100000
    Call Get_ISOTOP_COORDINATES(internal,internal0,XDIM,'./input.dat')
  end do
  call cpu_time(finish)
  
  print '("Time = ",f6.3," seconds.")',finish-start

  
END PROGRAM TESTING
