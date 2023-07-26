PROGRAM TESTING

        use helperFunc
        use testingFunc
 
        IMPLICIT NONE
        integer :: natom1,natom2,natom,XDIM,i,k, test_failed(5),ft,fileOutNumber
        real*8, allocatable :: ref1(:),ref2(:),mass(:),mass0(:),internal(:),internal0(:)
        real*8 :: pii
        real*8::distance(3)
        INTEGER DATE_TIME (8)
        CHARACTER (LEN = 10) BIG_BEN (3)


        real*8, allocatable :: CasesRand(:,:)
        real*8 :: test_dist(5)
        integer :: counterCase, test_number
        real*8::maxerr


        test_number=100000

        fileOutNumber=12

        ! open(unit=10,file='input.dat',status='old',action='read')

        ! read(10,*) XDIM
        ! read(10,*) natom1! number of atoms in fragment1
        ! read(10,*) natom2! number of atoms in fragment2

!  natom = natom1+natom2 ! total number of atoms

!  allocate(ref1(3*natom1),ref2(3*natom2),mass(natom),mass0(natom),internal(XDIM),internal0(XDIM))

!  if(natom1>1)then
!    do i=1,natom1
!      read(10,*)mass0(i),(ref1(k),k=i*3-2,i*3)
!    enddo
!  elseif(natom1==1)then
!    read(10,*)mass0(1)
!    ref1=0.d0
!  endif

!  if(natom2>1)then
!    do i=1,natom2
!      read(10,*)mass0(natom1+i),(ref2(k),k=i*3-2,i*3)
!    enddo
!  elseif(natom2==1)then
!    read(10,*)mass0(natom1+1)
!    ref2=0.d0
!  endif
!  close(10)
 
!  mass=mass0
!  mass(2)=mass(2)*2d0

        CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), &
        BIG_BEN (3), DATE_TIME)

        call FileChecking('./output/Testing_Output.txt',fileOutNumber)
        REWIND(fileOutNumber)

        write(fileOutNumber,*)"******************************************************************************"
        write(fileOutNumber,*)  "Test Day and Time Record"
        write(fileOutNumber,*)  "Month / Day / Year: ",DATE_TIME(2),"/",DATE_TIME(3),"/",DATE_TIME(1) 
        write(fileOutNumber,*)  "Hr    / Min / Sec : ",DATE_TIME(5),":",DATE_TIME(6),":",DATE_TIME(7) 
        
        call Test_Cart_to_Autosurf('./test/Cartesian/input_H2O.dat',"./test/Cartesian/2b.xyz",fileOutNumber)



        ! call Testing_InteratomicDistances("2DCase/input.dat",2, ft,fileOutNumber)
        ! call Testing_InteratomicDistances("3DCase/input.dat",1, ft,fileOutNumber)
        ! call Testing_InteratomicDistances("3DCase/input.dat",2, ft,fileOutNumber)
        ! call Testing_InteratomicDistances("4DCase/input.dat",2, ft,fileOutNumber)
        ! call Testing_InteratomicDistances("5DCase/input.dat",2, ft,fileOutNumber)
        ! call Testing_InteratomicDistances("6DCase/input.dat",2, ft,fileOutNumber)


        
        ! call Testing_AxisChanges("2DCase/input.dat",2,ft,fileOutNumber)
        ! call Testing_AxisChanges("3DCase/input.dat",2,ft,fileOutNumber)
        ! call Testing_AxisChanges("4DCase/input.dat",2,ft,fileOutNumber)
        ! call Testing_AxisChanges("5DCase/input.dat",2,ft,fileOutNumber)
        ! call Testing_AxisChanges("6DCase/input.dat",2,ft,fileOutNumber)

        close(fileOutNumber)

!  XDIM = 5


!  ALLOCATE(CasesRand(XDIM,test_number + 5000))
 
 

!  call random_seed()
!  call random_number(CasesRand)

 


 
!  Write(*,*)DACOS(1d0 - 1d-16),datan2(1d-16,1d0)

!  Call Same_Masses_Test(CasesRand,XDIM,test_number,counterCase,test_failed,maxerr)



!         write(*,*)
!         write(*,*)"TEST NAME: Same Masses same angle conversion"
!         write(*,*)"Int2Cart = ",0
!         write(*,*)"Number of Tests", counterCase
!         write(*,*)"MAX err_test",maxerr
!         write(*,*)"Failed Tests: ", test_failed, " out of",test_number 
!         write(*,*)
        

   
 !CALL CosineLaw_Test(CasesRand,test_number,mass,mass0,natom1,natom2,ref1,ref2,XDIM,ifun,counterCase,test_failed,maxerr)

 ! write(*,*)
 !       write(*,*)"Cosine Law in the a Molecule A Axis"
 !       write(*,*)"System Features : Dimension/ Number of Atoms", XDIM,"/",natom1,natom2
 !       write(*,*)"Int2Cart = ",ifun
 !       write(*,*)"Number of Tests", counterCase
 !       write(*,*)"MAX err_test",maxerr
 !       write(*,*)"Failed Tests: ", test_failed, " out of",test_number 
 !       write(*,*)

END PROGRAM TESTING











! SUBROUTINE Same_Masses_Test(CasesRand,XDIM,test_number,counterCase,test_failed,maxerr)
!   use helperFunc
!   IMPLICIT NONE
  
!   Integer,INTENT(IN) ::test_number,XDIM
!   integer :: nc
!   real*8 ::internal(XDIM),internal0(XDIM),internal0_BiSph(XDIM)
 
!   real*8 :: pii,passingThrough
!   real*8,INTENT(IN) :: CasesRand(XDIM,test_number+5000)
!   real*8 :: Max_test_dist,test_dist(5),err
!   real*8::err_test,threshold
!   integer,INTENT(INOUT) ::counterCase,test_failed
!   real*8,INTENT(INOUT):: maxerr
!   pii=dacos(-1d0) 
!   err = 1d-7
!   counterCase = 0
!   maxerr = 0d0
!   !threshold = 1d0 - 1d-12
!   nc = 0

!   test_failed = 0
!   Max_test_dist = 0
  
  
  
  
!   do while (counterCase < test_number)
         
!     nc=nc+1
   
    
  
!     internal(1) = CasesRand(1,nc)*10d0+5d0
!     internal(2) = CasesRand(2,nc)*pii
!     internal(3) = CasesRand(3,nc)*pii      
!     internal(4) = CasesRand(4,nc)*2d0*pii
!     internal(5) = CasesRand(5,nc)*2d0*pii 
  
!     if (internal(4) > pii)then
!       internal(4) = internal(4) -  2d0*pii 
!     end if 
!     if (internal(4) < -pii)then
!       internal(4) = internal(4) + 2d0*pii 
!     end if 
!     if (internal(5) > pii)then
!       internal(5) = internal(5) -  2d0*pii 
!     end if 
!     if (internal(5) < -pii)then
!       internal(5) = internal(5) + 2d0*pii 
!     end if
  
!     Call Get_ISOTOP_COORDINATES(internal,internal0,XDIM,'./input_tSM.dat')
   

  
!     counterCase = counterCase + 1
!                 internal0_BiSph = internal0
  
!                 err_test = sum(DABS(internal0_BiSph -internal))
                
!                 if (maxerr < err_test ) then 
!                   maxerr = err_test 
!                 endif
                
!                 if (err_test > err)Then
!                   test_failed = test_failed + 1
!                   write(*,*)
!                   write(*,*)"Internal           ", internal
!                   write(*,*)"Internal 0 BiSph   ", internal0_BiSph
!                   write(*,*)
!                   write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
!                   write(*,*)
               
!                   write (*,*)"R_test" , internal0_BiSph(1) ,internal(1) 
!                   write (*,*) "Error R",DABS(internal0_BiSph(1) -internal(1))
!                   write (*,*)"----------------------------------------------"
!                   write (*,*)"Theta1_Test" , internal0_BiSph(2) ,internal(2) 
!                   write (*,*) "Error Theta1",DABS(internal0_BiSph(2) -internal(2))
!                   write (*,*)"----------------------------------------------"
!                   write (*,*)"Theta2_Test" , internal0_BiSph(3) ,internal(3) 
!                   write (*,*) "Error Theta2",DABS(internal0_BiSph(3) -internal(3))
!                   write (*,*)"----------------------------------------------"
!                   write (*,*)"Phi1_Test" , internal0_BiSph(4) ,internal(4) 
!                   write (*,*) "Error Phi1",DABS(internal0_BiSph(4) -internal(4))
!                   write (*,*)"----------------------------------------------"
!                   write (*,*)"Phi2_Test" , internal0_BiSph(5) ,internal(5) 
!                   write (*,*) "Error Phi2",DABS(internal0_BiSph(5) -internal(5))
!                   write (*,*)"----------------------------------------------"
!                 end if 
        
  
  
!   end do

! END SUBROUTINE  Same_Masses_Test

!
!SUBROUTINE CosineLaw_Test(CasesRand,test_number,mass,mass0,natom1,natom2,ref1,ref2,XDIM,ifun,counterCase,test_failed,maxerr)
!  use helperFunc
!  IMPLICIT NONE
!  
!  Integer,INTENT(IN) ::ifun,test_number,natom1,natom2,XDIM
!  integer :: nc,internalFunction
!  real*8 ::internal(XDIM),internal0(XDIM),internal0_BiSph(XDIM)
!  real*8,INTENT(IN):: ref1(3*natom1),ref2(3*natom2),mass(natom1+natom2),mass0(natom1+natom2)
!  real*8 :: pii,passingThrough
!  real*8,INTENT(IN) :: CasesRand(XDIM,test_number+5000)
!  real*8 :: Max_test_dist,test_dist(5),err
!  real*8::err_test,threshold
!  integer,INTENT(INOUT) ::counterCase,test_failed
!  real*8,INTENT(INOUT):: maxerr
!
!  real*8::R_err,R0_err,theta1_err,theta10_err
!  real*8::cmA(3),cmA0(3),cmB(3),cmB0(3),O(3),O0(3),response1(6),cart(6*3)
!
!
!  pii=acos(-1d0) 
!  err = 10d0**(-8)
!  counterCase = 0
!  maxerr = 0d0
!  internalFunction=ifun
!  threshold = 0.9993d0
!  nc = 0
!
!  test_failed = 0
!  Max_test_dist = 0
!  
!  do while (counterCase < test_number)
!            
!    nc=nc+1
!    passingThrough =0d0
!    
!
!    internal(1) = CasesRand(1,nc)*10d0+5d0
!    internal(2) = CasesRand(2,nc)*pii
!    internal(3) = CasesRand(3,nc)*pii      
!    internal(4) = CasesRand(4,nc)*2d0*pii
!    internal(5) = CasesRand(5,nc)*2d0*pii 
!
!    if (internal(4) > pii)then
!      internal(4) = internal(4) -  2d0*pii 
!    end if 
!    if (internal(4) < -pii)then
!      internal(4) = internal(4) + 2d0*pii 
!    end if 
!    if (internal(5) > pii)then
!      internal(5) = internal(5) -  2d0*pii 
!    end if 
!    if (internal(5) < -pii)then
!      internal(5) = internal(5) + 2d0*pii 
!    end if
!
!    CALL convert_isotopic_coordinates(internal,internal0,mass,mass0,natom1,natom2,ref1,ref2,XDIM,internalFunction,test_dist)
!
!
!            if (dabs(internal0(2))< threshold  .and. test_dist(5)>-1 ) then
!                passingThrough = 1d0
!                if (size(internal0)>3)then
!                    if (dabs(internal0(3)) < threshold)then
!                        passingThrough = 1d0
!                    else
!                        passingThrough = 0d0
!                    end if 
!                end if
!            
!
!            endif
!
!
!            if (passingThrough>0d0  ) then 
!
!    
!          
!             
!                counterCase = counterCase + 1
!                Call ZYZ_to_BiSpherical(internal0,internal0_BiSph)
!
!                call Int_to_Cart_User(internal,cart)
!
!                call cmass2(cart(1:3*natom1),cmA,mass(1:natom1),natom1)
!                call cmass2(cart(1:3*natom1),cmA0,mass0(1:natom1),natom1)
!                call cmass2(cart(1+3*natom1:3*(natom1+natom2)),cmB,mass(1+natom1:natom1+natom2),natom2)
!                call cmass2(cart(1+3*natom1:3*(natom1+natom2)),cmB0,mass0(1+natom1:natom1+natom2),natom2)
!                                  
!                                      
!                O  =   cmA
!                O0 =   cmA0
!                                      
!
!                
!                call CosineLaw(O,O0,cmB,response1)
!
!                R_err  = dabs(response1(2)-internal(1))
!                R0_err = dabs(response1(3)-internal0_BiSph(1))
!                theta10_err = dabs(dabs(response1(5))-dabs(dcos(internal0_BiSph(2))))
!                theta1_err = dabs(dabs(response1(4))-dabs(dcos(internal(2))))
!
!
!                
!                
!
!                err_test = R_err + R0_err + theta10_err + theta1_err
!
!          
!                if (maxerr < err_test ) then 
!                  maxerr = err_test 
!                endif
!                
!                if (err_test > err)Then
!                  test_failed = test_failed + 1
!                  write(*,*)
!                  write(*,*)"Cosine law : ",response1
!                  write(*,*)"R_test  / R:   ",response1(3),internal0_BiSph(1),R_err
!                  write(*,*)"R0_test / R0 : ",response1(2),internal(1),R_err
!                  write(*,*)"Theta1_test / Theta1 : ",dabs(response1(4)),dcos(internal(2)),theta10_err
!                  write(*,*)"Theta10_test / Theta10 : ",dabs(response1(5)),dcos(internal0_BiSph(2)),theta10_err
!                  write (*,*)"----------------------------------------------"
!                end if 
!
!              endif
!        
!
!
!  end do
!
!END SUBROUTINE  CosineLaw_Test
!    
    
