PROGRAM TESTING

 
 IMPLICIT NONE
 integer :: XDIM


 real*8, allocatable :: CasesRand(:,:)
 real*8 :: test_dist(5)
 integer :: counterCase, test_number
 Integer::test_failed
 real*8::maxerr





 test_number=100000000
 
  write(*,*)


 call Test_Cart_to_Autosurf('input_H2O.dat',"2b.xyz")

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



SUBROUTINE Test_Cart_to_Autosurf(systPath,dataPath)
  use mathFunc
  use testingFunc
  use helperFunc

  IMPLICIT NONE
  character(*),INTENT(IN) :: systPath,dataPath
  real*8::pii,cart(18),cart_model(18),energies(4)
  character(len=1)::Atom_label
  real*8 :: internal0(6),td(4),err_test,err_tolerance,maxerr_1,maxerr_2
  integer :: stat2,stat3,i,natoms,int_to_cart_func,XDIM,ptos(3)
  integer :: failedTest_1,counterCase_1,failedTest_2,counterCase_2
  real*8::response1_model(6),response2_model(6),response1(6),response2(6)

  XDIM=6
  
  
  open(unit=100,file=dataPath,status='old',action='read')
    OPEN(unit=200, FILE="./output/Paesani/coord_H2O_Dymer.txt",ACTION='write',IOSTAT=stat2,STATUS='OLD')
    if (stat2 .ne. 0) then
            open(unit=200,file="./output/Paesani/coord_H2O_Dymer.txt",status='new',action='write')
    end if

    OPEN(unit=300, FILE="./output/Paesani/coord_H2O_Dymer_filtered.txt",ACTION='write',IOSTAT=stat3,STATUS='OLD')
    if (stat3 .ne. 0) then
            open(unit=300,file="./output/Paesani/coord_H2O_Dymer_filtered.txt",status='new',action='write')
    end if

    failedTest_1 = 0
    counterCase_1 = 0
    failedTest_2 = 0
    counterCase_2 = 0
    maxerr_1 = -1d0
    maxerr_2 = -1d0
    err_tolerance=1d-7

    ptos = (/1,2,3 /)


       do i=1,1!,42508  !42508
       
          read(100,*)natoms
          read(100,*)energies(1:4)
          read(100,*)Atom_label,cart(1:3)
          read(100,*)Atom_label,cart(4:6)
          read(100,*)Atom_label,cart(7:9)
          read(100,*)Atom_label,cart(10:12)
          read(100,*)Atom_label,cart(13:15)
          read(100,*)Atom_label,cart(16:18)

          int_to_cart_func = -1 !-1 is for cartesian input

         if (i==1)then
              cart_model=cart
         end if

         Call Get_ISOTOP_COORDINATES(cart,size(cart),internal0,6, int_to_cart_func ,systPath,td)



!Reporting Error

          
          Call Report_Cartesian_Error(counterCase_1,failedTest_1,maxerr_1,td,err_tolerance,natoms,XDIM,cart,internal0)
          call Check_Cartesian_Frames(ptos,cart_model,cart,3,3,counterCase_2,failedTest_2,maxerr_2&
                                      ,err_tolerance,XDIM,internal0)
! End Error Testing 
          
          write(200, *) i , internal0, energies(2)
          if (internal0(1)>5) then 
           write(300, *) i , internal0, energies(2)
          
          endif


   
      enddo
      
      write(*,*) "Number of failed Cartesian coordinates: ",failedTest_1, " out of ", counterCase_1, "maxErr: ",maxerr_1,' ' 
      write(*,*) "Number of failed Check_Cartesian_Frames: ",failedTest_2, " out of ", counterCase_2,"maxErr: ",maxerr_2 ,' '
      
   close(100)
   close(200)
   close(300)

End SUBROUTINE Test_Cart_to_Autosurf





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
    
