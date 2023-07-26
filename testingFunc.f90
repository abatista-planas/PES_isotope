module testingFunc
  
    contains

    SUBROUTINE Check_Cartesian_Frames(ptos,cart_model,cart,natom1,natom2,counterCase,test_failed,maxerr, &
                                        err_tolerance,XDIM,internal0)
        use mathFunc
        IMPLICIT NONE
        integer,INTENT(INOUT) ::  counterCase,test_failed
        integer,INTENT(IN) ::     XDIM,ptos(3),natom1,natom2
        real*8,INTENT(IN) ::      err_tolerance,internal0(XDIM),cart_model((natom1+natom2)*3),cart((natom1+natom2)*3)
        real*8,INTENT(INOUT) ::   maxerr
        integer:: i,j,k
        real*8::err_test,modelA(6),modelB(6),A(6),B(6)

        counterCase = counterCase + 1

        i =   ptos(1)
        j =   ptos(2)
        k =   ptos(3)

        Call CosineLaw(cart_model(1+3*(i-1):3*i),cart_model(1+3*(j-1):3*j),cart_model(1+3*(k-1):3*k),modelA)
        Call CosineLaw(cart_model(3*natom1+1+3*(i-1):3*natom1+3*i),cart_model(3*natom1+1+3*(j-1):3*natom1+3*j)&
                        ,cart_model(3*natom1+1+3*(k-1):3*natom1+3*k),modelB)

        Call CosineLaw(cart(1+3*(i-1):3*i),cart(1+3*(j-1):3*j),cart(1+3*(k-1):3*k),A)
        Call CosineLaw(cart(3*natom1+1+3*(i-1):3*natom1+3*i),cart(3*natom1+1+3*(j-1):3*natom1+3*j)&
                        ,cart(3*natom1+1+3*(k-1):3*natom1+3*k),B) 


        err_test = max(Norm2(A-modelA),Norm2(B-modelB))

                if (maxerr <  err_test) then 
                    maxerr =  err_test
                endif               
        
            
              if(err_test > err_tolerance) then
                 
                  test_failed = test_failed + 1

                 
                  write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                  write(*,*)"Cartesian coordinates           : ", counterCase
                  write(*,*)"Error           : ", err_test
                  write(*,*)"Norm2 Molec A: ",Norm2(A-modelA),"Norm2 Molec B: ",Norm2(B-modelB),&
                            "Norm2 AB : ",Norm2(A-B)
                            
                  write(*,*)"Model Molec A: ",modelA
                  write(*,*)"Molec A: ",A
                  write(*,*)"Model Molec B: ",modelB
                  write(*,*)"Molec B: ",B
                  write(*,*)"Dimension :",XDIM 
                write(*,*) internal0
                write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                write(*,*)
                            
              end if
              
         
    END SUBROUTINE
    !Reporting Error

    SUBROUTINE Report_Cartesian_Error(counterCase,test_failed,maxerr,td,err_tolerance,natoms,XDIM,cart,internal0)
        
        IMPLICIT NONE
        integer,INTENT(INOUT) ::  counterCase,test_failed
        integer,INTENT(IN) ::     natoms,XDIM
        real*8,INTENT(IN) ::      td(4),err_tolerance,cart(natoms*3),internal0(XDIM)
        real*8,INTENT(INOUT) ::   maxerr
        integer:: j
        real*8::err_test
          counterCase = counterCase + 1
          
          err_test = 0;
          err_test = sum(DABS(td))
          
          if (maxerr < err_test ) then 
            maxerr = err_test 
          endif
          
          if (err_test > err_tolerance)Then
            test_failed = test_failed + 1
            write(*,*)
            write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
            write(*,*)"Cartesian coordinates           : ", counterCase
            write(*,*)"Error           : ", err_test,td
            do j=1,natoms
                write(*,*)"Atom ",j," :", cart(1+(j-1)*3:j*3)
                write (*,*)"----------------------------------------------"
            end do

            write(*,*)"Dimension :",XDIM 
            write(*,*) internal0
            write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
            write(*,*)
          end if 
    END SUBROUTINE Report_Cartesian_Error

    


    ! SUBROUTINE Testing_InteratomicDistances(filename,internalFunction,test_failed,fileOutputNumber)
    !     use helperFunc
    !     use mathFunc
    !     IMPLICIT NONE
    !     integer :: natom1,natom2,natom,XDIM,i,k,nc
    !     real*8, allocatable :: ref1(:),ref2(:),mass(:),mass0(:),internal(:),internal0(:)
    !     real*8 :: pii,threshold,maxdist,passingThrough
    !     real*8, allocatable :: CasesRand(:,:)
    !     real*8 :: Max_test_dist,test_dist(5),err
    !     integer :: counterCase, test_number
    !     character(*),INTENT(IN)::filename
    !     Integer,INTENT(OUT)::test_failed
    !     Integer,INTENT(IN)::internalFunction,fileOutputNumber
        
    !     test_number=10000
    !     open(unit=10,file=filename,status='old',action='read')

    !     read(10,*) XDIM
    !     read(10,*) natom1! number of atoms in fragment1
    !     read(10,*) natom2! number of atoms in fragment2

    !     ALLOCATE(CasesRand(XDIM,test_number + 2000))
    !     pii=acos(-1d0) 
    !     Max_test_dist =0
    !     err = 10d0**(-8)
    !     counterCase = 0

    !     call random_seed()
    !     call random_number(CasesRand)

    !     natom=natom1+natom2 ! total number of atoms

    !     allocate(ref1(3*natom1),ref2(3*natom2),mass(natom),mass0(natom),internal(XDIM),internal0(XDIM))

    !     if(natom1>1)then
    !     do i=1,natom1
    !         read(10,*)mass0(i),(ref1(k),k=i*3-2,i*3)
    !     enddo
    !     elseif(natom1==1)then
    !     read(10,*)mass0(1)
    !     ref1=0.d0
    !     endif

    !     if(natom2>1)then
    !     do i=1,natom2
    !         read(10,*)mass0(natom1+i),(ref2(k),k=i*3-2,i*3)
    !     enddo
    !     elseif(natom2==1)then
    !     read(10,*)mass0(natom1+1)
    !     ref2=0.d0
    !     endif
    !     close(10)
        
    !     ! mass=mass0
    !     ! mass(1)=mass(1)+1d0

    !     call random_number(mass)
    !     mass = mass*50d0

    !         threshold = 0.9993d0
    !         nc = 0

    !         counterCase = 0
    !         test_failed = 0
    !         Max_test_dist = 0
            
            

    !         do while (counterCase < test_number)
    !             nc=nc+1
    !             passingThrough =0d0

    !             internal(1)=CasesRand(1,nc)*10d0 + 5d0
    !             internal(2)=(CasesRand(2,nc)*2d0-1d0)*threshold
    !             if (XDIM > 2) then  
    !                 if (XDIM == 3) then 
    !                     internal(3)=CasesRand(3,nc)*2d0*pii 
    !                 else
    !                     internal(3)=(CasesRand(3,nc)*2d0-1d0)*threshold 
    !                 end if
    !             end if
    !             if (XDIM>3)then  
    !                 internal(4)=CasesRand(4,nc)*2d0*pii 
    !             end if
    !             if (XDIM>4)then  
    !                 internal(5)=CasesRand(5,nc)*2d0*pii 
    !             end if
    !             if (XDIM>5)then  
    !                 internal(6)=CasesRand(6,nc)*2d0*pii 
    !             end if

                
                
                

    !             CALL convert_isotopic_coordinates(internal,internal0,mass,mass0,natom1,natom2,ref1,ref2,XDIM,internalFunction,test_dist)

                
    !             passingThrough = 0d0
    !             if (dabs(internal0(2))< threshold  .and. test_dist(5)>-1 ) then
    !                 passingThrough = 1d0
    !                 if (size(internal0)>3)then
    !                     if (dabs(internal0(3)) < threshold)then
    !                         passingThrough = 1d0
    !                     else
    !                         passingThrough = 0d0
    !                     end if 
    !                 end if
                

    !             endif


    !             if (passingThrough>0d0  ) then 

                    
                    
                    
            
    !                 counterCase = counterCase + 1
    !                 maxdist = MAXVAL(test_dist);
    !                 if (maxdist > Max_test_dist ) then 
    !                     Max_test_dist = maxdist 
    !                 endif
                
    !                 if (maxdist > err  ) then 
    !                     write(fileOutputNumber,*)"--------------------------------------- "
    !                     write(fileOutputNumber,*)"Intermolecular Distances Test : ", test_dist
    !                     write(fileOutputNumber,*)"Internal : ", internal
    !                     write(fileOutputNumber,*)"Internal0 : ",internal0
    !                     write(fileOutputNumber,*)"----------------------------------------"
    !                     test_failed = test_failed + 1

    !                 endif
    !             endif
    !             !End Testing Section
    !     enddo 
    !         write(fileOutputNumber,*)
    !         write(fileOutputNumber,*)"TEST NAME: Testing_InteratomicDistances"
    !         write(fileOutputNumber,*)"System Features : Dimension/ Number of Atoms", XDIM,"/",natom1,natom2
    !         write(fileOutputNumber,*)"Max_test_dist : ", Max_test_dist
    !         write(fileOutputNumber,*)"Int2Cart = ",internalFunction
    !         write(fileOutputNumber,*)"Number of Tests", counterCase
    !         write(fileOutputNumber,*)"Failed Tests: ", test_failed, " out of",test_number 
    !         write(fileOutputNumber,*)
    

    ! END SUBROUTINE Testing_InteratomicDistances




    ! SUBROUTINE Testing_AxisChanges(filename,internalFunction,test_failed,fileOutputNumber)
    !     use mathFunc
    !     use helperFunc
    !     IMPLICIT NONE
    !     integer :: natom1,natom2,natom,XDIM,i,k,nc
    !     real*8, allocatable :: ref1(:),ref2(:),mass(:),mass0(:),internal(:),internal0(:)
    !     real*8 :: pii,threshold,maxdist,passingThrough
    !     real*8, allocatable :: CasesRand(:,:),MassRand(:)
    !     real*8 :: Max_test_dist,test_dist(5),err,dcmA(3),dcmB(3)
    !     integer :: counterCase, test_number
    !     character(*),INTENT(IN)::filename
    !     Integer,INTENT(OUT)::test_failed
    !     Integer,INTENT(IN)::internalFunction,fileOutputNumber
    !     real*8, allocatable ::ref1_temp0(:),ref2_temp0(:),cart(:)
    !     real*8::cmA(3),cmA0(3),cmB(3),cmB0(3)

    !     real*8:: R0,R_,cb1,cb10,cb2,cb20,R0_test,cb1_test,cb10_test,cb2_test,cb20_test
    !     real*8:: R_error,cb10_error,cb1_error,cb2_error,cb20_error,O(3),O0(3),B1(3)
    !     real*8::response1(6),response2(6),response3(6)
        
    !     test_number=10000
    !     open(unit=10,file=filename,status='old',action='read')

    !     read(10,*) XDIM
    !     read(10,*) natom1! number of atoms in fragment1
    !     read(10,*) natom2! number of atoms in fragment2

    !     ALLOCATE(CasesRand(XDIM,test_number + 5000),MassRand(test_number + 5000))
    !     pii=acos(-1d0) 
    !     Max_test_dist =0
    !     err = 10d0**(-8)
    !     counterCase = 0

    !     call random_seed()
    !     call random_number(CasesRand)
    !     call random_number(MassRand)
    !     natom=natom1+natom2 ! total number of atoms

    !     allocate(ref1(3*natom1),ref2(3*natom2),mass(natom),mass0(natom),internal(XDIM),internal0(XDIM))
    !     allocate(ref1_temp0(3*natom1),ref2_temp0(3*natom2),cart(natom*3))
    !     if(natom1>1)then
    !     do i=1,natom1
    !         read(10,*)mass0(i),(ref1(k),k=i*3-2,i*3)
    !     enddo
    !     elseif(natom1==1)then
    !     read(10,*)mass0(1)
    !     ref1=0.d0
    !     endif

    !     if(natom2>1)then
    !     do i=1,natom2
    !         read(10,*)mass0(natom1+i),(ref2(k),k=i*3-2,i*3)
    !     enddo
    !     elseif(natom2==1)then
    !     read(10,*)mass0(natom1+1)
    !     ref2=0.d0
    !     endif
    !     close(10)
        
    !     mass=mass0
        

    !         threshold = 0.9993d0
    !         nc = 0

    !         counterCase = 0
    !         test_failed = 0
    !         Max_test_dist = 0
            
            

    !         do while (counterCase < test_number)
                
    !             nc=nc+1
    !             passingThrough =0d0
    
    !             mass(1)=mass0(1)+MassRand(nc)*10d0

    !             internal(1)=CasesRand(1,nc)*10d0 + 5d0
    !             internal(2)=(CasesRand(2,nc)*2d0-1d0)*threshold

    !             if (XDIM > 2) then  
    !                 if (XDIM == 3) then 
    !                     internal(3)=CasesRand(3,nc)*2d0*pii 
    !                 else
    !                     internal(3)=(CasesRand(3,nc)*2d0-1d0)*threshold 
    !                 end if
    !             end if
    !             if (XDIM>3)then  
    !                 internal(4)=CasesRand(4,nc)*2d0*pii 
    !             end if
    !             if (XDIM>4)then  
    !                 internal(5)=CasesRand(5,nc)*2d0*pii 
    !             end if
    !             if (XDIM>5)then  
    !                 internal(6)=CasesRand(6,nc)*2d0*pii 
    !             end if

                
                
                

    !         CALL convert_isotopic_coordinates(internal,internal0,mass,mass0,natom1,natom2,ref1,ref2,XDIM,internalFunction,test_dist)



    !         passingThrough = 0d0
    !             if (dabs(internal0(2))< threshold  .and. test_dist(5)>-1 ) then
    !                 passingThrough = 1d0
    !                 if (size(internal0)>3)then
    !                     if (dabs(internal0(3)) < threshold)then
    !                         passingThrough = 1d0
    !                     else
    !                         passingThrough = 0d0
    !                     end if 
    !                 end if
                

    !             endif


    !             if (passingThrough>0d0  ) then 

                    
                    
                    
            
                    
    !                 maxdist = MAXVAL(test_dist);
    !                 if (maxdist > Max_test_dist ) then 
    !                     Max_test_dist = maxdist 
    !                 endif
                
    !                 if (maxdist < err  ) then 
    !                     counterCase = counterCase + 1

    !                     R_error = 0d0
    !                     cb1_error = 0d0
    !                     cb10_error = 0d0
    !                     cb2_error = 0d0
    !                     cb20_error = 0d0

    !                                 ! make temporary copies of the input data
    !                         ref1_temp0=ref1! reference vector for frag1 (original frame)
    !                         ref2_temp0=ref2! reference vector for frag2 (original frame)
                            
    !                         ! make sure the original Cartesian frame of reference is at the CM of each fragment 
    !                         if(natom1>1)call rm_cmass(ref1_temp0,mass0(1:natom1),natom1,natom1)
    !                         if(natom2>1)call rm_cmass(ref2_temp0,mass0(natom1+1:natom),natom2,natom2)  
                            
    !                         !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of 
    !                         !!!    reference (origin at the new_CM of frag. 1)
                            
    !                         ! subroutine INT_Cart (convert: Internal coordinates --> Cartesian coordinates)  INT_Cart(internal,cart)
    !                         ref1=ref1_temp0
    !                         ref2=ref2_temp0
    !                         ! (for CH3CN-He system: internal coordinates taken as spherical coords.)
                            

    !                         call Int_to_Cart(internal,mass,ref1,ref2,XDIM,natom1,natom2,internalFunction,cart)

    !                         call cmass2(cart(1:3*natom1),cmA,mass(1:natom1),natom1)
    !                         call cmass2(cart(1:3*natom1),cmA0,mass0(1:natom1),natom1)
    !                         call cmass2(cart(1+3*natom1:3*(natom1+natom2)),cmB,mass(1+natom1:natom1+natom2),natom2)
    !                         call cmass2(cart(1+3*natom1:3*(natom1+natom2)),cmB0,mass0(1+natom1:natom1+natom2),natom2)
                        
                            
    !                         O  =   cmA
    !                         O0 =   cmA0
                            

    !                         call CosineLaw(O,O0,cmB,response1)

                    


    !                         R0  = internal0(1)
    !                         R_ = internal(1)
    !                         cb10 = internal0(2)
    !                         cb1 = internal(2)
                            
                            
    !                         R0_test = response1(3)
                            
    !                         cb1_test = response1(4)
    !                         cb10_test = response1(5)
                
                            
    !                         R_error = dabs(R0-R0_test)
    !                         cb1_error = dabs(dabs(cb1)-dabs(cb1_test))
    !                         cb10_error = dabs(dabs(cb10)-dabs(cb10_test))

    !                         ! write(*,*)"R comparation: ", R0,R0_test,"error: ", R_error
    !                         ! write(*,*)"cb1 comparation: ", cb1,cb1_test,"error: ", cb1_error
    !                         ! write(*,*)"cb10 comparation: ", cb10,cb10_test,"error: ", cb10_error

                            

    !                         if (XDIM > 3)then

    !                                 cb20 = internal0(3)
    !                                 cb2 = internal(3)
    !                                 B1 =   cart(1+3*natom1:3+3*natom1)

    !                                 if (NORM2(B1-cmB)<err)then
    !                                     B1 =   cart(4+3*natom1:6+3*natom1)
    !                                 end if

    !                                 call CosineLaw(O0,cmB,B1,response2)
    !                                 cb20_test =-1d0*response2(5)
    !                                 cb20_error = dabs(dabs(cb20)-dabs(cb20_test))

    !                                 !write(*,*)"cb20 comparation: ",cb20,cb20_test,"error: ", cb20_error

    !                                 call CosineLaw(O,cmB,B1,response3)
    !                                 cb2_test =-1d0*response3(5)
    !                                 cb2_error = dabs(dabs(cb2)-dabs(cb2_test))

    !                                 !write(*,*)"cb2 comparation: ",cb2,cb2_test,"error: ", cb2_error
    !                         end if

    !                     if (R_error > err .or. cb1_error > err .or. cb10_error > err &
    !                     .or. cb2_error > err .or. cb20_error > err) then 
    !                                 write(fileOutputNumber,*)"--------------------------------------- "
    !                                 write(fileOutputNumber,*)"TestDistance :",test_dist
    !                                 write(fileOutputNumber,*)"R comparation: ", R0,R0_test,"error: ", R_error
    !                                 write(fileOutputNumber,*)"cb1 comparation: ", cb1,cb1_test,"error: ", cb1_error
    !                                 write(fileOutputNumber,*)"cb10 comparation: ", cb10,cb10_test,"error: ", cb10_error
    !                                 write(fileOutputNumber,*)"cb20 comparation: ",cb20,cb20_test,"error: ", cb20_error
    !                                 write(fileOutputNumber,*)"cb2 comparation: ",cb2,cb2_test,"error: ", cb2_error
    !                                 write(fileOutputNumber,*)"Error Tests R,cb1,cb10,cb2,cb20: ", &
    !                                 R_error,cb1_error,cb10_error,cb2_error,cb20_error
    !                                 write(fileOutputNumber,*)"Internal : ", internal
    !                                 write(fileOutputNumber,*)"Internal0 : ",internal0
    !                                 write(fileOutputNumber,*)"----------------------------------------"
    !                                 test_failed = test_failed + 1
                
    !                     endif
    !                 endif
    !             endif




        
        
    !             !End Testing Section
    !     enddo 
    !         write(fileOutputNumber,*)
    !         write(fileOutputNumber,*)"TEST NAME: Testing_AxisChanges"
    !         write(fileOutputNumber,*)"System Features : Dimension/ Number of Atoms", XDIM,"/",natom1,natom2
    !         write(fileOutputNumber,*)"Int2Cart = ",internalFunction
    !         write(fileOutputNumber,*)"Number of Tests", counterCase
    !         write(fileOutputNumber,*)"Failed Tests: ", test_failed, " out of",test_number 
    !         write(fileOutputNumber,*)
    

    ! END SUBROUTINE Testing_AxisChanges








end module testingFunc
    