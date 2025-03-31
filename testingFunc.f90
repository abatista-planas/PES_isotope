module testingFunc
  
    contains

    SUBROUTINE Check_Cartesian_Frames(ptos,cart_model,cart,natom1,natom2,counterCase,test_failed,maxerr, &
                                        err_tolerance,XDIM,internal0,fileOutputNumber)
        use mathFunc
        IMPLICIT NONE
        integer,INTENT(INOUT) ::  counterCase,test_failed
        integer,INTENT(IN) ::     XDIM,ptos(3),natom1,natom2
        real*8,INTENT(IN) ::      err_tolerance,internal0(XDIM),cart_model((natom1+natom2)*3),cart((natom1+natom2)*3)
        real*8,INTENT(INOUT) ::   maxerr
        integer:: i,j,k
        real*8::err_test,modelA(6),modelB(6),A(6),B(6)
        integer,optional:: fileOutputNumber
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

                 if (present(fileOutputNumber))then
                  write(fileOutputNumber,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                  write(fileOutputNumber,*)"Cartesians Rigid Molec Checker            : ", counterCase
                  write(fileOutputNumber,*)"Error           : ", err_test
                  write(fileOutputNumber,*)"Norm2 Molec A: ",Norm2(A-modelA),"Norm2 Molec B: ",Norm2(B-modelB),&
                            "Norm2 AB : ",Norm2(A-B)
                            
                  write(fileOutputNumber,*)"Model Molec A: ",modelA
                  write(fileOutputNumber,*)"Molec A: ",A
                  write(fileOutputNumber,*)"Model Molec B: ",modelB
                  write(fileOutputNumber,*)"Molec B: ",B
                  write(fileOutputNumber,*)"Dimension :",XDIM 
                  write(fileOutputNumber,*) internal0
                  write(fileOutputNumber,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                  write(fileOutputNumber,*)
                 
                 else
                  write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                  write(*,*)"Cartesians Rigid Molec Checker           : ", counterCase
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
                  
                            
              end if
              
         
    END SUBROUTINE
    !Reporting Error

    SUBROUTINE Report_Cartesian_Error(counterCase,test_failed,maxerr,td,err_tolerance,natoms,XDIM,cart,internal0&
                                      ,fileOutputNumber)
        
        IMPLICIT NONE
        integer,INTENT(INOUT) ::  counterCase,test_failed
        integer,INTENT(IN) ::     natoms,XDIM
        real*8,INTENT(IN) ::      td(4),err_tolerance,cart(natoms*3),internal0(XDIM)
        real*8,INTENT(INOUT) ::   maxerr
        integer,optional:: fileOutputNumber
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

            if (present(fileOutputNumber))then
              write(fileOutputNumber,*)
              write(fileOutputNumber,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
              write(fileOutputNumber,*)"Cartesians-Autosurf  InterAtomic-Distance Test          : ", counterCase
              write(fileOutputNumber,*)"Error           : ", err_test,td
              do j=1,natoms
                  write(fileOutputNumber,*)"Atom ",j," :", cart(1+(j-1)*3:j*3)
                  write (fileOutputNumber,*)"----------------------------------------------"
              end do

              write(fileOutputNumber,*)"Dimension :",XDIM 
              write(fileOutputNumber,*) internal0
              write(fileOutputNumber,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
              write(fileOutputNumber,*)

            else

              write(*,*)
              write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
              write(*,*)"Cartesians-Autosurf  InterAtomic-Distance Test           : ", counterCase
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
          end if 
    END SUBROUTINE Report_Cartesian_Error

    
    SUBROUTINE Test_Cart_to_Autosurf(systPath,dataPath,fileOutputNumber)
      use mathFunc
      use helperFunc

      IMPLICIT NONE
      character(*),INTENT(IN) :: systPath,dataPath
      integer,optional:: fileOutputNumber

      real*8::cart(18),cart_model(18),energies(4)
      character(len=1)::Atom_label
      real*8 :: td(4),err_test,err_tolerance,maxerr_1,maxerr_2
      integer :: stat2,stat3,i,natoms,int_to_cart_func,XDIM,ptos(3)
      integer :: failedTest_1,counterCase_1,failedTest_2,counterCase_2
      real*8::response1_model(6),response2_model(6),response1(6),response2(6)
      integer :: stat
      real*8,allocatable:: internal0(:)

      XDIM=6
      
      
      open(unit=100,file=dataPath,status='old',action='read')
        ! OPEN(unit=200, FILE="./output/Paesani/coord_H2O_Dymer.txt",ACTION='write',IOSTAT=stat2,STATUS='OLD')
        ! if (stat2 .ne. 0) then
        !         open(unit=200,file="./output/Paesani/coord_H2O_Dymer.txt",status='new',action='write')
        ! end if

        ! OPEN(unit=300, FILE="./output/Paesani/coord_H2O_Dymer_filtered.txt",ACTION='write',IOSTAT=stat3,STATUS='OLD')
        ! if (stat3 .ne. 0) then
        !         open(unit=300,file="./output/Paesani/coord_H2O_Dymer_filtered.txt",status='new',action='write')
        ! end if

        failedTest_1 = 0
        counterCase_1 = 0
        failedTest_2 = 0
        counterCase_2 = 0
        maxerr_1 = -1d0
        maxerr_2 = -1d0
        err_tolerance=1d-7
        
        ptos = (/1,2,3 /)
        i=0

          do !i=1,2510!42508  !42508
            
            i=i+1
            
      
          
              read(100,*,iostat=stat)natoms
              if (stat /= 0 .or. i>2512) exit

              read(100,*)energies(1:4)
              read(100,*)Atom_label,cart(1:3)
              read(100,*)Atom_label,cart(4:6)
              read(100,*)Atom_label,cart(7:9)
              read(100,*)Atom_label,cart(10:12)
              read(100,*)Atom_label,cart(13:15)
              read(100,*)Atom_label,cart(16:18)

              !-1 is for cartesian input

            if (i==1)then
                  cart_model=cart
            end if

      
            Call Get_ISOTOP_COORDINATES(cart,size(cart),internal0,6, "Cartesian", "Autosurf" ,systPath,testArr_Errors = td)
  !  internal,internalLength,internal0_,XDIM,inputCoord,outputCoord,filePath,&
  !                                     testArr_Errors,newflag


        !Reporting Error

              
              Call Report_Cartesian_Error(counterCase_1,failedTest_1,maxerr_1,td,err_tolerance,natoms,XDIM,cart,internal0,&
              fileOutputNumber)
              call Check_Cartesian_Frames(ptos,cart_model,cart,3,3,counterCase_2,failedTest_2,maxerr_2&
                                        ,err_tolerance,XDIM,internal0,fileOutputNumber)
        !End Error Testing 
              
            !   write(200, *) i , internal0, energies(2)
            !   if (internal0(1)>5) then 
            !    write(300, *) i , internal0, energies(2)
              
            !   endif


      
          enddo
          
            if (present(fileOutputNumber))then
                    write(fileOutputNumber,*)" "
                    write(fileOutputNumber,*) "%%%%%%%%%%%%%%%    Cartesian Test    %%%%%%%%%%%%%%%"
                    write(fileOutputNumber,*) "Number of failed Cartesian coordinates: ",failedTest_1, " out of ",&
                                            counterCase_1, ", maxErr: ",maxerr_1,' ' 
                    write(fileOutputNumber,*) "Number of failed Check_Cartesian_Frames: ",failedTest_2, " out of ",&
                                            counterCase_2,", maxErr: ",maxerr_2 ,' '
                    write(fileOutputNumber,*)" "
            endif




        call Print_Test_Results("Cartesians-Autosurf  InterAtomic-Distance Test",maxerr_1,"Cartesian", "Autosurf",&
                                counterCase_1,counterCase_1 - failedTest_1)
        call Print_Test_Results("Cartesians Rigid Molec Checker                ",maxerr_2,"Cartesian", "Autosurf",&
                                counterCase_2,counterCase_2 - failedTest_2)

          
      close(100)
      !  close(200)
      !  close(300)

    End SUBROUTINE Test_Cart_to_Autosurf

  
    SUBROUTINE Print_Test_Results(testName,maxErr,inputCoord,outputCoord,ntest,success_test)

        IMPLICIT NONE
        character(*),INTENT(IN) :: testName,inputCoord,outputCoord
        integer,INTENT(IN)::ntest,success_test
        character (len =25 ) ::f_res1,f_res2
        character (len =7) :: str,str1
        integer::failed,ptg
        real*8,INTENT(IN)::maxErr


        f_res1 = achar(27)//'[32m TEST PASSED '//achar(27)//'[0m'
        f_res2 = achar(27)//'[31m TEST FAILED :('//achar(27)//'[0m'

        ptg = FLOOR(success_test*100d0/ntest)
        failed = ntest-success_test
        write (str1, '(I4,A)') ptg,"%"

        write(*,  *)
        if (ntest > success_test)then
    
          write(*, '(A, I5, A, I1, A, I5, A)') f_res2//"("//achar(27)//'[31m'//str1//achar(27)//'[0m)'//" ----  "&
                  //achar(27)//'[35m'//testName//achar(27)//'[0m)'
        else
          write(*,  *) f_res1//"("//achar(27)//'[32m 100% '//achar(27)//'[0m)'//" ----  "&
          //achar(27)//'[35m'//testName//achar(27)//'[0m)'                 
        end if
        write(*,  '(A, I8, A, I8, A, I8)') "                                Num. Test: ",ntest,"     ---> success: "&
                                          ,success_test," / failed : ",failed
        write(*,  '(A, F21.15)') "                                Error Max: ",maxErr                                 
        write(*,  *) "                               Input/Output Coordinates: "//achar(27)//'[33m'//inputCoord//" / "&
                            //outputCoord//achar(27)//'[0m)'
        write(*,  *)

    end SUBROUTINE Print_Test_Results


<<<<<<< HEAD
    SUBROUTINE Generate_RandomData(internal,szi,rms,new_mass,mass,mass0,natom1,natom2,ref1_0,&
                                    ref2_0,XDIM,inputCoord,ntest)
        use helperFunc
        use coordinateTransf

        IMPLICIT NONE
        Integer,INTENT(IN)::XDIM,rms,natom1,natom2,ntest,szi ! rm = random mass
        real*8, INTENT(OUT) :: internal(szi,ntest)
        real*8, INTENT(OUT) :: new_mass(natom1+natom2,ntest)
        real*8 :: pii,threshold,maxdist,passingThrough,iter(XDIM),rmin,rmax
        real*8, INTENT(IN) :: mass(natom1+natom2),mass0(natom1+natom2),ref1_0(natom1*3),ref2_0(natom2*3)
        character(*),INTENT(IN)::inputCoord
        real*8, ALLOCATABLE::CasesRand(:,:),cart(:)
        real*8::ref1_temp0(natom1*3),ref2_temp0(natom2*3),nm(natom1+natom2,ntest),Autosurf(XDIM,ntest)
             
        integer::i,natom,sz,nseeds
    

        natom = natom1+natom2
        pii=acos(-1d0) 
        sz = rms*(natom )+XDIM
        rmin=5d0
        rmax=15d0
        threshold = 0.9993d0

        if (rms==1)then
            
            sz = rms*(natom )+XDIM
        else
            sz = XDIM
        end if 
        ! make temporary copies of the input data
        ref1_temp0=ref1_0! reference vector for frag1 (original frame)
        ref2_temp0=ref2_0! reference vector for frag2 (original frame)

        ALLOCATE(CasesRand(sz,ntest))

        call random_seed()
        call random_number(CasesRand)
        
        nseeds = ntest+0

        call MassGenerator(rms,CasesRand,sz,nseeds,XDIM,natom,mass,mass0,nm)

        new_mass = nm
        
   
        do i = 1,ntest
            !write(*,*)(sz-XDIM)+1,sz,CasesRand(sz-XDIM+1:sz,i)
            !iter = CasesRand(sz-XDIM+1:sz,i)
            Autosurf(1,i)=CasesRand(sz-XDIM+1,i)*10d0 + 5d0
            Autosurf(2,i)=(CasesRand(sz-XDIM+2,i)*2d0-1d0)*threshold
          if (XDIM > 2) then  
              if (XDIM == 3) then 
                  Autosurf(3,i)=CasesRand(sz-XDIM+3,i)*2d0*pii 
              else
                  Autosurf(3,i)=(CasesRand(sz-XDIM+3,i)*2d0-1d0)*threshold 
              end if
          end if
                if (XDIM>3)then  
                    Autosurf(4,i)=CasesRand(sz-XDIM+4,i)*2d0*pii 
                end if
                if (XDIM>4)then  
                    Autosurf(5,i)=CasesRand(sz-XDIM+5,i)*2d0*pii 
                end if
                if (XDIM>5)then  
                    Autosurf(6,i)=CasesRand(sz-XDIM+6,i)*2d0*pii 
                end if
        end do 

        
        if ( inputCoord=="Autosurf" )then
            internal = Autosurf
        elseif ( inputCoord == "Cartesian")then 
          do i = 1,ntest
            call Int_to_Cart(Autosurf(:,i),XDIM,new_mass(:,i),ref1_temp0,ref2_temp0,XDIM,natom1,natom2,"Autosurf",cart)
            internal(:,i) = cart
          end do
        endif

       deallocate(CasesRand)
        

    end Subroutine Generate_RandomData

    Subroutine MassGenerator(rm,seedArr,sz1,sz2,XDIM,natom,mass,mass0,nmass)
        IMPLICIT NONE
        Integer,INTENT(IN)::sz1,sz2,natom,rm,XDIM
        real*8, INTENT(IN) :: seedArr(sz1,sz2),mass(natom),mass0(natom)
        real*8, INTENT(INOUT) :: nmass(natom,sz2)
        real*8::   rand(natom),fact_mass
        integer::i,j



        
        fact_mass = 0.2d0 ! maximun percentage of mass changing

          do i = 1,sz2
           
            if (rm==1)then
              
                rand = seedArr(1:natom,i)*2d0 - 1d0
                do j = 1,natom
                  nmass(j,i) = mass0(j)*(1d0+rm*rand(j)*fact_mass)
                end do
            elseif (rm==0)then
              
                nmass(1:natom,i) = mass0(1:natom)
            elseif (rm==-1)then
        
                nmass(:,i) = mass                 
            end if 
          end do

    end Subroutine MassGenerator

=======
    SUBROUTINE Generate_RandomData(internal,rm,mass,mass0,natom1,natom2,ref1_0,ref2_0,XDIM,inputCoord,ntest)
        use helperFunc

        IMPLICIT NONE
        Integer,INTENT(IN)::XDIM,rm,natom1,natom2,ntest ! rm = random mass
        real*8, ALLOCATABLE,INTENT(OUT) :: internal(:)
        real*8 :: pii,threshold,maxdist,passingThrough
        real*8, INTENT(IN) :: mass(natom1+natom2),mass0(natom1+natom2),ref1_0(natom1*3),ref2_0(natom2*3)
        character(*),INTENT(IN)::inputCoord
        real*8, ALLOCATABLE::CasesRand(:,:)

    

   
        ALLOCATE(CasesRand(rm*(natom1+natom2)+XDIM,ntest + 2000))

        if (inputCoord == "Cartesian")then
            ALLOCATE(internal(3*(natom1+natom2),ntest + 2000))
        else
            ALLOCATE(internal(XDIM,ntest + 2000))
        end if
        

        
        pii=acos(-1d0) 
   

        call random_seed()
        call random_number(CasesRand)


    end Subroutine Generate_RandomData

>>>>>>> 19b107b8a1b92c01564ed2550ae60ef0d3bec0e5

    

    SUBROUTINE Testing_InteratomicDistances(filename,XDIM,inputCoord,outputCoord,test_failed,fileOutputNumber)
        use helperFunc

        IMPLICIT NONE
        integer :: i,k,nc
        real*8, ALLOCATABLE :: internal(:),internal0(:)
        real*8 :: pii,threshold,maxdist,passingThrough
        real*8, allocatable :: CasesRand(:,:)
        real*8 :: Max_test_dist,testArr_Errors(4),err
        integer :: counterCase, test_number
        character(*),INTENT(IN)::filename,inputCoord,outputCoord
        Integer,INTENT(OUT)::test_failed
        Integer,INTENT(IN)::XDIM
        integer,optional:: fileOutputNumber
        character(len=50)::testName
    

        test_number=10
   

        ALLOCATE(CasesRand(XDIM,test_number + 2000))
        pii=acos(-1d0) 
        Max_test_dist =0
        err = 10d0**(-8)
        counterCase = 0

        call random_seed()
        call random_number(CasesRand)


        allocate(internal(XDIM))!,internal0(XDIM))



            threshold = 0.9993d0
            nc = 0

            counterCase = 0
            test_failed = 0
            Max_test_dist = 0
            
            

            do while (counterCase < test_number)
                nc=nc+1
                passingThrough =0d0

                internal(1)=CasesRand(1,nc)*10d0 + 5d0
                internal(2)=(CasesRand(2,nc)*2d0-1d0)*threshold
                if (XDIM > 2) then  
                    if (XDIM == 3) then 
                        internal(3)=CasesRand(3,nc)*2d0*pii 
                    else
                        internal(3)=(CasesRand(3,nc)*2d0-1d0)*threshold 
                    end if
                end if
                if (XDIM>3)then  
                    internal(4)=CasesRand(4,nc)*2d0*pii 
                end if
                if (XDIM>4)then  
                    internal(5)=CasesRand(5,nc)*2d0*pii 
                end if
                if (XDIM>5)then  
                    internal(6)=CasesRand(6,nc)*2d0*pii 
                end if


                Write(*,*)internal,XDIM,inputCoord,outputCoord,filename
              
                
                call Get_ISOTOP_COORDINATES(internal,size(internal),internal0,XDIM,inputCoord,outputCoord,filename,testArr_Errors)

                
                passingThrough = 0d0
                if (dabs(internal0(2))< threshold  ) then
                    passingThrough = 1d0
                    if (size(internal0)>3)then
                        if (dabs(internal0(3)) < threshold)then
                            passingThrough = 1d0
                        else
                            passingThrough = 0d0
                        end if 
                    end if
                

                endif


                if (passingThrough>0d0  ) then 

                    
                    
                    
            
                    counterCase = counterCase + 1
                    maxdist = MAXVAL(testArr_Errors);
                    if (maxdist > Max_test_dist ) then 
                        Max_test_dist = maxdist 
                    endif
                
                    if (maxdist > err  ) then 
                        write(fileOutputNumber,*)"--------------------------------------- "
                        write(fileOutputNumber,*)"Intermolecular Distances Test : ", testArr_Errors
                        write(fileOutputNumber,*)"Internal : ", internal
                        write(fileOutputNumber,*)"Internal0 : ",internal0
                        write(fileOutputNumber,*)"----------------------------------------"
                        test_failed = test_failed + 1

                    endif
                endif
                !End Testing Section
        enddo 

        !testName = "Interatomic Distances("XDIM//"D), ifun: "//internalFunction
        write(testName,'(A(I1)A)') "Interatomic Distances(",XDIM,"D)"
        testName = trim(testName);

        call Print_Test_Results(testName,Max_test_dist,inputCoord,outputCoord,counterCase,counterCase - test_failed)
            write(fileOutputNumber,*)
            write(fileOutputNumber,*)"TEST NAME: ", testName
            write(fileOutputNumber,*)"System Features : Dimension/ Number of Atoms", XDIM
            write(fileOutputNumber,*)"Max_test_dist : ", Max_test_dist
            write(fileOutputNumber,*)"Input/Output Coordinates: ",inputCoord," / ",outputCoord
            write(fileOutputNumber,*)"Number of Tests", counterCase
            write(fileOutputNumber,*)"Failed Tests: ", test_failed, " out of",test_number 
            write(fileOutputNumber,*)
      

    END SUBROUTINE Testing_InteratomicDistances

<<<<<<< HEAD
    SUBROUTINE Testing_InteratomicDistances_v2(filename,XDIM,inputCoord,outputCoord,test_failed,rm,fileOutputNumber)
=======
    SUBROUTINE Testing_InteratomicDistances_v2(filename,XDIM,inputCoord,outputCoord,test_failed,fileOutputNumber)
>>>>>>> 19b107b8a1b92c01564ed2550ae60ef0d3bec0e5
        use helperFunc
        use mathFunc
        use coordinateTransf
  
      
        IMPLICIT NONE
<<<<<<< HEAD

     

        
        character(*),INTENT(IN)::filename,inputCoord,outputCoord
        Integer,INTENT(OUT)::test_failed
        Integer,INTENT(IN)::XDIM
        Integer,INTENT(IN)::rm! random mass switcher  rm =1 random mass, rm =0 new mass = mass0, rm =-1 new_mass=mass of file
        integer,optional:: fileOutputNumber

        integer :: natom1,natom2,natom
        real*8, allocatable :: ref1_0(:),ref2_0(:),new_mass(:,:),nmass(:),mass(:),mass0(:),&
                                ref1_temp0(:),ref2_temp0(:),cart(:)
        real*8 ::  R_ZYZ(7)
        Integer :: Xdim_file,internalLength
        integer :: i,k,nc
        real*8, ALLOCATABLE :: internal(:,:),internal0(:)
        real*8 :: pii,threshold,maxdist,passingThrough
        real*8 :: Max_test_dist,testArr_Errors(4),err
        character(len=50)::testName
        character(len=20)::massStatement
        real*8 ::  td(4),internal_i(XDIM)
        integer :: counterCase, ntest,internal0_length,szi
        real*8,allocatable:: internal0_(:),new_mass_(:,:)

        

        natom = natom1 + natom2
        threshold = 0.9993d0
        nc = 0
        pii=acos(-1d0)
        err = 10d0**(-8)
        counterCase = 0
        test_failed = 0
        Max_test_dist = 0

        

        Call Read_File(filename,mass,mass0,natom1,natom2,ref1_0,ref2_0,Xdim_file)
        


        
        ntest=1000000

        allocate(new_mass(natom1+natom2,ntest),new_mass_(natom1+natom2,ntest))

        if (inputCoord == "Cartesian")then
            szi = 3*(natom1+natom2)
        else
            szi = XDIM
        endif

        ALLOCATE(internal(szi,ntest))

     


        call Generate_RandomData(internal,szi,rm,new_mass,mass,mass0,natom1,natom2,ref1_0,ref2_0&
                                    ,XDIM,inputCoord,ntest)   

                                   

        allocate(ref1_temp0(natom1*3),ref2_temp0(natom2*3),nmass(natom))

        ! make temporary copies of the input data
        ref1_temp0=ref1_0! reference vector for frag1 (original frame)
        ref2_temp0=ref2_0! reference vector for frag2 (original frame)
    
        ! make sure the original Cartesian frame of reference is at the CM of each fragment 
        if(natom1>1)call rm_cmass(ref1_temp0,mass0(1:natom1),natom1,natom1)
        if(natom2>1)call rm_cmass(ref2_temp0,mass0(natom1+1:natom),natom2,natom2) 


        if (outputCoord == "Cartesian")then
              internal0_length = 3*(natom1+natom2)
          else
              internal0_length = XDIM
        endif

       allocate(internal0_(internal0_length),internal0(internal0_length))

   

        do while (nc < ntest)
       
                 nc=nc+1
                 passingThrough =0d0
                 internal_i  = internal(:,nc)
                 nmass = new_mass(:,nc)
                 call Int_to_Cart(internal_i,szi,nmass,ref1_temp0,ref2_temp0,XDIM,natom1,natom2,inputCoord,cart)

                 call convert_isotopic_coordinates(cart,R_ZYZ,nmass,mass0,natom1,natom2,ref1_temp0,ref2_temp0,XDIM,td,1)
                 
                 call ZYZ_to_OutPut(internal_i,internal0_,internal0_length,R_ZYZ,ref1_temp0,ref2_temp0,XDIM,&
                            natom1,natom2,outputCoord)

                  internal0_ = internal0

                  passingThrough = 0d0

                  if (outputCoord =="Cartesian")then
                  else
                        
                    if (dabs(internal0(2))< threshold  ) then
                    passingThrough = 1d0
                    if (size(internal0)>3)then
                        if (dabs(internal0(3)) < threshold)then
                            passingThrough = 1d0
                        else
                            passingThrough = 0d0
                        end if 
                    end if
                

                endif


                if (passingThrough>0d0  ) then 
=======
        integer,INTENT(IN) :: XDIM,szi
        character(*),INTENT(IN) :: filePath,inputCoord,outputCoord
        real*8, INTENT(IN):: internal(szi)
        real*8,allocatable, INTENT(OUT):: internal0_(:)
        real*8, optional :: testArr_Errors(4)
     

        integer :: natom1,natom2,natom
        real*8, allocatable :: ref1_0(:),ref2_0(:),mass(:),mass0(:),ref1_temp0(:),ref2_temp0(:),cart(:)
        real*8 ::  R_ZYZ(7)
        Integer :: Xdim_file,internal0_length,internalLength
        
     
        integer :: i,k,nc
        real*8, ALLOCATABLE :: internal(:),internal0(:)
        real*8 :: pii,threshold,maxdist,passingThrough
        real*8, allocatable :: CasesRand(:,:)
        real*8 :: Max_test_dist,testArr_Errors(4),err
        integer :: counterCase, ntest
        character(*),INTENT(IN)::filename,inputCoord,outputCoord
        Integer,INTENT(OUT)::test_failed
        Integer,INTENT(IN)::XDIM
        integer,optional:: fileOutputNumber
        character(len=50)::testName
        

        real*8 ::  td(4)


        natom = natom1 + natom2

   

         Call Read_File(filename,mass,mass0,natom1,natom2,ref1_0,ref2_0,Xdim_file)



    

        ntest=10

        call Generate_RandomData(internal,rm,mass,mass0,natom1,natom2,ref1_0,ref2_0,XDIM,inputCoord,ntest)   

        ! ALLOCATE(CasesRand(XDIM,test_number + 2000))
        ! pii=acos(-1d0) 
        ! Max_test_dist =0
        ! err = 10d0**(-8)
        ! counterCase = 0

        ! call random_seed()
        ! call random_number(CasesRand)


        ! allocate(internal(XDIM))!,internal0(XDIM))



        !     threshold = 0.9993d0
        !     nc = 0

        !     counterCase = 0
        !     test_failed = 0
        !     Max_test_dist = 0
            
            

        !     do while (counterCase < test_number)
        !         nc=nc+1
        !         passingThrough =0d0

        !         internal(1)=CasesRand(1,nc)*10d0 + 5d0
        !         internal(2)=(CasesRand(2,nc)*2d0-1d0)*threshold
        !         if (XDIM > 2) then  
        !             if (XDIM == 3) then 
        !                 internal(3)=CasesRand(3,nc)*2d0*pii 
        !             else
        !                 internal(3)=(CasesRand(3,nc)*2d0-1d0)*threshold 
        !             end if
        !         end if
        !         if (XDIM>3)then  
        !             internal(4)=CasesRand(4,nc)*2d0*pii 
        !         end if
        !         if (XDIM>4)then  
        !             internal(5)=CasesRand(5,nc)*2d0*pii 
        !         end if
        !         if (XDIM>5)then  
        !             internal(6)=CasesRand(6,nc)*2d0*pii 
        !         end if


        !         Write(*,*)internal,XDIM,inputCoord,outputCoord,filename
              
                
        !         call Get_ISOTOP_COORDINATES(internal,size(internal),internal0,XDIM,inputCoord,outputCoord,filename,testArr_Errors)

                
        !         passingThrough = 0d0
        !         if (dabs(internal0(2))< threshold  ) then
        !             passingThrough = 1d0
        !             if (size(internal0)>3)then
        !                 if (dabs(internal0(3)) < threshold)then
        !                     passingThrough = 1d0
        !                 else
        !                     passingThrough = 0d0
        !                 end if 
        !             end if
                

        !         endif


        !         if (passingThrough>0d0  ) then 
>>>>>>> 19b107b8a1b92c01564ed2550ae60ef0d3bec0e5

                    
                    
                    
            
<<<<<<< HEAD
                    counterCase = counterCase + 1
                    maxdist = MAXVAL(td);
                    if (maxdist > Max_test_dist ) then 
                        Max_test_dist = maxdist 
                    endif
                
                    if (maxdist > err  ) then 
                        write(fileOutputNumber,*)"--------------------------------------- "
                        write(fileOutputNumber,*)testName, td
                        write(fileOutputNumber,*)"Internal : ", internal
                        write(fileOutputNumber,*)"Internal0 : ",internal0
                        write(fileOutputNumber,*)"----------------------------------------"
                        test_failed = test_failed + 1

                    endif
                endif
                !End Testing Section

                  endif

        end do

        if (rm==-1)then
              massStatement = "-File Masses"
        elseif(rm==0)then
              massStatement = "-Same Masses"
        else
              massStatement = "-Random Masses"
        endif

        write(testName,'(A(I1)A)') "Interatomic Distances(",XDIM,"D)"//massStatement
        testName = trim(testName);
        call Print_Test_Results(testName,Max_test_dist,inputCoord,outputCoord,counterCase,counterCase - test_failed)

        write(fileOutputNumber,*)
        write(fileOutputNumber,*)"TEST NAME: ", testName
        write(fileOutputNumber,*)"System Features : Dimension/ Number of Atoms", XDIM
        write(fileOutputNumber,*)"Max_test_dist : ", Max_test_dist
        write(fileOutputNumber,*)"Input/Output Coordinates: ",inputCoord," / ",outputCoord
        write(fileOutputNumber,*)"Number of Tests", counterCase
        write(fileOutputNumber,*)"Failed Tests: ", test_failed, " out of",ntest 
        write(fileOutputNumber,*)
      
     
        deallocate(ref1_0,ref2_0,new_mass,mass,mass0,ref1_temp0,ref2_temp0,cart,internal0_,nmass)
=======
        !             counterCase = counterCase + 1
        !             maxdist = MAXVAL(testArr_Errors);
        !             if (maxdist > Max_test_dist ) then 
        !                 Max_test_dist = maxdist 
        !             endif
                
        !             if (maxdist > err  ) then 
        !                 write(fileOutputNumber,*)"--------------------------------------- "
        !                 write(fileOutputNumber,*)"Intermolecular Distances Test : ", testArr_Errors
        !                 write(fileOutputNumber,*)"Internal : ", internal
        !                 write(fileOutputNumber,*)"Internal0 : ",internal0
        !                 write(fileOutputNumber,*)"----------------------------------------"
        !                 test_failed = test_failed + 1

        !             endif
        !         endif
        !         !End Testing Section
        ! enddo 

        ! !testName = "Interatomic Distances("XDIM//"D), ifun: "//internalFunction
        ! write(testName,'(A(I1)A)') "Interatomic Distances(",XDIM,"D)"
        ! testName = trim(testName);

        ! call Print_Test_Results(testName,Max_test_dist,inputCoord,outputCoord,counterCase,counterCase - test_failed)
        !     write(fileOutputNumber,*)
        !     write(fileOutputNumber,*)"TEST NAME: ", testName
        !     write(fileOutputNumber,*)"System Features : Dimension/ Number of Atoms", XDIM
        !     write(fileOutputNumber,*)"Max_test_dist : ", Max_test_dist
        !     write(fileOutputNumber,*)"Input/Output Coordinates: ",inputCoord," / ",outputCoord
        !     write(fileOutputNumber,*)"Number of Tests", counterCase
        !     write(fileOutputNumber,*)"Failed Tests: ", test_failed, " out of",test_number 
        !     write(fileOutputNumber,*)
      
>>>>>>> 19b107b8a1b92c01564ed2550ae60ef0d3bec0e5

    END SUBROUTINE Testing_InteratomicDistances_v2




! SUBROUTINE Testing_AxisChanges(filename,internalFunction,test_failed,fileOutputNumber)
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

            
            
            

!            CALL convert_isotopic_coordinates(internal,internal0,mass,mass0,natom1,natom2,ref1,ref2,XDIM,internalFunction,test_dist)



!            passingThrough = 0d0
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
    