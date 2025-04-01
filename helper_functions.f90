module helper_functions
    use mod_types, only : int32, real64
    contains

    SUBROUTINE FileChecking(fileName,num)
      implicit none
      character(*), intent(in)::fileName
      integer(int32), intent(in)::num
      logical :: exist

      inquire(file=fileName, exist=exist)
      if (exist) then
        open(num, file=fileName, status="old", position="append", action="write")
      else
        open(num, file=fileName, status="new", action="write")
      end if

    END SUBROUTINE FileChecking

    SUBROUTINE Read_File(filename,mass,mass0,natom1,natom2,ref1,ref2,xdim)
      implicit none
      integer(int32),INTENT(INOUT) :: natom1,natom2,xdim
      real(real64), allocatable,INTENT(INOUT) :: ref1(:),ref2(:),mass0(:),mass(:)
      integer(int32) :: n_atom,i,k
     
     
      character(*), intent(in)::filename
    
      open(unit=10,file=filename,status='old',action='read')
    
      read(10,*) xdim
      read(10,*) natom1! number of atoms in fragment1
      read(10,*) natom2! number of atoms in fragment2

      n_atom=natom1+natom2 ! total number of atoms
    
    
      allocate(ref1(3*natom1),ref2(3*natom2),mass0(n_atom),mass(n_atom))
    
    
    
      if(natom1>1)then
      do i=1,natom1
          read(10,*)mass(i),mass0(i),(ref1(k),k=i*3-2,i*3)
      enddo
      elseif(natom1==1)then
      read(10,*)mass(1),mass0(1)
      ref1=0.d0
      endif
    
      if(natom2>1)then
      do i=1,natom2
          read(10,*)mass(natom1+i),mass0(natom1+i),(ref2(k),k=i*3-2,i*3)
      enddo
      elseif(natom2==1)then
      read(10,*)mass(natom1+1),mass0(natom1+1)
      ref2=0.d0
      endif
      close(10)
    END SUBROUTINE  Read_File

    SUBROUTINE MolecularDistance_Comapare(arr1,ref1,ref2,N1,N2,internal0,xdim,Distance)
        use math_functions,only: interatomic_distance,molecular_rotation_zyz ,compared_distance
        use coordinateTransf,only:internal0_to_angles_zyz
        implicit none
        integer(int32), intent(in) ::  N1,N2,xdim
        real(real64), intent(in) ::  internal0(xdim),arr1(3,N1+N2),ref1(3,N1),ref2(3,N2)
        real(real64),INTENT(INOut) ::  Distance(2)
        real(real64) ::  d1((N1+N2)*(N1+N2-1)/2),d0((N1+N2)*(N1+N2-1)/2),ref1_temp(3,N1),ref2_temp(3,N2)
        real(real64):: R,rotarr1(3,N1),rotarr2(3,N2),arr0(3,N1+N2),PII,angles(6)
        real(real64):: a1,b1,c1,a2,b2,c2,internal0_temp(xdim)
        
  
        PII=dacos(-1d0) 
        R=internal0(1)
  
        internal0_temp = internal0
        ref1_temp=ref1
        ref2_temp=ref2
  
  
        call internal0_to_angles_zyz(internal0_temp,xdim,angles)
  
  
            
            a1 = angles(1)
            b1 = angles(2)
            c1 = angles(3)
            a2 = angles(4)
            b2 = angles(5)
            c2 = angles(6)
  
            ! write(*,*)"a1,b1,c1 ",a1,b1,c1
            ! write(*,*)"a2,b2,c2 ",a2,b2,c2
  
  
            call molecular_rotation_zyz(ref1_temp,N1,a1,b1,c1,rotarr1)
            call molecular_rotation_zyz(ref2_temp,N2,a2,b2,c2,rotarr2)
        
        
            arr0(1:3,1:N1) = rotarr1
            arr0(1:2,N1+1:N2) = rotarr2(1:2,1:N2)
            arr0(3,N1+1:N2) = rotarr2(3,1:N2)+R
  
  
        call interatomic_distance(arr1,(N1+N2),d1)
        call interatomic_distance(arr0,(N1+N2),d0)
  
    
  
        call compared_distance(d1,d0,(N1+N2)*((N1+N2)-1)/2,Distance)
        
    
    END SUBROUTINE MolecularDistance_Comapare
  
    !input  is in cartesian
    !output is R and the 6 euler angles in ZYZ

    SUBROUTINE convert_isotopic_coordinates(cart,R_ZYZ,mass,mass0,natom1,natom2,ref1_0,ref2_0,xdim,testErr,doTest)
       
        use math_functions,only: find_euler_angles, remove_center_of_mass,vec_to_mat2,mat_to_vec2, rotate_molecule_v2,&
                                 center_of_mass
        use coordinateTransf,only: euler_angles_to_autosurf

        implicit none
  
        integer(int32), intent(in) :: natom1,natom2,xdim
        real(real64), intent(in) :: cart((natom1+natom2)*3),&
                                    mass(natom1+natom2),mass0(natom1+natom2),&
                                    ref1_0(natom1*3),ref2_0(natom2*3)
        real(real64),INTENT(INOUT) :: R_ZYZ(7)
        real(real64),INTENT(INOUT):: testErr(4)
        integer(int32) , optional :: doTest
        
        real(real64) :: ref1_temp0(natom1*3),ref2_temp0(natom2*3),cm(3)
        real(real64) :: ref1(natom1*3),ref2(natom1*3),cart_ref1(3,natom1),cart_mat1(3,natom1)
        real(real64) :: cart_ref2(3,natom2),cart_mat2(3,natom2),cart_frag2(natom2*3)
        real(real64) :: cart_mat(3,natom1+natom2)
        real(real64) :: U_rot(3,3)
        real(real64) :: gamma1,gamma2,beta1,beta2,alpha1,alpha2,theta,phi,cos_theta,sin_theta
        real(real64) :: tci(3,natom1+natom2),tcf(3,natom1+natom2)
        real(real64) :: td(4),distance_arr_test(2),inter0_test(xdim)
        real(real64), parameter :: PII = dacos(-1d0)
        integer(int32) :: n_atom
        
       ! make temporary copies of the input data
       ref1_temp0=ref1_0! reference vector for frag1 (original frame)
       ref2_temp0=ref2_0! reference vector for frag2 (original frame)
   
       ! make sure the original Cartesian frame of reference is at the CM of each fragment 
       if(natom1>1)call remove_center_of_mass(ref1_temp0,mass0(1:natom1),natom1,natom1)
       if(natom2>1)call remove_center_of_mass(ref2_temp0,mass0(natom1+1:n_atom),natom2,natom2)

        n_atom=natom1+natom2
   
        ref1=ref1_temp0
        ref2=ref2_temp0

       
         if(present(doTest))then
            call vec_to_mat2(cart,tci,n_atom) 
          endif

   
  
         !!! 2) transfer the origin of the reference frame to the original CM of frag1 
          call remove_center_of_mass(cart,mass0,natom1+natom2,natom1)
  
          !!! 3) align Z and Z' frames
          ! find new relative position of the original CM of frag2 (with respect to the original CM of frag1)
  
          call center_of_mass(cart,cm,mass0,natom1+natom2,natom2)

  
          ! find spherical coordinates of the c.m. of frag2
          R_ZYZ(1) = dsqrt(cm(1)**2+cm(2)**2+cm(3)**2)! "R"! Norm2
  
   
  
          theta = dacos(cm(3)/R_ZYZ(1))! "TH"
          cos_theta = cm(3)/R_ZYZ(1)
          sin_theta = DSQRT(1d0 - cos_theta**2)
  
          if(dabs(dcos(theta))>1d0-1d-12)then! "PHI"
            !td(5) = -1
            phi=0d0
          else
            phi=datan2(cm(2),cm(1))
          endif
  
  
        ! construct rotation matrix
  
       U_rot=0d0
       U_rot(1,1)=dcos(phi)*cos_theta
       U_rot(2,2)=dcos(phi)
       U_rot(3,3)=cos_theta
       U_rot(1,2)=dsin(phi)*cos_theta
       U_rot(1,3)=-1d0*sin_theta
       U_rot(2,1)=-1d0*dsin(phi)
       U_rot(3,1)=dcos(phi)*sin_theta
       U_rot(3,2)=dsin(phi)*sin_theta
  
  
  
       ! transform vector "cart[1:3*n_atom]" into a matrix: "cart_mat[1:3,1:n_atom]"
       call vec_to_mat2(cart,cart_mat,natom1+natom2)
       ! rotate the initial geometry to align (make "Z" axis contain the c.m. of both frags!)
       call rotate_molecule_v2(natom1+natom2,cart_mat,cart_mat,U_rot)! Now "cart_mat" is in the proper Frame
  
  
       ! transform back matrix "cart_mat[1:3,1:n_atom]" into a vector: "cart[1:3*n_atom]"
       call mat_to_vec2(cart_mat,cart,natom1+natom2)! Now "cart" is in the proper/original Frame: where the 
       ! PES was initially fitted, with its origin at the CM of frag1 and the Z-axis containing both CMs
       
       if(present(doTest))then
              call vec_to_mat2(cart,tcf,n_atom) 
       endif
       
       !!! ----------------------------------------------------------------------------------
       !!! 4) once the geometry is in the proper frame, find the corresponding Euler angles:
       !!! ----------------------------------------------------------------------------------
  
       ! Fragment 1 
       ! -----------
       ! cart_mat1 = new (rotated/aligned) coordinates of frag 1 
       cart_mat1(1:3,1:natom1)=cart_mat(1:3,1:natom1)
       ! transform vector "ref1_temp0[1:3*natom1]" (original reference for frag1) into a matrix: "cart_ref1[1:3,1:natom1]"
       call vec_to_mat2(ref1_temp0,cart_ref1,natom1)
       ! find rotation matrix (U_rot) that transforms (rotates) coordinates "cart_ref2" into "cart_mat2"
  
  
  
       call find_euler_angles(natom1,cart_ref1,cart_mat1,mass0(1:natom1),alpha1, beta1,gamma1)
   
  
  
       ! Fragment 2 
       ! -----------
       cart_frag2=cart(3*natom1+1:3*(natom1+natom2))
       ! set center of mass of frag2 at origin
       call remove_center_of_mass(cart_frag2,mass0(natom1+1:natom1+natom2),natom2,natom2)
       ! transform vector "cart_frag2[1:3*natom2]" (new --rotated/aligned-- coord. of frag 1) into a matrix: "cart_mat2[1:3,1:natom2]"
       call vec_to_mat2(cart_frag2,cart_mat2,natom2)
       ! transform vector "ref2_temp0[1:3*natom2]" (original reference for frag2) into a matrix: "cart_ref2[1:3,1:natom2]"
       call vec_to_mat2(ref2_temp0,cart_ref2,natom2)
       ! find rotation matrix (U_rot2) that transforms (rotates) coordinates "cart_ref2" into "cart_mat2"
  
       call find_euler_angles(natom2,cart_ref2,cart_mat2,mass0(natom1+1:natom1+natom2),alpha2, beta2,gamma2)
  
  
       ! *** internal coordinates ***
       ! ----------------------------
  
       R_ZYZ(2) = alpha1
       R_ZYZ(3) = beta1
       R_ZYZ(4) = gamma1
       R_ZYZ(5) = alpha2
       R_ZYZ(6) = beta2
       R_ZYZ(7) = gamma2

     


      if(present(doTest))then
       
            if (doTest>0)then
              !!! Testing Unit
                
                    ! This function test if the cart coodinates given by int_to_cart in step 1
                    ! has the same interatomic distances than the bimolecular system in the internal coordiantes internal0
                    ! This test calculate an array d1 and d0 of all interatomic distance of all atoms taken in pairs for both 
                    ! cartesian configurations   and print the sum of abs(d1-d0)
              
                    call euler_angles_to_autosurf(xdim,R_ZYZ,inter0_test)
                    call  MolecularDistance_Comapare(tcf,cart_ref1,cart_ref2,natom1,natom2,&
                                                      inter0_test,xdim,distance_arr_test)
                    td(1:2) = distance_arr_test
                    call  MolecularDistance_Comapare(tci,cart_ref1,cart_ref2,natom1,natom2,&
                                                      inter0_test,xdim,distance_arr_test)
                    td(3:4) = distance_arr_test

                    testErr = td
              !!! End Testing Unit
            end if
        endif
     
     
      

  
    END SUBROUTINE convert_isotopic_coordinates

  
    SUBROUTINE get_pes_coordinates(  PES_internal_coordinate,&
                                        isotope_internal_coordinate,&
                                        isotope_coordinate_size,&
                                        xdim,&
                                        isotope_coordinate_format,&
                                        PES_coordinate_format,&
                                        path_to_file,&
                                        errors)
        
        use math_functions,only: remove_center_of_mass
        use coordinateTransf,only: int_to_cart,zyz_to_output
        

        implicit none
        real(real64),allocatable, intent(out):: PES_internal_coordinate(:)
        integer(int32), intent(in) :: xdim
        integer(int32),intent(in) :: isotope_coordinate_size
        real(real64),  intent(in):: isotope_internal_coordinate(isotope_coordinate_size)
        character(*,kind=1), intent(in) :: isotope_coordinate_format
        character(*,kind=1), intent(in) :: PES_coordinate_format
        character(*,kind=1), intent(in) :: path_to_file
        real(real64),dimension(4), optional :: errors ! only for testing
        ! local variables
        integer(int32) :: natom1,natom2,n_atom
        real(real64), allocatable :: ref1_0(:),ref2_0(:),mass(:),mass0(:),ref1_temp0(:),ref2_temp0(:),cart(:)
        real(real64) ::  R_ZYZ(7) ! general coordinate array in EUler-ZYZ
        integer(int32) :: xdim_file,PES_internal_coordinate_size
        
        integer(int32) :: initflag
        save initflag
        data initflag /1/
        save mass,mass0,natom1,natom2,ref1_0,ref2_0,xdim_file
        real(real64) ::  td(4)

        n_atom = natom1 + natom2
        ! Reading Cartesian Molecular frames from file
        if(initflag==1)then! initialize 
         Call Read_File(path_to_file,mass,mass0,natom1,natom2,ref1_0,ref2_0,xdim_file)
         initflag=2  
        endif
      
        allocate(ref1_temp0(natom1*3),ref2_temp0(natom2*3))

        ! make temporary copies of the input data
        ref1_temp0=ref1_0! reference vector for frag1 (original frame)
        ref2_temp0=ref2_0! reference vector for frag2 (original frame)
    
        ! make sure the original Cartesian frame of reference is at the CM of each fragment 
        if(natom1>1)call remove_center_of_mass(ref1_temp0,mass0(1:natom1),natom1,natom1)
        if(natom2>1)call remove_center_of_mass(ref2_temp0,mass0(natom1+1:n_atom),natom2,natom2)  

        if (PES_coordinate_format == "Cartesian")then
            PES_internal_coordinate_size = 3*(natom1+natom2)
          else
            PES_internal_coordinate_size = xdim
        endif

        allocate(PES_internal_coordinate(PES_internal_coordinate_size))

        call int_to_cart(   isotope_internal_coordinate,&
                            isotope_coordinate_size,&
                            mass,ref1_temp0,ref2_temp0,xdim,&
                            natom1,natom2,&
                            isotope_coordinate_format,&
                            cart)

        if (present(errors))then
         call convert_isotopic_coordinates(cart,R_ZYZ,mass,mass0,natom1,natom2,ref1_temp0,ref2_temp0,xdim,td,1)
         errors = td
        else
         call convert_isotopic_coordinates(cart,R_ZYZ,mass,mass0,natom1,natom2,ref1_temp0,ref2_temp0,xdim,td)
        end if 

    
        call zyz_to_output(     isotope_internal_coordinate,&
                                PES_internal_coordinate,&
                                PES_internal_coordinate_size,&
                                R_ZYZ,ref1_temp0,ref2_temp0,xdim,&
                                natom1,natom2,&
                                PES_coordinate_format)


    END SUBROUTINE get_pes_coordinates

end module helper_functions
    