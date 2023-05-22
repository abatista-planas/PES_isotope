SUBROUTINE Int_to_Cart(internal,mass,ref1,ref2,XDIM,natom1,natom2,internalFunction,cart)
    use helperFunc
    IMPLICIT NONE

    integer,INTENT(IN) :: natom1,natom2,XDIM,internalFunction
    real*8,INTENT(IN) :: internal(XDIM),mass(natom1+natom2),ref1(natom1*3),ref2(natom2*3)
    integer :: intfunc
    real*8,INTENT(OUT)  :: cart((natom1+natom2)*3)
    real*8::ref1_temp0(natom1*3),ref2_temp0(natom2*3),cart_temp((natom1+natom2)*3)
    
    intfunc = internalFunction
    !************** Subroutine definition 
        ! Place your custom function here! 
    !********************************************************
    ! make temporary copies of the input data
    ref1_temp0 = ref1! reference vector for frag1 (original frame)
    ref2_temp0 = ref2! reference vector for frag2 (original frame)

    if ( intfunc == 0)then
        call Int_to_Cart_User(internal,mass,cart_temp)
      elseif ( intfunc == 1)then
       ! (for CH3CN-He system: internal coordinates taken as spherical coords.)
       call Int_to_Cart_Spherical(internal,cart_temp,mass,natom1,natom2,ref1_temp0)
      elseif ( intfunc == 2)then
      call Int_to_Cart_ZYZ(internal,XDIM,cart_temp,mass,natom1,natom2, ref1_temp0, ref2_temp0)
    endif

    cart = cart_temp
        
END SUBROUTINE Int_to_Cart



SUBROUTINE Int_to_Cart_ZYZ(internal,XDIM,cart,mass,natom1,natom2,ref1_0,ref2_0)

  !************** Subroutine definition 
      ! subroutine Int_to_Cart_ZYZ (convert: Internal coordinates in ZYZ--> Cartesian coordinates)  
      ! This function considers local axis parallels to ref1_0 and ref2_0 (systems known by Autosurf) 
      ! and global axis with origin in Center Mass of A with z-axis passing through Center Mass of B

      ! * the second angle beta  (rotation around Y axis) is given by cos_beta
      ! * the otehr two, alpha and gamma are given in radians
  !********************************************************
    use helperFunc
    IMPLICIT NONE
    integer,INTENT(IN) :: natom1,natom2,XDIM
    real*8,INTENT(IN) :: internal(XDIM),mass(natom1+natom2),ref1_0(natom1*3),ref2_0(natom2*3)
    real*8,INTENT(OUT)::cart((natom1+natom2)*3)
    integer :: i,natom
    real*8 :: cm(3)
    real*8 :: ref1(natom1*3),rotarr1(3,natom1),ref_mat1(3,natom1),arr1(natom1*3)
    real*8 :: ref2(natom2*3),rotarr2(3,natom2),R_arr(3,natom2),ref_mat2(3,natom2),arr2(natom2*3)
    real*8 :: R,angles(6),internal_temp(XDIM)
    
    internal_temp=internal
    natom=natom1+natom2
    ref1=ref1_0
    ref2=ref2_0

    R = internal_temp(1)

    if (XDIM==2)then
        angles(1)=0d0
        angles(2)=dacos(internal_temp(2))
        angles(3)=0d0
        angles(4)=0d0
        angles(5)=0d0
        angles(6)=0d0
    elseif (XDIM==3)then
        angles(1)=0d0
        angles(2)=dacos(internal_temp(2))
        angles(3)=internal_temp(3)
        angles(4)=0d0
        angles(5)=0d0
        angles(6)=0d0
    elseif (XDIM==4)then
        angles(1)=internal_temp(4)
        angles(2)=dacos(internal_temp(2))
        angles(3)=0d0
        angles(4)=0d0
        angles(5)=dacos(internal_temp(3))
        angles(6)=0d0     
    elseif (XDIM==5)then
        angles(1)=internal_temp(4)
        angles(2)=dacos(internal_temp(2))
        angles(3)=internal_temp(5)
        angles(4)=0d0
        angles(5)=dacos(internal_temp(3))
        angles(6)=0d0 
    elseif (XDIM==6)then
        angles(1)=internal_temp(4)
        angles(2)=dacos(internal_temp(2))
        angles(3)=internal_temp(5)
        angles(4)=0d0
        angles(5)=dacos(internal_temp(3))
        angles(6)=internal_temp(6)                                
    endif

    

     !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of 
     !!!    reference (origin at the new_CM of frag. 1)
     
     
    
     ! set new CM of fragment 1 at the origin
     call rm_cmass(ref1,mass(1:natom1),natom1,natom1) 
     call rm_cmass(ref2,mass(natom1+1:natom1+natom2),natom2,natom2) 



    call vec_to_mat2(ref1,ref_mat1,natom1)
    call vec_to_mat2(ref2,ref_mat2,natom2)

    call MolecularRotation_ZYZ(ref_mat1,natom1,angles(1),angles(2),angles(3),rotarr1)
    call MolecularRotation_ZYZ(ref_mat2,natom2,angles(4),angles(5),angles(6),rotarr2)

    

    call mat_to_vec2(rotarr1,arr1,natom1)
    R_arr(1:2,:) = rotarr2(1:2,:)
    R_arr(3,:) = rotarr2(3,:)+R
    call mat_to_vec2(R_arr,arr2,natom2)
    

    do i=1,3*natom1
        cart(i)=arr1(i)! Cartesian coordinates for the atoms in fragment 1
    enddo
    do i=1,3*natom2
        cart(3*natom1+i)=arr2(i)! Cartesian coordinates for the atoms in fragment 2
    enddo

  
END SUBROUTINE Int_to_Cart_ZYZ



SUBROUTINE Int_to_Cart_Spherical(internal,cart,mass,natom1,natom2,ref1_0)
    use helperFunc
    IMPLICIT NONE
    integer,INTENT(IN) :: natom1,natom2
    real*8,INTENT(IN) :: internal(3),mass(natom1+natom2),ref1_0(natom1*3)

    integer :: i,natom
    real*8 :: cm(3)
    real*8 :: ref1(natom1*3)
    real*8 :: ref2(natom2*3)
    real*8 :: cart((natom1+natom2)*3)
    real*8 :: sin_theta


    natom=natom1+natom2
    ref1=ref1_0

    !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of 
    !!!    reference (origin at the new_CM of frag. 1)
    
    ! subroutine INT_Cart (convert: Internal coordinates --> Cartesian coordinates)  INT_Cart(internal,cart)
    ! (for CH3CN-He system: internal coordinates taken as spherical coords.)

    ! set new CM of fragment 1 at the origin
    call rm_cmass(ref1,mass(1:natom1),natom1,natom1) 
    do i=1,3*natom1
    cart(i)=ref1(i)! Cartesian coordinates for the atoms in fragment 1
    enddo
    sin_theta=dsqrt(1.d0-internal(2)**2)
    cm(1)=internal(1)*sin_theta*dcos(internal(3))
    cm(2)=internal(1)*sin_theta*dsin(internal(3))
    cm(3)=internal(1)*internal(2)
    ref2=cm
    do i=1,3
    cart(3*natom1+i)=ref2(i)! Cartesian coordinates for the atoms in fragment 2
    enddo




END SUBROUTINE Int_to_Cart_Spherical 