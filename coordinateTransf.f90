module coordinateTransf
  
    contains


    ! convert: Internal coordinates in ZYZ--> BiSpherical for 5D ONLY
    SUBROUTINE ZYZ_to_BiSpherical(internal,internal0_ZYZ,internal0_BiSph)

      !************** Subroutine definition 
          ! subroutine ZYZ --> BiSpherical (convert: Internal coordinates in ZYZ--> BiSpherical for 5D)  
          ! internal0_ZYZ(1) = R
          ! internal0_ZYZ(2) = cos_b1
          ! internal0_ZYZ(3) = cos_b2
          ! internal0_ZYZ(4) = alpha
          ! internal0_ZYZ(5) = gamma

          ! internal0_BiSph(1) = R
          ! internal0_BiSph(2) = theta1
          ! internal0_BiSph(3) = theta2
          ! internal0_BiSph(4) = phi1
          ! internal0_BiSph(5) = phi2
      !********************************************************
          IMPLICIT NONE
          real*8,INTENT(IN) :: internal0_ZYZ(5) ,internal(5)     
          real*8,INTENT(OUT)::internal0_BiSph(5) 
          real*8:: pii,term1,term2,threshold,internal0p_ZYZ(5)  
          
          !write(*,*)"internal0_ZYZ",internal0_ZYZ
          pii=dacos(-1d0)
          threshold = 1d0 - 1d-12!0.9993d0
          internal0p_ZYZ = internal0_ZYZ

          internal0_BiSph(1) =  internal0_ZYZ(1) ! R

          if(internal0p_ZYZ(2)>1d0)then
            internal0p_ZYZ(2) = 1d0
          elseif(internal0p_ZYZ(2)<-1d0)then
            internal0p_ZYZ(2) = -1d0
          endif

          if(internal0p_ZYZ(3)>1d0)then
            internal0p_ZYZ(3) = 1d0
          elseif(internal0p_ZYZ(3)<-1d0)then
            internal0p_ZYZ(3) = -1d0
          endif


     ! Note : There is a discontinuity in Phi =0 for cos_b1 =0 
       ! Taylor serie for cos b1 
       ! Analyse b1 = 0 and b2 = 0


            
  
          if(DABS(internal0p_ZYZ(2))>threshold .and. DABS(internal0p_ZYZ(3))>threshold)then

              if(internal0p_ZYZ(2)>threshold)then
                internal0_BiSph(2) = 0d0
                internal0_BiSph(4) = internal(4)
              elseif( internal0p_ZYZ(2) <-1d0*threshold)then
                internal0_BiSph(2) = pii
                internal0_BiSph(4) = internal(4)
              else
                internal0_BiSph(2) = DACOS(internal0p_ZYZ(2))
                internal0_BiSph(4) = pii- internal0p_ZYZ(5)
              end if
              if(internal0p_ZYZ(3)>threshold)then
                internal0_BiSph(3) = 0d0
                internal0_BiSph(5) = internal(5)
              elseif( internal0p_ZYZ(3) <-1d0*threshold)then
                internal0_BiSph(3) = pii
                internal0_BiSph(5) = internal(5)
              end if

          else

              if(internal0p_ZYZ(2)>threshold)then
                internal0_BiSph(2) = 0d0
                internal0_BiSph(3) = DACOS(internal0p_ZYZ(3))
                internal0_BiSph(4) = internal(4)
                internal0_BiSph(5) = -1d0*(internal0p_ZYZ(4)+internal0p_ZYZ(5))
              elseif( internal0p_ZYZ(2) <-1d0*threshold)then
                internal0_BiSph(2) = pii
                internal0_BiSph(3) = DACOS(-1d0*internal0p_ZYZ(3))
                internal0_BiSph(4) = internal(4)
                internal0_BiSph(5) = internal0p_ZYZ(4)+internal0p_ZYZ(5)-pii
              elseif( internal0p_ZYZ(3) >threshold)then
                internal0_BiSph(2) = DACOS(internal0p_ZYZ(2))
                internal0_BiSph(3) = 0d0
                internal0_BiSph(4) =  pii- internal0p_ZYZ(5)
                internal0_BiSph(5) = internal(5)
              elseif( internal0p_ZYZ(3) <-1d0*threshold)then
                internal0_BiSph(2) = DACOS(internal0p_ZYZ(2))
                internal0_BiSph(3) = pii
                internal0_BiSph(4) = pii- internal0p_ZYZ(5)
                internal0_BiSph(5) = internal(5) 
              else

                  internal0_BiSph(2) = DACOS(internal0p_ZYZ(2))
                  internal0_BiSph(3) = DACOS(internal0p_ZYZ(2)*internal0p_ZYZ(3)+dcos(internal0p_ZYZ(4)) &
                                        *DSQRT((1d0-internal0p_ZYZ(2)**2)*(1d0-internal0p_ZYZ(3)**2)))

                  internal0_BiSph(4) = pii- internal0p_ZYZ(5)
              
                  term1 = DCos(internal0p_ZYZ(5))*DSIN(internal0p_ZYZ(4))*DSQRT(1d0-internal0p_ZYZ(3)**2) &
                          + DSIN(internal0p_ZYZ(5))*(DCOS(internal0p_ZYZ(4))*internal0p_ZYZ(2)*DSQRT(1d0-internal0p_ZYZ(3)**2) &
                          - DSQRT(1d0-internal0p_ZYZ(2)**2)*internal0p_ZYZ(3) )
                          
              
              
                  term2 = DSQRT(1d0-internal0p_ZYZ(3)**2)*(DCOS(internal0p_ZYZ(4))* & 
                          internal0p_ZYZ(2)*DCos(internal0p_ZYZ(5)) - DSIN(internal0p_ZYZ(4))*DSIN(internal0p_ZYZ(5))) &
                          - internal0p_ZYZ(3)*DCos(internal0p_ZYZ(5))*DSQRT(1d0-internal0p_ZYZ(2)**2)
              
                  internal0_BiSph(5) = datan2(-1d0*term1,term2)
              end if
          end if
          
          
       
     

          

          !Output Alpha Check
         
          if (internal0_BiSph(4) > pii)then
            internal0_BiSph(4) = internal0_BiSph(4) -  2d0*pii 
          end if 
          if (internal0_BiSph(4) < -pii)then
            internal0_BiSph(4) = internal0_BiSph(4) + 2d0*pii 
          end if 
         
    END SUBROUTINE ZYZ_to_BiSpherical

    ! transformation from Autosurf internal angles coordinates to ZYZ angles
    SUBROUTINE Internal0_to_anglesZYZ(internal0,XDIM,angles)
        
        IMPLICIT NONE
        integer,INTENT(IN) ::  XDIM
        real*8,INTENT(IN) ::  internal0(XDIM)
        real*8,INTENT(INOut) ::  angles(6)
        
        real*8::pii
  
        pii = dacos(-1d0)
  
        
    
  
        if (XDIM==2)then
            angles(1)=0d0
            angles(2)=dacos(internal0(2))
            angles(3)=0d0
            angles(4)=0d0
            angles(5)=0d0
            angles(6)=0d0
        elseif (XDIM==3)then
            angles(1)=0d0
            angles(2)=dacos(internal0(2))
            angles(3)=pii - internal0(3)
            angles(4)=0d0
            angles(5)=0d0
            angles(6)=0d0
        elseif (XDIM==4)then
            angles(1)=internal0(4)
            angles(2)=dacos(internal0(2))
            angles(3)=0d0
            angles(4)=0d0
            angles(5)=dacos(internal0(3))
            angles(6)=0d0     
        elseif (XDIM==5)then
            angles(1)=internal0(4)
            angles(2)=dacos(internal0(2))
            angles(3)=internal0(5)
            angles(4)=0d0
            angles(5)=dacos(internal0(3))
            angles(6)=0d0 
        elseif (XDIM==6)then
            angles(1)=internal0(4)
            angles(2)=dacos(internal0(2))
            angles(3)=internal0(5)
            angles(4)=0d0
            angles(5)=dacos(internal0(3))
            angles(6)=internal0(6)                             
        endif
  
    
    END SUBROUTINE Internal0_to_anglesZYZ

    subroutine EulerAngles_2_BiSpherical(inter,R,a1,b1,g1,a2,b2,g2,internal0)
            ! *** internal coordinates ***
      ! ----------------------------
      
      !************** Subroutine definition 
          ! subroutine ZYZ --> BiSpherical (convert: Internal coordinates in ZYZ--> BiSpherical for 5D)  

          ! Input : R and Euler Angles in radians

          ! Output : R and cos_theta1,cos_theta2, phi_1,phi_2
          
      !********************************************************
      IMPLICIT NONE
      real*8,INTENT(IN):: R,a1,b1,g1,a2,b2,g2,inter(5)
      real*8,INTENT(OUT):: internal0(5)
      real*8::pii  ,alpha1,beta1,gamma1,alpha2,beta2,gamma2
      real*8:: term1,term2,threshold,alpha

      pii=dacos(-1d0) 
      
      alpha1=a1
      beta1=b1
      gamma1=g1
      alpha2=a2
      beta2=b2
      gamma2=g2

      
      alpha=alpha1-alpha2
      if(alpha>pii)then
        alpha=alpha-2d0*pii
      endif
      if(alpha<-pii)then
        alpha=alpha+2d0*pii
      endif
      
      
      !write(*,*)"internal0_ZYZ",internal0_ZYZ
  
      threshold = 1d-8

      internal0(1) = R
      internal0(2) = b1
      

      If (DABS(b1)<threshold  .or. DABS(b1-pii)< threshold)Then
          
          
          internal0(4) = inter(4)

          If (DABS(b1)<threshold)Then
            internal0(2) = 0d0
            internal0(3) = b2
          else
            internal0(2) = pii
            internal0(3) = pii - b2
          end if


          If (DABS(b2)<threshold .or. DABS(b2-pii)<threshold)Then        
            internal0(5) = inter(5)
          else 
            If (DABS(b1)<threshold)Then
              internal0(5) = -1d0*(alpha + gamma1)
            else
              internal0(5) = alpha + gamma1-pii  
            end if 
          end if

      else
        If (DABS(b1)<2d-6  .or. DABS(b1-pii)< 2d-6)Then
          internal0(4) = inter(4)
        else
          internal0(4) = pii- gamma1
        end if 
        

        
          
          If (DABS(b2)<threshold .or. DABS(b2-pii)<threshold)Then   
            internal0(5) = inter(5)     
            If (DABS(b2)<threshold)Then
              internal0(3) = 0d0
            else
              internal0(3) = pii
            end if
          else 
            internal0(3) = DACOS(DCOS(b1)*DCOS(b2)+dcos(alpha)*DSIN(b1)*DSIN(b2))


            term1 = DCos(gamma1)*DSIN(alpha)*DSIN(b2)+ DSIN(gamma1)*(DCOS(alpha)*Dcos(b1)*DSIN(b2) &
            - DSIN(b1)*DCOS(b2) )
            


            term2 = DSIN(b2)*(DCOS(alpha)*DCos(b1)*DCos(gamma1) - DSIN(alpha)*DSIN(gamma1)) &
                    - DCOS(b2)*DCos(gamma1)*DSIN(b1)

            internal0(5) = datan2(-1d0*term1,term2)
          end if

      end if


      

      !Output Alpha Check
     
      if (internal0(4) > pii)then
        internal0(4) = internal0(4) -  2d0*pii 
      end if 
      if (internal0(4) < -pii)then
        internal0(4) = internal0(4) + 2d0*pii 
      end if 

    
  
      return
    end subroutine EulerAngles_2_BiSpherical

    subroutine EulerAngles_2_Autosurf(XDIM,R,a1,b1,g1,a2,b2,g2,internal0)
      IMPLICIT NONE
      Integer,INTENT(IN)::XDIM
      real*8,INTENT(IN):: R,a1,b1,g1,a2,b2,g2
      real*8,INTENT(OUT):: internal0(XDIM)
      real*8::pii  ,alpha1,beta1,gamma1,alpha2,beta2,gamma2
  
      pii=dacos(-1d0) 
      internal0(1) = R
      alpha1=a1
      beta1=b1
      gamma1=g1
      alpha2=a2
      beta2=b2
      gamma2=g2
      ! *** internal coordinates ***
      ! ----------------------------
      if(XDIM==2)then
        internal0(2)=dcos(beta1)
      elseif(XDIM==3)then
        internal0(2)=dcos(beta1)
        ! ----------------------------
        internal0(3)=pii-gamma1
      elseif(XDIM==4)then
        internal0(2)=dcos(beta1)
        ! ----------------------------
        internal0(3)=dcos(beta2)
        ! ----------------------------
        internal0(4)=alpha1-alpha2
        if(internal0(4)>pii)then
          internal0(4)=internal0(4)-2d0*pii
        endif
        if(internal0(4)<-pii)then
          internal0(4)=internal0(4)+2d0*pii
        endif
      elseif(XDIM==5)then
        internal0(5)=gamma1
        ! ----------------------------
        internal0(2)=dcos(beta1)
        ! ----------------------------
        internal0(3)=dcos(beta2)
        ! ----------------------------
        internal0(4)=alpha1-alpha2
        if(internal0(4)>pii)then
          internal0(4)=internal0(4)-2d0*pii
        endif
        if(internal0(4)<-pii)then
          internal0(4)=internal0(4)+2d0*pii
        endif
      elseif(XDIM==6)then
        internal0(5)=gamma1
        ! ----------------------------
        internal0(6)=gamma2
        ! ----------------------------
        internal0(2)=dcos(beta1)
        ! ----------------------------
        internal0(3)=dcos(beta2)
        ! ----------------------------
        internal0(4)=alpha1-alpha2
        if(internal0(4)>pii)then
          internal0(4)=internal0(4)-2d0*pii
        endif
        if(internal0(4)<-pii)then
          internal0(4)=internal0(4)+2d0*pii
        endif
      endif
  
    
  
      return
    end subroutine EulerAngles_2_Autosurf
  

    SUBROUTINE Cart_to_Euler(cart,natom,ref1_0,angles,doTest)
     
        IMPLICIT NONE

        !cart and ref1_0 must have the same origin
  
        integer,INTENT(IN) :: natom
        real*8,INTENT(IN) ::  ref1_0(natom*3),cart(natom*3)
        real*8 :: ref1_temp0(natom*3),mass0(natom)
        real*8 :: ref1(natom*3),cart_ref1(3,natom)
        real*8 :: cart_mat(3,natom)
        real*8 :: gamma1,beta1,alpha1,pii
        real*8 :: eps,sum_norms_ref0,sum_norms_cart,norms(natom)
        real*8,INTENT(INOUT) :: angles(3)
        integer,INTENT(IN) ::doTest ! doTest = -1 will not do the distance test, any other value will 
        integer::i
  

        pii=dacos(-1d0) 

        ! make temporary copies of the input data
        ref1_temp0=ref1_0! reference vector for frag1 (original frame)

        eps = 1d-8

          !Check the same origin

          sum_norms_ref0 = 0d0
          sum_norms_cart = 0d0
          do i=1,natom
            sum_norms_ref0 = sum_norms_ref0 + Norm2(ref1_0(1+(i-1)*3:i*3))
            sum_norms_cart = sum_norms_cart + Norm2(cart(1+(i-1)*3:i*3))
            norms(i) = Norm2(ref1_0(1+(i-1)*3:i*3))-Norm2(cart(1+(i-1)*3:i*3))

          
          end do



        if (sum(Dabs(norms)) < eps)then

                !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of 
                !!!    reference (origin at the new_CM of frag. 1)
            
                ! subroutine INT_Cart (convert: Internal coordinates --> Cartesian coordinates)  INT_Cart(internal,cart)
                  ref1=ref1_temp0
                
            

            
                ! transform back matrix "cart_mat[1:3,1:natom]" into a vector: "cart[1:3*natom]"
                call mat_to_vec2(cart_mat,cart,natom)! Now "cart" is in the proper/original Frame: where the 
                ! PES was initially fitted, with its origin at the CM of frag1 and the Z-axis containing both CMs
            
            
                !!! ----------------------------------------------------------------------------------
                !!! 4) once the geometry is in the proper frame, find the corresponding Euler angles:
                !!! ----------------------------------------------------------------------------------
            
                ! Fragment 1 
                ! -----------
                ! transform vector "ref1_temp0[1:3*natom1]" (original reference for frag1) into a matrix: "cart_ref1[1:3,1:natom1]"
                call vec_to_mat2(ref1_temp0,cart_ref1,natom)
                ! find rotation matrix (U_rot) that transforms (rotates) coordinates "cart_ref2" into "cart_mat2"
            
            
                mass0(1:natom) =1d0
                call Find_EulerAngles(natom,cart_ref1,cart_mat,mass0,alpha1, beta1,gamma1)
            
              
        else
              gamma1 = 0d0
              beta1  = 0d0
              alpha1 = 0d0
          write(*,*)"The data input does not have the same origin as the reference"
        end if

        angles(1) = alpha1
        angles(2) = beta1
        angles(3) = gamma1


    END SUBROUTINE Cart_to_Euler

    SUBROUTINE Cart_to_Autosurf(systPath,dataPath)

      IMPLICIT NONE
      character(*),INTENT(IN) :: systPath,dataPath
      real*8::pii,cart(18),energies(4),angles1(3),angles2(3)
      integer::natoms,i,Xdim,natom1,natom2
      character(len=1)::Atom_label

      real*8, allocatable :: ref1_0(:),ref2_0(:),mass(:),mass0(:),internal0(:)
      integer :: initflag
      save initflag
      character(len=25) :: i0Type ! it defines Internal0 coodinate system in the output
      data initflag /1/
      save mass,mass0,natom1,natom2,ref1_0,ref2_0,Xdim,i0Type
            
            ! i0Type options: "BiSpherical","Autosurf","Cartesian"
            
          
            IF(initflag==1)THEN! initialize 
            Call Read_File(systPath,mass,mass0,natom1,natom2,ref1_0,ref2_0,Xdim,i0Type)
            initflag=2  
            ENDIF
          
            write(*,*) "Test_Cart_to_Euler"

      ALLOCATE(internal0(Xdim)) 
      
      
      open(unit=100,file=dataPath,status='old',action='read')
      ! open(unit=200,file="coord_H2O_Dymer.txt",status='old',action='write')
      ! open(unit=300,file="coord_H2O_Dymer_filtered.txt",status='old',action='write')
      
      
          do i=1,42508
          
              read(100,*)natoms
              read(100,*)energies(1:4)
              read(100,*)Atom_label,cart(1:3)
              read(100,*)Atom_label,cart(4:6)
              read(100,*)Atom_label,cart(7:9)
              read(100,*)Atom_label,cart(10:12)
              read(100,*)Atom_label,cart(13:15)
              read(100,*)Atom_label,cart(16:18)
              
              call Cart_to_Euler(cart(1:natom1*3),natom1,ref1_0,angles1,0)
              call Cart_to_Euler(cart(1+natom1*3:natom2*3),natom2,ref1_0,angles2,0)
              
          
              
              call EulerAngles_2_Autosurf(Xdim,10,angles1(1),angles1(2),angles1(3)&
                                          ,angles2(1),angles2(2),angles2(3),internal0)
          
              
              ! write(200, *) i , internal0, energies(2)
              ! if (internal0(1)>5) then 
              !  write(300, *) i , internal0, energies(2)
              
              ! endif
          enddo
          
          
          
        close(100)
        !  close(200)
        !  close(300)

    End SUBROUTINE Cart_to_Autosurf


end module coordinateTransf
    