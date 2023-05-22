! In thy function is defined the sys

SUBROUTINE Get_ISOTOP_COORDINATES(internal,internal0,XDIM,ifun,filePath)
  use helperFunc

  IMPLICIT NONE
  integer,INTENT(IN) :: XDIM,ifun
  character(*),INTENT(IN) :: filePath
  real*8, INTENT(IN):: internal(XDIM)
  real*8, INTENT(OUT):: internal0(XDIM)
  integer :: natom1,natom2
  real*8, allocatable :: ref1_0(:),ref2_0(:),mass(:),mass0(:)
  real*8 ::  int0_AS(XDIM),int0_BiSph(XDIM),td(5)
  Integer :: ifun_temp,Xdim_file
  integer :: initflag
  save initflag
  character(len=25) :: i0Type ! it defines what EulerAngle_to_FinalInternal0 has to be defined 
  data initflag /1/
  save mass,mass0,natom1,natom2,ref1_0,ref2_0,Xdim_file
  
  i0Type = "BiSpherical" ! options: "BiSpherical","AutoSurf"
  i0Type = trim(i0Type)

  IF(initflag==1)THEN! initialize 
   Call Read_File(filePath,mass,mass0,natom1,natom2,ref1_0,ref2_0,Xdim_file)
   initflag=2  
  ENDIF



  ifun_temp =ifun


 Call convert_isotopic_coordinates(internal,int0_BiSph,mass,mass0,natom1,natom2,ref1_0,ref2_0,XDIM,ifun_temp,td,-1,i0Type)



 
 internal0 = int0_BiSph
 
 
 


END SUBROUTINE Get_ISOTOP_COORDINATES
    
    
    

    
    
    
