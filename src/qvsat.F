REAL FUNCTION POLYSVP_func (T,TYPE)

  !-------------------------------------------
  
  !  COMPUTE SATURATION VAPOR PRESSURE
  
  !  POLYSVP RETURNED IN UNITS OF PA.
  !  T IS INPUT IN UNITS OF K.
  !  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
  
  ! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM MARAT KHROUTDINOV
  
  IMPLICIT NONE
  
  REAL DUM
  REAL T
  INTEGER TYPE
  ! ice
  real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i 
  data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
       6.11147274, 0.503160820, 0.188439774e-1, &
       0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
       0.385852041e-9, 0.146898966e-11, 0.252751365e-14/	
  
  ! liquid
  real a0,a1,a2,a3,a4,a5,a6,a7,a8 
  data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
       6.105851, 0.4440316, 0.1430341e-1, &
       0.2641412e-3, 0.2995057e-5, 0.2031998e-7, &
       0.6936113e-10, 0.2564861e-13,-0.3704404e-15/
  real dt
  
  ! ICE
  
  IF (TYPE.EQ.1) THEN
     
     !         POLYSVP = 10.**(-9.09718*(273.16/T-1.)-3.56654*                &
     !          LOG10(273.16/T)+0.876793*(1.-T/273.16)+						&
     !          LOG10(6.1071))*100.
     
     
     dt = max(-80.,t-273.16)
     polysvp_func = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt))))))) 
     polysvp_func = polysvp_func*100.
     
  END IF
  
  ! LIQUID
  
  IF (TYPE.EQ.0) THEN
     
     dt = max(-80.,t-273.16)
     polysvp_func = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
     polysvp_func = polysvp_func*100.
     
     !         POLYSVP = 10.**(-7.90298*(373.16/T-1.)+                        &
     !             5.02808*LOG10(373.16/T)-									&
     !             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+				&
     !             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+				&
     !             LOG10(1013.246))*100.
     
  END IF
  
  
END FUNCTION POLYSVP_func

subroutine latent_heat(xxlv,xxls,cpm,t3d,qv3d)
  implicit none
  real,intent(out) :: xxlv,xxls,cpm
  real :: t3d,qv3d,cp,r,rv

  R = 287.15
  RV = 461.5
  CP = 1005.
  
! LATENT HEAT OF VAPORATION
  
  XXLV = 3.1484E6-2370.*T3D

! LATENT HEAT OF SUBLIMATION

  XXLS = 3.15E6-2370.*T3D+0.3337E6
  
  CPM = CP*(1.+0.887*QV3D)
end subroutine latent_heat
            
subroutine saturation_new(t3d,pres,qvs,qvi)
  implicit none
  real :: evs, eis, polysvp_func, t3d, pres
  real, intent(out) :: qvs, qvi
  EVS = POLYSVP_func(T3D,0)   ! PA
  EIS = POLYSVP_func(T3D,1)   ! PA

  IF (EIS.GT.EVS) EIS = EVS
  QVS = .622*EVS/(PRES-EVS)
  QVI = .622*EIS/(PRES-EIS)
!  print*, t3d,pres,evs,qvs    
end subroutine saturation_new
