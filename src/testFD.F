subroutine fdeep_mp_model(v1,v2,v3,v4,o1,o2)
 real :: v1,v2,v3,v4
 real, intent(out) :: o1,o2

 call mp_model(v1,v2,v3,v4,o1,o2)

end subroutine fdeep_mp_model
