parameter(mnz=200, mnx=600, mnp=4*mnx )

common /geom_arrays/ cord(mnz+1,mnx+1,2),area(4,mnz,mnx),phasez(mnz,mnx), &
    rmarker(mnz,mnx),dvol(4,mnz,mnx)

common /all_arrays/ vel(mnz+1,mnx+1,2),strainr(3,4,mnz,mnx),strain(3,mnz,mnx), &
    force(mnz+1,mnx+1,2),balance(mnz+1,mnx+1,2),rmass(mnz+1,mnx+1),amass(mnz+1,mnx+1), &
    stress0(4,4,mnz,mnx),fmelt(mnz,mnx)

common /rheol_arrays/ aps(mnz,mnx),visn(mnz,mnx),irheol_fl(mnz,mnx),e2sr(mnz,mnx)

common /temp_arrays/ temp(mnz+1,mnx+1),temp0(mnz+1,mnx+1),source(mnz,mnx),shrheat(mnz,mnx)
    
common /boundc/ bc(mnz+1,mnx+1,2),ibc(mnz+1,mnx+1),ncod(mnz+1,mnx+1,2), &
   bcstress((mnz+mnx)*4,3), nopbou((mnz+mnx)*4,4), ncodbou((mnz+mnx)*4,3)

common /particles/ xp(mnp), zp(mnp)

! total of 57 arrays
! approx size = 57*mnz*mnz
