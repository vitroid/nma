!******************************************************************************
!*                       Code generated with sympy 1.3                        *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                         This file is part of 'NMA'                         *
!******************************************************************************

REAL*8 function dxixi(phi_r, phi_rr, r, r_x)
implicit none
REAL*8, intent(in) :: phi_r
REAL*8, intent(in) :: phi_rr
REAL*8, intent(in) :: r
REAL*8, intent(in) :: r_x

dxixi = phi_r*(1 - r_x**2/r**2)/r + phi_rr*r_x**2/r**2

end function

REAL*8 function dxiyi(phi_r, phi_rr, r, r_x, r_y)
implicit none
REAL*8, intent(in) :: phi_r
REAL*8, intent(in) :: phi_rr
REAL*8, intent(in) :: r
REAL*8, intent(in) :: r_x
REAL*8, intent(in) :: r_y

dxiyi = -phi_r*r_x*r_y/r**3 + phi_rr*r_x*r_y/r**2

end function

REAL*8 function dxiyk(phi_r, phi_rr, r, r_x, r_y)
implicit none
REAL*8, intent(in) :: phi_r
REAL*8, intent(in) :: phi_rr
REAL*8, intent(in) :: r
REAL*8, intent(in) :: r_x
REAL*8, intent(in) :: r_y

dxiyk = phi_r*r_x*r_y/r**3 - phi_rr*r_x*r_y/r**2

end function

REAL*8 function daiai(phi_r, phi_rr, r, r_x, r_y, r_z, s_aax, s_aay, s_aaz, &
      s_ax, s_ay, s_az)
implicit none
REAL*8, intent(in) :: phi_r
REAL*8, intent(in) :: phi_rr
REAL*8, intent(in) :: r
REAL*8, intent(in) :: r_x
REAL*8, intent(in) :: r_y
REAL*8, intent(in) :: r_z
REAL*8, intent(in) :: s_aax
REAL*8, intent(in) :: s_aay
REAL*8, intent(in) :: s_aaz
REAL*8, intent(in) :: s_ax
REAL*8, intent(in) :: s_ay
REAL*8, intent(in) :: s_az

daiai = (phi_r*(r**2*(r_x*s_aax + r_y*s_aay + r_z*s_aaz + s_ax**2 + s_ay &
      **2 + s_az**2) - (r_x*s_ax + r_y*s_ay + r_z*s_az)**2) + phi_rr*r* &
      (r_x*s_ax + r_y*s_ay + r_z*s_az)**2)/r**3

end function

REAL*8 function daibi(phi_r, phi_rr, r, r_x, r_y, r_z, s_abx, s_aby, s_abz, &
      s_ax, s_ay, s_az, s_bx, s_by, s_bz)
implicit none
REAL*8, intent(in) :: phi_r
REAL*8, intent(in) :: phi_rr
REAL*8, intent(in) :: r
REAL*8, intent(in) :: r_x
REAL*8, intent(in) :: r_y
REAL*8, intent(in) :: r_z
REAL*8, intent(in) :: s_abx
REAL*8, intent(in) :: s_aby
REAL*8, intent(in) :: s_abz
REAL*8, intent(in) :: s_ax
REAL*8, intent(in) :: s_ay
REAL*8, intent(in) :: s_az
REAL*8, intent(in) :: s_bx
REAL*8, intent(in) :: s_by
REAL*8, intent(in) :: s_bz

daibi = (phi_r*(r**2*(r_x*s_abx + r_y*s_aby + r_z*s_abz + s_ax*s_bx + &
      s_ay*s_by + s_az*s_bz) - (r_x*s_ax + r_y*s_ay + r_z*s_az)*(r_x* &
      s_bx + r_y*s_by + r_z*s_bz)) + phi_rr*r*(r_x*s_ax + r_y*s_ay + &
      r_z*s_az)*(r_x*s_bx + r_y*s_by + r_z*s_bz))/r**3

end function

REAL*8 function daibk(phi_r, phi_rr, r, r_x, r_y, r_z, s_ax, s_ay, s_az, &
      t_bx, t_by, t_bz)
implicit none
REAL*8, intent(in) :: phi_r
REAL*8, intent(in) :: phi_rr
REAL*8, intent(in) :: r
REAL*8, intent(in) :: r_x
REAL*8, intent(in) :: r_y
REAL*8, intent(in) :: r_z
REAL*8, intent(in) :: s_ax
REAL*8, intent(in) :: s_ay
REAL*8, intent(in) :: s_az
REAL*8, intent(in) :: t_bx
REAL*8, intent(in) :: t_by
REAL*8, intent(in) :: t_bz

daibk = -(phi_r*(r**2*(s_ax*t_bx + s_ay*t_by + s_az*t_bz) - (r_x*s_ax + &
      r_y*s_ay + r_z*s_az)*(r_x*t_bx + r_y*t_by + r_z*t_bz)) + phi_rr*r &
      *(r_x*s_ax + r_y*s_ay + r_z*s_az)*(r_x*t_bx + r_y*t_by + r_z*t_bz &
      ))/r**3

end function

REAL*8 function dxiai(phi_r, phi_rr, r, r_x, r_y, r_z, s_ax, s_ay, s_az)
implicit none
REAL*8, intent(in) :: phi_r
REAL*8, intent(in) :: phi_rr
REAL*8, intent(in) :: r
REAL*8, intent(in) :: r_x
REAL*8, intent(in) :: r_y
REAL*8, intent(in) :: r_z
REAL*8, intent(in) :: s_ax
REAL*8, intent(in) :: s_ay
REAL*8, intent(in) :: s_az

dxiai = (phi_r*(r**2*s_ax - r_x*(r_x*s_ax + r_y*s_ay + r_z*s_az)) + &
      phi_rr*r*r_x*(r_x*s_ax + r_y*s_ay + r_z*s_az))/r**3

end function

REAL*8 function dxiak(phi_r, phi_rr, r, r_x, r_y, r_z, t_ax, t_ay, t_az)
implicit none
REAL*8, intent(in) :: phi_r
REAL*8, intent(in) :: phi_rr
REAL*8, intent(in) :: r
REAL*8, intent(in) :: r_x
REAL*8, intent(in) :: r_y
REAL*8, intent(in) :: r_z
REAL*8, intent(in) :: t_ax
REAL*8, intent(in) :: t_ay
REAL*8, intent(in) :: t_az

dxiak = -(phi_r*(r**2*t_ax - r_x*(r_x*t_ax + r_y*t_ay + r_z*t_az)) + &
      phi_rr*r*r_x*(r_x*t_ax + r_y*t_ay + r_z*t_az))/r**3

end function
