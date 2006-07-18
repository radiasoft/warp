!     Last change:  JLV   3 Jun 2004    1:09 pm
!    version numero 2   :02/10/2000

subroutine emi2d_init()

USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain
implicit none


! iboundary =1 systeme isole
! iboundary =2 systeme refroidi sur 3 faces
! iboundary =3 systeme avec compensation du courant sur 3 faces
! iboundary =4 systeme avec refroidissement des particules au dela d un rayon
! ou d'une certaine zone x,y

call constantes

nplan=max(nplani,4*nplane)

nxmain=nx
nymain=ny
nplanmain=nplan
call gchange('EMImain',0)

if(cyl) then
call setrep_cyl(frho_cyl)
elseif(trapeze) then
call setrep_trap
else
call setrep(frho)
endif

nspectr_max=(lx_spectr+1)*(ly_spectr+1)*(nenerg+1)

call setpot

if(.not.esirkepov) then
  call deprho()
endif

end subroutine emi2d_init

subroutine deprho()

USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain

implicit none
rhoaux=0. ;	field%rhoi=0.
call depose_ion_rho(rhoaux,field,.false.)

do i=1,nplani
  field%rhoi(:,:)=field%rhoi(:,:)+rhoaux(:,:,i)
enddo

END subroutine deprho

subroutine propalaser()

USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain

implicit none

WRITE(0,*) 'Enter propalaser'

itime=0 ;	time=0.

if(ylbound==periodic) then
  field%Ex(0:field%nx+1,0)          = field%Ex(0:field%nx+1,field%ny+1)
  field%Ex(0:field%nx+1,field%ny+2) = field%Ex(0:field%nx+1,1)
end if
if(xlbound==periodic) then
  field%Ey(0,0:field%ny+1)          = field%Ey(field%nx+1,0:field%ny+1)
  field%Ey(field%nx+2,0:field%ny+1) = field%Ey(1,0:field%ny+1)
end if

if(cyl) then
 write(nprt,*) 'avant tpropag xcentre,rayon_ext = ',xcentre,rayon_ext
tpropag=(xcentre-rayon_ext)*dx-deupi
elseif(trapeze) then
tpropag=xrh1_trap-deupi
else
if(abs(angle) < 1.e-3) then
if(rho2.eq. 0.) then
tpropag=xrh2-deupi
else
tpropag=xrh1-deupi
end if
else
tpropag=xrhmn_angle*dx
end if
end if

npropag=tpropag/dt
write(nprt,*) 'npropag,tpropag = ',npropag,tpropag

do i=1,npropag
call maxwell(base,time)
itime=itime+1
time=itime*dt
enddo

IF(l_copyfields) then
  field%excopy=field%ex; field%eycopy=field%ey; field%bzcopy=field%bz
END if

WRITE(0,*) 'Exit propalaser'
end subroutine propalaser

subroutine emi2d_init2()

USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain

implicit none
integer(ISZ):: istepe

WRITE(0,*) 'Enter emi2d_init2'

if(ylbound==periodic) then
  field%Ex(0:field%nx+1,0)          = field%Ex(0:field%nx+1,field%ny+1)
  field%Ex(0:field%nx+1,field%ny+2) = field%Ex(0:field%nx+1,1)
end if
if(xlbound==periodic) then
  field%Ey(0,0:field%ny+1)          = field%Ey(field%nx+1,0:field%ny+1)
  field%Ey(field%nx+2,0:field%ny+1) = field%Ey(1,0:field%ny+1)
end if

call griuni(base)
if (.not. l_nodeposition) then
  rhoauxfmpolesi=0.
  rhoauxfmpolese=0.
end if

do istepe=1,nspes2

  IF(l_noacceleration) then
    call move_posonly()
  else
   IF(l_residue .and. type_residue==2) then
    call move_residue(field%ex,field%ey,field%bz,exfsum,eyfsum,bzfsum,)
   else
    call move(field%ex,field%ey,field%bz,exfsum,eyfsum,bzfsum,transition_zone)
   END if
  END if

  call exeybzall()

  itime=itime+1 ;		time=itime*dt

enddo

 
field%exi=0. ;	field%eyi=0. ;	field%bzi=0.

write(nentet) field%nx,field%ny,xlong,ylong,real(dx,4),real(dy,4),real(dt,4),real(misme,4),	&
              real(vte,4),real(tiste,4),real(dnc,4),	 			&
              real(xgauss,4),real(ampfac,4),real(teta,4),nbpart, 		&
              real(xrh1,4),real(xrh2,4),real(xrh3,4),real(xrh4,4),real(xrh5,4),	&
              real(rho1,4),real(rho2,4),real(rho3,4),real(rho4,4),real(rho5,4),	&
              real(yrh1,4),real(yrh2,4),					&
              nxcham,nycham,ixcham,jycham
close(nentet)
WRITE(0,*) 'Exit emi2d_init2'

end subroutine emi2d_init2

subroutine move_window_particles()
use mod_ion
implicit none

if (l_elaser_out_plane) then
  ion%ryi(idebi:ifini)=ion%ryi(idebi:ifini)-1.
else
  ion%rxi(idebi:ifini)=ion%rxi(idebi:ifini)-1.
end if

end subroutine move_window_particles

subroutine exeybzall()

USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain

implicit none

logical :: cond

  IF(l_verbose) WRITE(0,*) 'call exeybzall'

  call exeybz(time)
  if(l_moving_window .and. (mod(itime,ndelta_t)==0)) then
     cond = .true.
    if(cond .and. time>tmin_moving_main_window) then
      call move_window_particles()
      if (l_add_particles) then
        if (l_elaser_out_plane) then
          call add_particles(npcell_max,xleft,xright,ytop-2.,ytop-1.) 
        else
          call add_particles(npcell_max,xright-2.,xright-1.,ybot,ytop) 
        end if
      end if
      call boundary_ions(isort_left,isort_right,isort_bot,isort_top)
    end if
  end if

END subroutine exeybzall

subroutine dep_ion()

USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain

implicit none

    IF(l_verbose) WRITE(0,*) 'call dep_ion'

! we assume nplane = 1
   call depose_ion_esirkepov_linear_serial(rhoaux,rhoauxf)
   call depose_ion_jz (rhoaux,rhoauxf)
   field%rhojxjy(:,:,2:4)=field%rhojxjy(:,:,2:4)+rhoaux(:,:,2:4)

end subroutine dep_ion

subroutine move_ions()
USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain
implicit none

  IF(l_verbose) WRITE(0,*) 'call movi'
  IF(l_noacceleration) then
    call movi_posonly()
  else
    IF(l_residue .and. type_residue==2) then
      call movi_residue(field%exi,field%eyi,field%bzi,exfsumi,eyfsumi,bzfsumi)
    else
      if (nvect==1) then
        call movi_linear_serial_3d(field%exi,field%eyi,field%ezi,field%bxi,field%byi,field%bzi, &
                                   exfsumi,eyfsumi,ezfsumi,bxfsumi,byfsumi,bzfsumi,transition_zone)
      else
        call movi(field%exi,field%eyi,field%bzi,exfsumi,eyfsumi,bzfsumi,transition_zone)
      end if
    END if
  END if

  if(iboundary.eq.1 .or. iboundary.eq.4) then
    IF(l_verbose) WRITE(0,*) 'call boundary_ions'
    call boundary_ions(isort_left,isort_right,isort_bot,isort_top)
  endif

  return
end subroutine move_ions


subroutine zero_ion_fields()
USE mod_field
USE EMImain
implicit none

    field%exi=0.  ;  field%eyi=0.  ;  field%ezi=0.
    field%bxi=0.  ;  field%byi=0.  ;  field%bzi=0.

return
end subroutine zero_ion_fields

subroutine emi2d_step()

USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain

implicit none

INTEGER :: it

IF(l_verbose) WRITE(0,*) 'enter emi2d_step'
do it=1, nsteps			! on avance les ions

  call move_ions()

  if (.not. l_nodeposition) then
    rhoaux=0.
    field%rhojxjy=0.

    call dep_ion_part2()
    call dep_ion()
    call zero_ion_fields()
  end if
  
    call exeybzall()

    itime=itime+1
    time=itime*dt

    IF(l_savefields) then
      IF(l_verbose) WRITE(0,*) 'call save_fields'
      call save_fields(itime,isave_freq)
      IF(l_verbose) WRITE(0,*) 'done'
    END if

enddo

IF(l_verbose) WRITE(0,*) 'exit emi2d_step'
return
!*****************************************************************

end subroutine emi2d_step

subroutine dep_rho()

USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain

implicit none

rhoaux=0. ;	field%rhoi=0.
call depose_ion_rho(rhoaux,field,.false.)

do i=1,nplani
  field%rhoi(:,:)=field%rhoi(:,:)+rhoaux(:,:,i)
enddo

 do i=1,nplani
  field%rhojxjy(:,:,1)=field%rhojxjy(:,:,1)+rhoaux(:,:,i)
 enddo
  field%rhojxjy(:,:,1)=anormr*(field%rhoi+field%rhojxjy(:,:,1))

if(ksmth) then
  call smooth(field,field%rhojxjy(:,:,1),field%nx,field%ny)
endif

end subroutine dep_rho

subroutine check_divemrho()

USE mod_system
USE mod_ion
USE mod_field
USE mod_bnd
use EMImain

implicit none

call dep_rho()

! calcul de dive=rho

call grimaxall()
do k=2,field%ny
do j=2,field%nx
field%rhojxjy(j,k,1)=(field%Ex(j,k) - field%Ex(j-1,k))*field%dxi &
               +(field%Ey(j,k) - field%Ey(j,k-1))*field%dyi &
                          - field%rhojxjy(j,k,1)
enddo
enddo
write(nprt,*) ' *** ecart a la divergence max : ',maxval(field%rhojxjy(ixleft+2:ixright-2,iybot+2:iytop-2,1))
write(nprt,*) ' *** ecart a la divergence min : ',minval(field%rhojxjy(ixleft+2:ixright-2,iybot+2:iytop-2,1))
write(nprt,*) ' '

call griuniall()

END subroutine check_divemrho

subroutine grimaxall()
USE mod_field
call grimax(base)
end subroutine grimaxall

subroutine griuniall()
USE mod_field
call griuni(base)
end subroutine griuniall

      subroutine  isort1 (ia,la)
        
!     function   : isort1 -sort arrays by algebraic value
!     parameters : ia     - on input, contains the array to be sorted
!                  la      - input variable containing the number of
!                            elements in the array to be sorted
      integer t,tt,ia,la,iu,il,m,i,j,k,ij,l
      real(8) :: r
!      dimension ia(1),iu(21),il(21)
      dimension ia(la),iu(21),il(21)

      m=1
      i=1
      j=la
      r=.375

   10 continue
      if (i .eq. j) go to 55
   15 continue
      if (r .gt. .5898437) go to 20
      r=r+3.90625e-2
      go to 25
   20 continue
      r=r-.21875
   25 continue
      k=i
!                                  select a central element of the
!                                  array and save yt in location t
      ij=i+(j-i)*r
      t=ia(ij)
!                                  if first element of array is greater
!                                  than t, interchange with t
      if (ia(i) .le. t) go to 30
      ia(ij)=ia(i)
      ia(i)=t
      t=ia(ij)
   30 continue
      l=j
!                                  if last element of array is less than
!                                  t, interchange with t
      if (ia(j) .ge. t) go to 40
      ia(ij)=ia(j)
      ia(j)=t
      t=ia(ij)
!                                   if first element of array is greater
!                                  than t, interchange with t
      if (ia(i) .le. t) go to 40
      ia(ij)=ia(i)
      ia(i)=t
      t=ia(ij)
      go to 40
   35 continue

      tt=ia(l)
      ia(l)=ia(k)
      ia(k)=tt
!                                  find an element in the second half of
!                                  the array which is smaller than t
   40 continue
      l=l-1
      if (ia(l) .gt. t) go to 40
!                                  find an element in the first half of
!                                  the array which is greater than t
   45 continue
      k=k+1
      if (ia(k) .lt. t) go to 45
!                                  interchange these elements
      if (k .le. l) go to 35
!                                  save upper and lower subscripts of
!                                  the array yet to be sorted
      if (l-i .le. j-k) go to 50
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 60
   50 continue
      il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 60
!                                  begin again on another portion of
!                                  the unsorted array
   55 continue
      m=m-1

      if (m .eq. 0) return 


      i=il(m)
      j=iu(m)
   60 continue

      if (j-i .ge. 11) go to 25
      if (i .eq. 1) go to 10
      i=i-1
   65 continue

      i=i+1
      if (i .eq. j) go to 55
      t=ia(i+1)
      if (ia(i) .le. t) go to 65
      k=i

   70 continue
      ia(k+1)=ia(k)
      k=k-1
      if (t .lt. ia(k)) go to 70
      ia(k+1)=t
      go to 65
      end subroutine  isort1

       subroutine zbrent2(sol,func,x1,x2,x3,tol,ier)

      integer :: itmax=100
      integer :: ier,iter
      real(KIND=8) :: eps=1.e-11
      real(kind=8) :: x1,x2,x3,tol,sol
      real(kind=8) :: a,b,c,d,e,xxx,fa,fb,fc
      real(kind=8) ::  tol1,xm,p,q,r,s
      real(8),external :: func


!     modif jca a zbrent pour accepter une fonction a 2 arguments
 
      a=x1
      b=x2
      xxx=x3
  
      fa=func(a,xxx)
      fb=func(b,xxx)

      if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.)) then
      write(6,*) ' zbrent2 probleme dans l intervalle choisi ',	&
                   ' x1,x2,x3,fa,fb  : ',x1,x2,x3,fa,fb
      sol=0
      ier=1
      return
      endif
!
      c=b
      fc=fb

      do 11 iter=1,itmax
      if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
      c=a
      fc=fa
      d=b-a
      e=d
      endif

      if(abs(fc).lt.abs(fb)) then
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
      endif

      tol1=2.*eps*abs(b)+0.5*tol
      xm=.5*(c-b)
      if(abs(xm).le.tol1 .or. fb.eq.0.)then
      sol=b
      ier=0
      return
      endif

      if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
      s=fb/fa

      if(a.eq.c) then
      p=2.*xm*s
      q=1.-s
      else
      q=fa/fc
      r=fb/fc
      p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
      q=(q-1.)*(r-1.)*(s-1.)
      endif

      if(p.gt.0.) q=-q
      p=abs(p)
      if(2*p .lt. min(3*xm*q-abs(tol1*q),abs(e*q))) then
      e=d
      d=p/q
      else
      d=xm
      e=d
      endif

      else
      d=xm
      e=d
      endif
      a=b
      fa=fb
      if(abs(d) .gt. tol1) then
      b=b+d
      else
      b=b+sign(tol1,xm)
      endif
      fb=func(b,xxx)
11    continue

      write(6,*) 'zbrent2 nombre maximum d iterations depasse'
      ier=2
      sol=b
      return
      end subroutine zbrent2
      double precision function erf(z)
      REAL(8) :: z
      WRITE(0,*) "La fonction erf n'est pas definie!"
      erf = 0.
      stop
      end
      double precision function func(z,zmoy)

!     le zero de func correspond a la valeur du decalage d'une maxwellienne 
!     telle que sa moyenne de zero a l'infini soit egal a vmoy
!     doit etre appellee avec zmoy positif
!     la solution aussi
      real(kind=8) :: z,zmoy
      real(8), external:: erf
      func=zmoy-z-0.5641895*exp(-z**2)/(1+erf(z))
      return
      end

subroutine getex(it)
USE mod_system
USE mod_field
use EMImain
implicit none
INTEGER, INTENT(IN) :: it
CHARACTER(6) :: nnum
   write (nnum, '(i6)') 100000+it
   open(10, file = trim(adjustl(datapath))//'ex_'//nnum//'.dat', status = 'unknown')
   READ(10,*) field%ex
   CLOSE(10)
end subroutine getex
subroutine getey(it)
USE mod_system
USE mod_field
use EMImain
implicit none
INTEGER, INTENT(IN) :: it
CHARACTER(6) :: nnum
   write (nnum, '(i6)') 100000+it
   open(10, file = trim(adjustl(datapath))//'ey_'//nnum//'.dat', status = 'unknown')
   READ(10,*) field%ey
   CLOSE(10)
end subroutine getey
subroutine getbz(it)
USE mod_system
USE mod_field
use EMImain
implicit none
INTEGER, INTENT(IN) :: it
CHARACTER(6) :: nnum
   write (nnum, '(i6)') 100000+it
   open(10, file = trim(adjustl(datapath))//'bz_'//nnum//'.dat', status = 'unknown')
   READ(10,*) field%bz
   CLOSE(10)
end subroutine getbz
subroutine getrhojxjy(it)
USE mod_system
USE mod_field
use EMImain
implicit none
INTEGER, INTENT(IN) :: it
CHARACTER(6) :: nnum
   write (nnum, '(i6)') 100000+it
   open(10, file = trim(adjustl(datapath))//'rhojxjy_'//nnum//'.dat', status = 'unknown')
   READ(10,*) field%rhojxjy
   CLOSE(10)
end subroutine getrhojxjy



