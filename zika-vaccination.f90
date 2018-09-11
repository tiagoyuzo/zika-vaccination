!========================================================================
!
      module mgeometry
!
      implicit none
!
      integer, parameter :: nsd =2
      integer, parameter :: nf  =1 ! numero de fontes
!
      integer :: nnodes
      integer :: nelem
      integer :: nedges
      integer :: ndof
      integer :: nmax
      integer :: npos
!
      integer :: nx = 16!2**2
      integer :: ny = 16!2**2
      integer :: nt = 560!252*2!2**8
      integer :: expnum = 1
      integer, dimension(9) :: ind
!
      integer, parameter :: nen   = 3 ! number of element nodes (quadrilateral element)
      integer, parameter :: nee   = 4 ! number of element edges (quadrilateral element)
!
      real(8) :: dt, tt, hx, hy, nt1, nt2
      integer :: it, itt, iswp
!
      real(8) :: xl = 1.d0
      real(8) :: yl = 0.6d0!1.d0
      real(8) :: tf = 140.d0!140.d0!42.d0!
!
!     types (structures)
!
!     edge
!
      type :: tedge
!
      integer, dimension(4) :: n
      integer, dimension(2) :: e
!
      end type tedge
!
!     quadrilateral
!
      type :: tquad
!
      integer, dimension(nen)   :: n
      integer, dimension(nee)   :: e
!
      real(8), dimension(nee)   :: ngb ! neighbors
!
      integer, dimension(nen)   :: id
!
      end type tquad
!
!     global arrays
!
!     coordinates
!
      real(8), dimension(:,:), allocatable :: coord
!
!     edges
!
      type(tedge), dimension(:), allocatable :: ed
!
!     elements
!
      type(tquad), dimension(:), allocatable :: el
!
      end module
!
!----------------------------------------------------------------------
!
      module mgauss
!
      use mgeometry
!
      implicit none
!
      integer, parameter :: nint1dmax=5
      integer, parameter :: nintmax  =25
      real(8), dimension(nintmax,nintmax,nsd) :: xi
      real(8), dimension(nintmax,nintmax)     :: wg
      real(8), dimension(nint1dmax,nint1dmax) :: xig
      real(8), dimension(nint1dmax,nint1dmax) :: wgl
      integer :: nint
!
      end module
!
!----------------------------------------------------------------------
!
      module mcoeficientes
!
      implicit none
!
      integer, parameter :: nsim = 200
! 
      ! real(8), dimension(5) :: alpha = (/6.305d-5,6.305d-5,6.305d-5,3.1525d-5,3.1525d-5/)
      real(8), dimension(5) :: alpha = (/3.1525d-5,3.1525d-5,3.1525d-5,1.9861d-6,1.9861d-6/)
      real(8) :: beta = 4.3834d-1 !3.1225d-1
      real(8) :: sig = 0.d0
      real(8) :: lamb = 0.d0
      real(8), parameter :: kp = 1.d0
      real(8) :: ikp = 0.d0
      real(8) :: mui = 0.d0
      real(8) :: mu = 0.d0
      real(8) :: delta = 0.06667d0 !0.2d0
      real(8) :: betav = 2.8865d-3 !4.2048d-3
      real(8) :: lambv = 0.07143d0!0.2d0
      real(8) :: kpv = 5.d0
      real(8) :: muv = 0.07143d0!0.2d0
      real(8) :: betas = 0.05d0/64.2853d0
!
      real(8), parameter :: a = 1.d0
      real(8) :: c1 = 1.d0
      real(8) :: c2 = 1.5d0
      real(8) :: c3 = 0.2333d0
      real(8) :: umin = 0.d0
      real(8) :: umax = 0.005d0 !0.001d0
!
      real(8) :: s0 = kp
      real(8) :: i0 = 0.01d0*kp
      real(8) :: r0 = 0.d0
      real(8) :: sv0 = 5.d0
      real(8) :: iv0 = 0.d0
!
      real(8), dimension(20) :: params = 0.d0
!
      real(8), parameter :: eps = 1.d-2!2

      real(8) :: i01 = 1.9806d-3!0.d0
      real(8) :: i02 = 0.1142d0
      real(8) :: i03 = 0.8401d-1
!
      real(8), dimension(20) :: paramlhs
      integer, dimension(12) :: indpar
!
!----------------------------------------------------------------------
!
      contains
!
!----------------------------------------------------------------------
! 
      subroutine iniparams
! 
      implicit none     
!
!     default parameters
!
      params(1)  = alpha(1)
      params(2)  = alpha(2)
      params(3)  = alpha(4)
      params(4)  =     beta 
      params(5)  =     lamb
      params(6)  =      mui
      params(7)  =       mu
      params(8)  =    delta 
      params(9)  =    betav
      params(10) =    lambv
      params(11) =      muv 
      params(12) =       c1 
      params(13) =       c2 
      params(14) =     umin
      params(15) =     umax
      params(16) =    betas 
      params(17) =      kpv
      params(18) =       c3 
!
      ! params(1)  = alpha(1)
      ! params(2)  = alpha(2)
      ! params(3)  = alpha(3)
      ! params(4)  = alpha(4)
      ! params(5)  = alpha(5)
      ! params(6)  =     beta 
      ! params(7)  =     lamb
      ! params(8)  =      mui
      ! params(9)  =       mu
      ! params(10) =    delta 
      ! params(11) =    betav
      ! params(12) =    lambv
      ! params(13) =      muv 
      ! params(14) =       c1 
      ! params(15) =       c2 
      ! params(16) =     umin
      ! params(17) =     umax 
! 
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine updtparam(param)
!
      implicit none
!
      real(8), dimension(*) :: param
!
      alpha(1) = param(1) !1.d-3
      alpha(2) = param(2) !5.d-4
      alpha(3) = param(1) !1.d-3
      alpha(4) = param(3) !1.d-4
      alpha(5) = param(3) !1.d-4
          beta = param(4) !1.28d-2
          lamb = param(5) !4.233d-5
           mui = param(6) !3.d-4
            mu = param(7) !1.849d-5
         delta = param(8) !1.15d-3
         betav = param(9) !4.d-2
         lambv = param(10) !8.881d-2
           muv = param(11) !3.73d-2
            c1 = param(12) !0.25d0
            c2 = param(13) !100.d0
          umin = param(14) !0.0d0
          umax = param(15) !0.001d0
         betas = param(16) !0.05d0
           kpv = param(17) !5.d0
            c3 = param(18) !1.d0   
!
      end subroutine
!
!--------------------------------------------------------------------------
!
      end module
!
!----------------------------------------------------------------------
!
      module marrays
!
      use mgeometry
!
      implicit none
!
!     states
!
      real(8), dimension(:,:), allocatable :: suh,inh,reh,suv,inv
      real(8), dimension(:,:), allocatable :: shold,ihold,rhold,svold,ivold
      real(8), dimension(:,:), allocatable :: intsuh,intinh,intreh,intsuv,intinv
      real(8), dimension(:,:), allocatable :: intctl,intj,intnnc,intvac
      real(8), dimension(:,:), allocatable :: suh0,inh0,reh0,suv0,inv0
      real(8), dimension(:,:), allocatable :: obsinh
!
!     adjoint
!
      real(8), dimension(:,:), allocatable :: l1,l2,l3,l4,l5
      real(8), dimension(:,:), allocatable :: l1old,l2old,l3old,l4old,l5old
      real(8), dimension(:,:), allocatable :: intl1,intl2,intl3,intl4,intl5
      real(8), dimension(:,:), allocatable :: l10,l20,l30,l40,l50
!
!     control
!
      real(8), dimension(:,:), allocatable :: u, uold
!
!     lapack: band structure
!
      integer :: kl,ku,lda,banda,ld
!
!     ilin jcol structure
!
      integer, dimension(:), allocatable :: ilin,jcol      
!
!     state matrices
!
      real(8), dimension(:,:), allocatable :: mesh,meih,merh,mesv,meiv
      real(8), dimension(:,:), allocatable :: mdsh,mdih,mdrh,mdsv,mdiv
      real(8), dimension(:,:), allocatable :: meshnl,meihnl,merhnl
      real(8), dimension(:,:), allocatable :: mesvnl,meivnl
      real(8), dimension(:,:), allocatable :: mdshnl,mdihnl,mdrhnl
      real(8), dimension(:,:), allocatable :: mdsvnl,mdivnl
      real(8), dimension(:,:), allocatable :: fsh,fih,frh,fsv,fiv
!
!     adjoint matrices
!
      real(8), dimension(:,:), allocatable :: mel1,mel2,mel3,mel4,mel5
      real(8), dimension(:,:), allocatable :: mdl1,mdl2,mdl3,mdl4,mdl5
      real(8), dimension(:,:), allocatable :: mel1nl,mel2nl,mel3nl
      real(8), dimension(:,:), allocatable :: mel4nl,mel5nl
      real(8), dimension(:,:), allocatable :: mdl1nl,mdl2nl,mdl3nl
      real(8), dimension(:,:), allocatable :: mdl4nl,mdl5nl
      real(8), dimension(:,:), allocatable :: fl1,fl2,fl3,fl4,fl5
!
      end module
!
!----------------------------------------------------------------------
!
      module merro
!
      use mcoeficientes, only: nsim
!
      implicit none
!
      integer :: indice
      real(8), dimension(nsim,5) :: errol2, errol2d, orders, ordersd
!
      end module
!
!----===============-------------------------------------------------------
      module msolver
!----===============-------------------------------------------------------
!
      implicit none
!
      real(8), parameter :: tol = 1.e-10!!12
      integer, parameter :: maxit = 100!50      
!
!--------------------------------------------------------------------------
      contains
!--------------------------------------------------------------------------
!
!----===============================================-----------------------
      subroutine diagonal(A,ilin,jcol,nelem,n,diag)
!----===============================================-----------------------
!     diagonal da matriz A
!--------------------------------------------------------------------------
!     declaracao das variaveis
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: i, n, nelem
      integer, dimension(*) :: ilin, jcol
      real(8), dimension(n) :: diag
      real(8), dimension(*) :: A
!--------------------------------------------------------------------------
!
      do i = 1, nelem
         if ( ilin(i).eq.jcol(i) ) then
            diag(ilin(i)) = A(i)
         end if
      end do
!
!--------------------------------------------------------------------------
      end subroutine
!--------------------------------------------------------------------------
!
!----============================================--------------------------
      subroutine prodax(A,ilin,jcol,nelem,n,x,y)
!----============================================--------------------------
!     produto interno y = A*x
!--------------------------------------------------------------------------
!     declaracao das variaveis
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: i, j, n, nelem
      integer, dimension(*) :: ilin, jcol
      real(8), dimension(n) :: x, y
      real(8), dimension(*) :: A
!--------------------------------------------------------------------------
!
      y = 0.0d0
!
      do i = 1, nelem
            y(ilin(i)) = y(ilin(i)) + A(i)*x(jcol(i))
      end do
!
!--------------------------------------------------------------------------
      end subroutine
!--------------------------------------------------------------------------
!
!----================================--------------------------------------
      subroutine prodscal(x,y,n,ytx)
!----================================--------------------------------------
!     produto escalar ytx = y^t*x
!--------------------------------------------------------------------------
!     declaracao das variaveis
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: i, n
      real(8) :: ytx
      real(8), dimension(*) :: x, y
!--------------------------------------------------------------------------
!
      ytx = 0.0d0
!
      do i = 1, n
         ytx = ytx + x(i)*y(i)
      end do
!
!--------------------------------------------------------------------------
      end subroutine
!--------------------------------------------------------------------------
!
!----===================================================-------------------
      subroutine prodxtax(A,ilin,jcol,nelem,n,x,y,xtAx)
!----===================================================-------------------
!     produto xtAx = x^t*A*x
!--------------------------------------------------------------------------
!     declaracao das variaveis
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: n, nelem
      integer, dimension(*) :: ilin, jcol
      real(8) :: xtAx
      real(8), dimension(*) :: x, y
      real(8), dimension(*) :: A
!--------------------------------------------------------------------------
!
      call prodax(A,ilin,jcol,nelem,n,x,y)
!
      call prodscal(x,y,n,xtAx)
!
!--------------------------------------------------------------------------
      end subroutine
!--------------------------------------------------------------------------
!
!----=============================================-------------------------
      subroutine respcg(A,ilin,jcol,nelem,x,b,n,r)
!----=============================================-------------------------
!     residuo r = b - A*x
!--------------------------------------------------------------------------
!     declaracao das variaveis
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: n, nelem
      integer, dimension(*) :: ilin, jcol
      real(8), dimension(n) :: x, r, b
      real(8), dimension(*) :: A
!--------------------------------------------------------------------------
!
      call prodax(A,ilin,jcol,nelem,n,x,r)
!
      r = b - r
!
!--------------------------------------------------------------------------
      end subroutine
!--------------------------------------------------------------------------
!
!----============================================================----------
      subroutine pcgdiag(A,ilin,jcol,nelem,b,x,diag,n,tol,maxit)
!----============================================================----------
!     gradientes conjugados com pre condicionamento diagonal
!--------------------------------------------------------------------------
!     declaracao das variaveis
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: n, nelem, maxit, it
      integer, dimension(*) :: ilin, jcol
      real(8) :: step, beta, gama, erro, alfa, errorel, tol
      real(8), dimension(n) :: x, r, b, diag, Ap, p, z
      real(8), dimension(*) :: A
!--------------------------------------------------------------------------
!     iteracao inicial
!
      it = 0
!
      call respcg(A,ilin,jcol,nelem,x,b,n,r)
      z = r/diag
      p = z
      call prodscal(r,z,n,alfa)
!
!     ---------------------------------
      errorel = dsqrt(alfa)
      if ( errorel.le.tol ) then
!         write(*,200) it, errorel
         return
      end if
!
!     inicio das iteracoes
!
      do
!
      it = it + 1
!
!     calculo do step = (r^t*r) / (p^t*A*p)
!
      call prodxtax(A,ilin,jcol,nelem,n,p,Ap,gama)
      step = alfa/gama
!
!     nova aproximacao
!
      x = x + step*p
!
!     calculo do residuo
!
      r = r - step*Ap
      z = r/diag
      call prodscal(p,p,n,erro)
!
!     calculo da nova direcao A-conjugada
!
      beta = alfa
      call prodscal(r,z,n,alfa)
      beta = alfa/beta
      p = z + beta*p
!
      erro = dabs(step)*dsqrt(erro)
      errorel = dsqrt(alfa)
!
!     teste de convergencia
!
      if( (erro.le.tol.and.errorel.le.tol).or.(it.gt.maxit) ) then
!         write(*,200) it, erro, errorel
!200      format(i3,3x,e10.4,3x,e10.4)
         return
      end if
!
      end do
!--------------------------------------------------------------------------
      end subroutine
!--------------------------------------------------------------------------
!
!----==============================================================--------
      subroutine pcgsqdiag(A,ilin,jcol,nelem,b,x,diag,n,tol,maxit)
!----==============================================================--------
!     residuo r = b - A*x
!--------------------------------------------------------------------------
!      declaracao das variaveis
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: n, nelem, maxit, it, i
      integer, dimension(*) :: ilin, jcol
      real(8) :: step,beta,gama,sigma,erro,alfa,errorel,tol,aux,xnew
      real(8), dimension(n) :: x, b, diag, Ap, p, q, r, r0, u
      real(8), dimension(*) :: A
!--------------------------------------------------------------------------
!     iteracao inicial
!
      it = 0
      gama = 1.0d0
!
      call respcg(A,ilin,jcol,nelem,x,b,n,r)
      q = 0.0d0
      p = 0.0d0
      r = r/diag
      r0 = r
!
!     inicio das iteracoes
!
      do
!
      it = it + 1
!
      call prodscal(r,r0,n,beta)
      alfa = beta
!
      if( dabs(gama).le.tol ) then
         return
      end if
!
      beta = beta/gama
!
      u = r + beta*q
      p = u + beta*( q + beta*p )
!
      call prodax(A,ilin,jcol,nelem,n,p,q)
!
      q = q/diag
!
      call prodscal(r0,q,n,sigma)
!
      gama = alfa
!
      if( dabs(gama).le.tol ) then
!         write(*,*) gama
         return
      end if
!
      alfa = alfa/sigma
!
!     nova aproximacao
!
      erro = 0.0d0
!
      do i = 1, n
         q(i) = u(i) - alfa*q(i)
         xnew = x(i) + alfa*( u(i) + q(i) )
         aux = dabs(xnew - x(i))
         if( erro.lt.aux ) then
         erro = aux
         end if
         x(i) = xnew
      end do
!
!     calculo do residuo
!
      call respcg(A,ilin,jcol,nelem,x,b,n,r)
!
      r = r/diag
!
      call prodscal(r,r,n,errorel)
      errorel = dsqrt(errorel)
!
!     teste de convergencia
!
      if( (erro.le.tol.and.errorel.le.tol).or.(it.gt.maxit) ) then
!         write(*,200) it, erro, errorel
200      format(i2,3x,e10.4,3x,e10.4)
         return
      end if
!
      end do
!--------------------------------------------------------------------------
      end subroutine
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
      end module
!----------------------------------------------------------------------
!
      module mpdes
!
      ! use mgeometry
      ! use marrays
      ! use mcoeficientes
!
      implicit none
!
      contains
!
!----------------------------------------------------------------------
!
      subroutine mosquitoes
!
      use mgeometry
      use marrays
      use mcoeficientes
      use merro
!
      implicit none
!
      integer :: i,j,ii,k
      real(8) :: cflag
!
!     allocating
!
!     non linear matrices
      allocate(meshnl(nmax,1))
      allocate(meihnl(nmax,1))
      allocate(mesvnl(nmax,1))
!
!     rhs vectors
!
      allocate(fsh(ndof,1))
      allocate(fih(ndof,1))
      allocate(frh(ndof,1))
      allocate(fsv(ndof,1))
      allocate(fiv(ndof,1))
!
      allocate(suh0(ndof,1))
      allocate(inh0(ndof,1))
      allocate(reh0(ndof,1))
      allocate(suv0(ndof,1))
      allocate(inv0(ndof,1))
!
!     initial conditions
!
      ! call initial
!
      ! call grafsurf(suh,'suh')
      ! call grafsurf(inh,'inh')
      ! call grafsurf(reh,'reh')
      ! call grafsurf(suv,'suv')
      ! call grafsurf(inv,'inv')
      ! call graftemp(suh,suht,'suh',intsuh)
      ! call graftemp(inh,inht,'inh',intinh)
      ! call graftemp(reh,reht,'reh',intreh)
      ! call graftemp(suv,suvt,'suv',intsuv)
      ! call graftemp(inv,invt,'inv',intinv)
!
!     temporal iterations
!
      do it=nt1,nt2!1,nt
!
      tt=(it-1)*dt
!
      ! do ii=1,4
      cflag = -1.d0
      ii = 0
      do while ((cflag<0.d0).and.(ii<20))
!
      ii = ii + 1
!
!     assemble non linear matrices
!
      call assemblenl
!
!     solving the global system
!
      call solve
!
      call pctest(cflag,ii)
!
      ! end do !ii
      end do !while
!
!........................................
!       do j=1,3
! !
! !       k=int((3)*7/dt+1)
! !       inh(el(156)%n(j),k) = 2.d0/4.3646d-3/10195921.63d0/0.2d0 !2
! !       inv(el(156)%n(j),k) = 5.d0*inh(el(156)%n(j),k)
! !       ! suh(el(156)%n(j),k) = suh(el(156)%n(j),k) - inh(el(156)%n(j),k)
! !       ! suv(el(156)%n(j),k) = suv(el(156)%n(j),k) - inv(el(156)%n(j),k)
! ! !
! !       ! k=int((6)*7/dt+1)
! !       ! inh(el(316)%n(j),int((6)*7/dt+1)) = 1.d0/4.3646d-3/10195921.63d0/0.2d0 !1
! !       ! inv(el(316)%n(j),k) = 5.d0*inh(el(316)%n(j),k)
! !       ! suh(el(316)%n(j),k) = suh(el(316)%n(j),k) - inh(el(316)%n(j),k)
! !       ! suv(el(316)%n(j),k) = suv(el(316)%n(j),k) - inv(el(316)%n(j),k)
! ! ! 
! !       k=int((6)*7/dt+1)
! !       inh(el(247)%n(j),int((6)*7/dt+1)) = 2.d0/4.3646d-3/10195921.63d0/0.2d0 !2
! !       inv(el(247)%n(j),k) = 5.d0*inh(el(247)%n(j),k)
! !       ! suh(el(247)%n(j),k) = suh(el(247)%n(j),k) - inh(el(247)%n(j),k)
! !       ! suv(el(247)%n(j),k) = suv(el(247)%n(j),k) - inv(el(247)%n(j),k)
! ! 
!       k=int((8)*7/dt+1)
!       inh(el(156)%n(j),k) = 25.d0/4.3646d-3/10195921.63d0/0.2d0 !2
!       inv(el(156)%n(j),k) = 5.d0*inh(el(156)%n(j),k)
!       ! suh(el(156)%n(j),k) = suh(el(156)%n(j),k) - inh(el(156)%n(j),k)
!       ! suv(el(156)%n(j),k) = suv(el(156)%n(j),k) - inv(el(156)%n(j),k)
! ! 
!       k=int((10)*7/dt+1)
!       inh(el(247)%n(j),k) = 87.d0/4.3646d-3/10195921.63d0/0.2d0 !2
!       inv(el(247)%n(j),k) = 5.d0*inh(el(247)%n(j),k)
!       ! suh(el(247)%n(j),k) = suh(el(247)%n(j),k) - inh(el(247)%n(j),k)
!       ! suv(el(247)%n(j),k) = suv(el(247)%n(j),k) - inv(el(247)%n(j),k)
! !
!       end do !j
!........................................
!
      end do !it
!
      ! call integral
!
      ! call grafsurf(suh,'suh')
      ! call grafsurf(inh,'inh')
      ! call grafsurf(reh,'reh')
      ! call grafsurf(suv,'suv')
      ! call grafsurf(inv,'inv')
      ! call grafsurf(l1,'l1')
      ! call grafsurf(l2,'l2')
      ! call grafsurf(l3,'l3')
      ! call grafsurf(l4,'l4')
      ! call grafsurf(l5,'l5')
!
!     calculating errors
!
      ! call norms
!
      ! open(unit=300,file='test.dat',status='unknown')
      ! write(300,*) suh
      ! close(unit=300)
!
!     deallocating
      deallocate(meshnl,meihnl,mesvnl)
      deallocate(fsh,fih,frh,fsv,fiv)
      deallocate(suh0,inh0,reh0,suv0,inv0)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine adjoint
!
      use mgeometry
      use marrays
      use mcoeficientes
      use merro
!
      implicit none
!
      integer :: i,j,ii
      real(8) :: cflag
!
      ndof = nnodes !2*nelem + 1
!
!     for lapack
!
      call calcband
!
!     allocating
!
!     system matrices
      allocate(mel1(nmax,1))
      allocate(mel2(nmax,1))
      allocate(mel3(nmax,1))
      allocate(mel4(nmax,1))
      allocate(mel5(nmax,1))
      allocate(mdl1(nmax,1))
      allocate(mdl2(nmax,1))
      allocate(mdl3(nmax,1))
      allocate(mdl4(nmax,1))
      allocate(mdl5(nmax,1))
!
!     rhs vectors
!
      allocate(fl1(ndof,1))
      allocate(fl2(ndof,1))
      allocate(fl3(ndof,1))
      allocate(fl4(ndof,1))
      allocate(fl5(ndof,1))
!
!     solutions
!
      allocate(intl1(nt+1,1))
      allocate(intl2(nt+1,1))
      allocate(intl3(nt+1,1))
      allocate(intl4(nt+1,1))
      allocate(intl5(nt+1,1))
!
      allocate(l10(ndof,1))
      allocate(l20(ndof,1))
      allocate(l30(ndof,1))
      allocate(l40(ndof,1))
      allocate(l50(ndof,1))
!
!     'initial' conditions
!
      call initialad
!
      ! call integralad(intl1,intl2,intl3,intl4,intl5)
!
      ! call graftemp(l1,l1t,'lb1',intl1)
      ! call graftemp(l2,l2t,'lb2',intl2)
      ! call graftemp(l3,l3t,'lb3',intl3)
      ! call graftemp(l4,l4t,'lb4',intl4)
      ! call graftemp(l5,l5t,'lb5',intl5)
!
!     temporal iterations
!
      do it=nt1,nt2!1,nt
!
!     backward time step
!
      itt = nt2+nt1+1-it!nt+2-it
!
      tt=(itt-1)*dt
!
!     assemble constant matrices
!
      call assemblead
!
      ! do ii=1,4
      cflag = -1.d0
      ii = 0
      do while ((cflag<0.d0).and.(ii<20))
!
      ii = ii + 1
!
!     assemble non linear matrices
!
      call assemblenlad
!
!     solving the global system
!
      call solvead
!
      call pctestad(cflag,ii)
!
      end do !while !ii
!
      ! call integralad(intl1,intl2,intl3,intl4,intl5)
!
      ! call graftemp(l1,l1t,'lb1',intl1)
      ! call graftemp(l2,l2t,'lb2',intl2)
      ! call graftemp(l3,l3t,'lb3',intl3)
      ! call graftemp(l4,l4t,'lb4',intl4)
      ! call graftemp(l5,l5t,'lb5',intl5)
! !
!       if (mod(it,nt/4).eq.0.or.it == nt/8) then
! !      if (mod(it,8).eq.0) then
!          call grafsurf(suh,it,'suh')
!          call grafsurf(inh,it,'inh')
!          call grafsurf(reh,it,'reh')
!          call grafsurf(suv,it,'suv')
!          call grafsurf(inv,it,'inv')
!       end if
!
      ! print *, expnum, iswp, 'adj', it
!
      end do !it
!
      ! call integralad
!
      ! call grafsurf(suh,'suh')
      ! call grafsurf(inh,'inh')
      ! call grafsurf(reh,'reh')
      ! call grafsurf(suv,'suv')
      ! call grafsurf(inv,'inv')
      ! call grafsurf(l1,'lb1')
      ! call grafsurf(l2,'lb2')
      ! call grafsurf(l3,'lb3')
      ! call grafsurf(l4,'lb4')
      ! call grafsurf(l5,'lb5')
!
!     calculating errors
!
      ! call normsad
!
!     deallocating
      deallocate(mel1,mel2,mel3,mel4,mel5)
      deallocate(mdl1,mdl2,mdl3,mdl4,mdl5)
      deallocate(fl1,fl2,fl3,fl4,fl5)
      deallocate(intl1,intl2,intl3,intl4,intl5)
      deallocate(l10,l20,l30,l40,l50)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine solve
!
!     resolve (mesh+meshnl)*shst = (mdsh)*suh + fsh
!                    (meih)*ihst = (mdsh)*inh + fih
!                    (merh)*rhst = (mdrh)*reh + frh
!             (mesv+mesvnl)*svst = (mdsv)*suv + fsv
!                    (meiv)*ivst = (mdiv)*inv + fiv
!
      use mgeometry
      use marrays
      use msolver
!
      implicit none
!
      integer :: i,info,ipiv(ndof),j
      real(8),dimension(ndof,1) :: ysh,yih,yrh,ysv,yiv
      real(8),dimension(ndof) :: diag
!
!     saving old solutions
!
      suh0(:,1) = suh(:,it+1)
      inh0(:,1) = inh(:,it+1)
      reh0(:,1) = reh(:,it+1)
      suv0(:,1) = suv(:,it+1)
      inv0(:,1) = inv(:,it+1)
!
!     ysh = (mdsh)*suh + fsh
!     yih = (mdih)*inh + fih
!     yrh = (mdrh)*reh + frh
!     ysv = (mdsv)*suv + fsv
!     yiv = (mdiv)*inv + fiv
!
      call prodax(mdsh,ilin,jcol,npos,ndof,suh(:,it),ysh)
      ysh = ysh + fsh
      call prodax(mdih,ilin,jcol,npos,ndof,inh(:,it),yih)
      yih = yih + fih
      call prodax(mdrh,ilin,jcol,npos,ndof,reh(:,it),yrh)
      yrh = yrh + frh
      call prodax(mdsv,ilin,jcol,npos,ndof,suv(:,it),ysv)
      ysv = ysv + fsv
      call prodax(mdiv,ilin,jcol,npos,ndof,inv(:,it),yiv)
      yiv = yiv + fiv
!
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,mdsh,lda,suh(:,it),1,0.d0,ysh,1)
      ! ysh = ysh + fsh
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,mdih,lda,inh(:,it),1,0.d0,yih,1)
      ! yih = yih + fih
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,mdrh,lda,reh(:,it),1,0.d0,yrh,1)
      ! yrh = yrh + frh
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,mdsv,lda,suv(:,it),1,0.d0,ysv,1)
      ! ysv = ysv + fsv
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,mdiv,lda,inv(:,it),1,0.d0,yiv,1)
      ! yiv = yiv + fiv
!
!     solving
!
      !  call dgbsv(ndof,kl,ku,1,(mesh+meshnl),lda,ipiv,ysh,ndof,info)
      !  call dgbsv(ndof,kl,ku,1,(meih),lda,ipiv,yih,ndof,info)
      !  call dgbsv(ndof,kl,ku,1,(merh),lda,ipiv,yrh,ndof,info)
      !  call dgbsv(ndof,kl,ku,1,(mesv+mesvnl),lda,ipiv,ysv,ndof,info)
      !  call dgbsv(ndof,kl,ku,1,(meiv),lda,ipiv,yiv,ndof,info)
      ! suh(:,it+1) = ysh(:,1)
      ! inh(:,it+1) = yih(:,1)
      ! reh(:,it+1) = yrh(:,1)
      ! suv(:,it+1) = ysv(:,1)
      ! inv(:,it+1) = yiv(:,1)
!
      call diagonal((mesh+meshnl),ilin,jcol,npos,ndof,diag)
      call pcgdiag((mesh+meshnl),ilin,jcol,npos,ysh,suh(:,it+1),diag,ndof,tol,maxit)

      call diagonal((meih+meihnl),ilin,jcol,npos,ndof,diag)
      call pcgdiag((meih+meihnl),ilin,jcol,npos,yih,inh(:,it+1),diag,ndof,tol,maxit)

      call diagonal(merh,ilin,jcol,npos,ndof,diag)
      call pcgdiag(merh,ilin,jcol,npos,yrh,reh(:,it+1),diag,ndof,tol,maxit)

      call diagonal((mesv+mesvnl),ilin,jcol,npos,ndof,diag)
      call pcgdiag((mesv+mesvnl),ilin,jcol,npos,ysv,suv(:,it+1),diag,ndof,tol,maxit)

      call diagonal(meiv,ilin,jcol,npos,ndof,diag)
      call pcgdiag(meiv,ilin,jcol,npos,yiv,inv(:,it+1),diag,ndof,tol,maxit)
!
      ! do i=1,ndof
      ! do j=1,nt+1
      ! if (suh(i,j) < 0.d0) then
      !    suh(i,j) = 0.d0
      ! end if
      ! if (inh(i,j) < 0.d0) then
      !    inh(i,j) = 0.d0
      ! end if
      ! if (reh(i,j) < 0.d0) then
      !    reh(i,j) = 0.d0
      ! end if
      ! if (suv(i,j) < 0.d0) then
      !    suv(i,j) = 0.d0
      ! end if
      ! if (inv(i,j) < 0.d0) then
      !    inv(i,j) = 0.d0
      ! end if
      ! end do
      ! end do
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine solvead
!
!     resolve (mel1)*l1st = (mdl1)*l1 + fl1
!             (mel2)*l2st = (mdl2)*l2 + fl2
!             (mel3)*l3st = (mdl3)*l3 + fl3
!             (mel4)*l4st = (mdl4)*l4 + fl4
!             (mel5)*l5st = (mdl5)*l5 + fl5
!
      use mgeometry
      use marrays
      use msolver
!
      implicit none
!
      integer :: i,info,ipiv(ndof)
      real(8),dimension(ndof,1) :: yl1,yl2,yl3,yl4,yl5
      real(8),dimension(ndof) :: diag
!
!     saving old solutions
!
      l10(:,1) = l1(:,itt-1)
      l20(:,1) = l2(:,itt-1)
      l30(:,1) = l3(:,itt-1)
      l40(:,1) = l4(:,itt-1)
      l50(:,1) = l5(:,itt-1)
!
!     yl1 = (mdl1)*l1 + fl1
!     yl2 = (mdl2)*l2 + fl2
!     yl3 = (mdl3)*l3 + fl3
!     yl4 = (mdl4)*l4 + fl4
!     yl5 = (mdl5)*l5 + fl5
!
      call prodax(mdl1,ilin,jcol,npos,ndof,l1(:,itt),yl1)
      yl1 = yl1 + fl1
      call prodax(mdl2,ilin,jcol,npos,ndof,l2(:,itt),yl2)
      yl2 = yl2 + fl2
      call prodax(mdl3,ilin,jcol,npos,ndof,l3(:,itt),yl3)
      yl3 = yl3 + fl3
      call prodax(mdl4,ilin,jcol,npos,ndof,l4(:,itt),yl4)
      yl4 = yl4 + fl4
      call prodax(mdl5,ilin,jcol,npos,ndof,l5(:,itt),yl5)
      yl5 = yl5 + fl5
!
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,(mdl1),lda,l1(:,itt),1,0.d0,yl1,1)
      ! yl1 = yl1 + fl1
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,(mdl2),lda,l2(:,itt),1,0.d0,yl2,1)
      ! yl2 = yl2 + fl2
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,(mdl3),lda,l3(:,itt),1,0.d0,yl3,1)
      ! yl3 = yl3 + fl3
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,(mdl4),lda,l4(:,itt),1,0.d0,yl4,1)
      ! yl4 = yl4 + fl4
      ! call dgbmv('n',ndof,ndof,kl,ku,1.d0,(mdl5),lda,l5(:,itt),1,0.d0,yl5,1)
      ! yl5 = yl5 + fl5
!
!     solving
!
      ! call dgbsv(ndof,kl,ku,1,(mel1),lda,ipiv,yl1,ndof,info)
      ! call dgbsv(ndof,kl,ku,1,(mel2),lda,ipiv,yl2,ndof,info)
      ! call dgbsv(ndof,kl,ku,1,(mel3),lda,ipiv,yl3,ndof,info)
      ! call dgbsv(ndof,kl,ku,1,(mel4),lda,ipiv,yl4,ndof,info)
      ! call dgbsv(ndof,kl,ku,1,(mel5),lda,ipiv,yl5,ndof,info)
      ! l1(:,itt-1) = yl1(:,1)
      ! l2(:,itt-1) = yl2(:,1)
      ! l3(:,itt-1) = yl3(:,1)
      ! l4(:,itt-1) = yl4(:,1)
      ! l5(:,itt-1) = yl5(:,1)
!
      call diagonal(mel1,ilin,jcol,npos,ndof,diag)
      call pcgdiag(mel1,ilin,jcol,npos,yl1,l1(:,itt-1),diag,ndof,tol,maxit)

      call diagonal(mel2,ilin,jcol,npos,ndof,diag)
      call pcgdiag(mel2,ilin,jcol,npos,yl2,l2(:,itt-1),diag,ndof,tol,maxit)

      call diagonal(mel3,ilin,jcol,npos,ndof,diag)
      call pcgdiag(mel3,ilin,jcol,npos,yl3,l3(:,itt-1),diag,ndof,tol,maxit)

      call diagonal(mel4,ilin,jcol,npos,ndof,diag)
      call pcgdiag(mel4,ilin,jcol,npos,yl4,l4(:,itt-1),diag,ndof,tol,maxit)

      call diagonal(mel5,ilin,jcol,npos,ndof,diag)
      call pcgdiag(mel5,ilin,jcol,npos,yl5,l5(:,itt-1),diag,ndof,tol,maxit)
!
      ! do i=1,ndof
      ! if (sst(i,1) < 1.e-8) then
      !    sst(i,1) = 0.d0
      ! end if
      ! if (ist(i,1) < 1.e-8) then
      !    ist(i,1) = 0.d0
      ! end if
      ! end do
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine pctest(cflag,ii)
!
!     test convergence for predictor corrector
!
      use mgeometry
      ! use mcoeficientes, only : eps
      use marrays
!
      implicit none
!
      integer :: i, j, ii
      real(8), dimension(5) :: ers, nrs, aux
      real(8) :: epsilon, cflag
!
      epsilon = 1d-6
!
!     max norms
!
      ers(1) = maxval(dabs(suh0(:,1)-suh(:,it+1)))
      ers(2) = maxval(dabs(inh0(:,1)-inh(:,it+1)))
      ers(3) = maxval(dabs(reh0(:,1)-reh(:,it+1)))
      ers(4) = maxval(dabs(suv0(:,1)-suv(:,it+1)))
      ers(5) = maxval(dabs(inv0(:,1)-inv(:,it+1)))
!
      nrs(1) = maxval(dabs(suh(:,it+1)))
      nrs(2) = maxval(dabs(inh(:,it+1)))
      nrs(3) = maxval(dabs(reh(:,it+1)))
      nrs(4) = maxval(dabs(suv(:,it+1)))
      nrs(5) = maxval(dabs(inv(:,it+1)))
!
!     2 norms
!
!      ers = 0.d0
!      nrs = 0.d0
!!
!      do i=1,ndof
!!
!      ers(1) = ers(1) + (dabs(suh(i,it)-suh(i,it+1)))**2
!      ers(2) = ers(2) + (dabs(inh(i,it)-inh(i,it+1)))**2
!      ers(3) = ers(3) + (dabs(reh(i,it)-reh(i,it+1)))**2
!      ers(4) = ers(4) + (dabs(suv(i,it)-suv(i,it+1)))**2
!      ers(5) = ers(5) + (dabs(inv(i,it)-inv(i,it+1)))**2
!!
!      nrs(1) = nrs(1) + (dabs(suh(i,it+1)))**2
!      nrs(2) = nrs(2) + (dabs(inh(i,it+1)))**2
!      nrs(3) = nrs(3) + (dabs(reh(i,it+1)))**2
!      nrs(4) = nrs(4) + (dabs(suv(i,it+1)))**2
!      nrs(5) = nrs(5) + (dabs(inv(i,it+1)))**2
!!
!      end do !i
! !
!       do i=1,5
!          ers(i) = dsqrt(ers(i))
!          nrs(i) = dsqrt(nrs(i))
!       end do
!
!     convergence criteria
!
      aux(1:5) = epsilon*nrs(1:5) - ers(1:5)
      ! aux(1:5) = epsilon - ers(1:5)
      cflag = minval(aux)
!
      ! print *, expnum, it, ii, 'ste', maxval(ers)
      ! write(*,'(i2,x,i2,x,i4,x,i2,x,e10.3,x,i2)') expnum,iswp,it,ii,maxval(ers),1
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine pctestad(cflag,ii)
!
!     test convergence for predictor corrector
!
      use mgeometry
      ! use mcoeficientes, only : eps
      use marrays
!
      implicit none
!
      integer :: i, j, ii
      real(8), dimension(5) :: ers, nrs, aux
      real(8) :: epsilon, cflag
!
      epsilon = 1d-6
!
!     max norms
!
      ers(1) = maxval(dabs(l10(:,1)-l1(:,itt-1)))
      ers(2) = maxval(dabs(l20(:,1)-l2(:,itt-1)))
      ers(3) = maxval(dabs(l30(:,1)-l3(:,itt-1)))
      ers(4) = maxval(dabs(l40(:,1)-l4(:,itt-1)))
      ers(5) = maxval(dabs(l50(:,1)-l5(:,itt-1)))
!
      nrs(1) = maxval(dabs(l1(:,itt-1)))
      nrs(2) = maxval(dabs(l2(:,itt-1)))
      nrs(3) = maxval(dabs(l3(:,itt-1)))
      nrs(4) = maxval(dabs(l4(:,itt-1)))
      nrs(5) = maxval(dabs(l5(:,itt-1)))
!
!     2 norms
!
!      ers = 0.d0
!      nrs = 0.d0
!!
!      do i=1,ndof
!!
!      ers(1) = ers(1) + (dabs(suh(i,it)-suh(i,it+1)))**2
!      ers(2) = ers(2) + (dabs(inh(i,it)-inh(i,it+1)))**2
!      ers(3) = ers(3) + (dabs(reh(i,it)-reh(i,it+1)))**2
!      ers(4) = ers(4) + (dabs(suv(i,it)-suv(i,it+1)))**2
!      ers(5) = ers(5) + (dabs(inv(i,it)-inv(i,it+1)))**2
!!
!      nrs(1) = nrs(1) + (dabs(suh(i,it+1)))**2
!      nrs(2) = nrs(2) + (dabs(inh(i,it+1)))**2
!      nrs(3) = nrs(3) + (dabs(reh(i,it+1)))**2
!      nrs(4) = nrs(4) + (dabs(suv(i,it+1)))**2
!      nrs(5) = nrs(5) + (dabs(inv(i,it+1)))**2
!!
!      end do !i
!
      ! do i=1,5
      !    ers(i) = dsqrt(ers(i))
      !    nrs(i) = dsqrt(nrs(i))
      ! end do
!
!     convergence criteria
!
     aux(1:5) = epsilon*nrs(1:5) - ers(1:5)
      ! aux(1:5) = epsilon - ers(1:5)
      cflag = minval(aux)
!
      ! print *, expnum, it, ii, 'adj', maxval(ers)
      ! write(*,'(i2,x,i2,x,i4,x,i2,x,e10.3,x,i2)') expnum,iswp,it,ii,maxval(ers),2
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine converg(cflag)
!
!     test convergence for state, adjoints and control
!
      use mgeometry
      use mcoeficientes, only : eps
      use marrays
!
      implicit none
!
      integer :: i, j
      real(8), dimension(5) :: ers, erl, nrs, nrl
      real(8), dimension(11) :: aux
      real(8) :: cflag, eru, nru
!
!     norms
!
      aux = 0.d0
      ers = 0.d0
      ! erl = 0.d0!
      eru = 0.d0
      nrs = 0.d0
      ! nrl = 0.d0!
      nru = 0.d0
!
      do i=1,ndof
      do j=1,(nt+1)
!
      ers(1) = ers(1) + dabs(suh(i,j)-shold(i,j))
      ers(2) = ers(2) + dabs(inh(i,j)-ihold(i,j))
      ers(3) = ers(3) + dabs(reh(i,j)-rhold(i,j))
      ers(4) = ers(4) + dabs(suv(i,j)-svold(i,j))
      ers(5) = ers(5) + dabs(inv(i,j)-ivold(i,j))
!
      ! erl(1) = erl(1) + dabs(l1(i,j)-l1old(i,j))!
      ! erl(2) = erl(2) + dabs(l2(i,j)-l2old(i,j))!
      ! erl(3) = erl(3) + dabs(l3(i,j)-l3old(i,j))!
      ! erl(4) = erl(4) + dabs(l4(i,j)-l4old(i,j))!
      ! erl(5) = erl(5) + dabs(l5(i,j)-l5old(i,j))!
!
      nrs(1) = nrs(1) + dabs(suh(i,j))
      nrs(2) = nrs(2) + dabs(inh(i,j))
      nrs(3) = nrs(3) + dabs(reh(i,j))
      nrs(4) = nrs(4) + dabs(suv(i,j))
      nrs(5) = nrs(5) + dabs(inv(i,j))
!
      ! nrl(1) = nrl(1) + dabs(l1(i,j)) !
      ! nrl(2) = nrl(2) + dabs(l2(i,j)) !
      ! nrl(3) = nrl(3) + dabs(l3(i,j)) !
      ! nrl(4) = nrl(4) + dabs(l4(i,j)) !
      ! nrl(5) = nrl(5) + dabs(l5(i,j)) !
!
      eru = eru + dabs(u(i,j)-uold(i,j))
      nru = nru + dabs(u(i,j))
!
      end do !j
      end do !i
!
!     convergence criteria
!
      aux(1:5) = eps*nrs(1:5) - ers(1:5)
      ! aux(6:10) = eps*nrl(1:5) - erl(1:5)!
      aux(11) = eps*nru - eru
!
      cflag = minval(aux)
   ! write(1,*) cflag
      write(*,'(i3,x,i2,x,x,e10.3)') expnum,iswp,cflag
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine updtctrl
!
!     update control
!
      use mgeometry
      use mcoeficientes
      use marrays
!
      implicit none
!
      integer :: i, j
!
!     norms
!
   do i=1,ndof
   do j=1,(nt+1)
!
      u(i,j) = (min(umax,max((l1(i,j)-l3(i,j)-c2)*suh(i,j)/(2.d0*c3)&
         ,umin))+uold(i,j))*0.5d0
      ! u(i,j) = min(umax,max((l1(i,j)-l3(i,j)-c2)*suh(i,j)/(2.d0*c3)&
      !  ,0.d0))!+0.1d0*uold(i,j)
      ! u(i,j) = (1.d0-(0.5d0**iswp))*min(umax,max((l1(i,j)-l3(i,j)-c2)*suh(i,j)/(2.d0*c3)&
      !       ,0.d0))+(0.5d0**iswp)*uold(i,j)
!
   end do !j
   end do !i
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine calcband
!
      use mgeometry
      use marrays
      implicit none
!
      select case(nen)
!
      case(3)
!
      banda=nx+2
!
      case(4)
!
      ! banda=nx+2
      banda = nnodes/2-1
!
      case(9)
!
      banda=2*(2*nx+1)+2
!
      case(8)
!
      banda=2*(2*nx+1)+2
!
      case(16)
!
      banda=3*(3*nx+1)+3
!
      case default
!
      write(*,*) 'element degree not defined'
      stop
!
      end select
!
      ku=banda
      kl=ku
!
      lda=ku+1+kl+kl
      ld =kl+ku+1  ! posicao da diagonal
!
      end subroutine
!
!----------------------------------------------------------------------
!
      function ibanda(i,j)
!
      use marrays
!
      implicit none
!
      integer             :: ibanda
      integer, intent(in) :: i,j
!
      ibanda = kl+ku+1+i-j
!
      end function
!
!----------------------------------------------------------------------
!
      function jbanda(i,j)
!
      use marrays
!
      implicit none
!
      integer             :: jbanda
      integer, intent(in) :: i,j
!
      jbanda = j
!
      end function
!
!----------------------------------------------------------------------
!
      subroutine genmesh
!
      use mgeometry
      implicit none
!
      integer :: i,j,l,k,l1,l2,l3,l4
      integer :: nnx,nny,nk,ni
!
      dt = tf/nt
!
!      nx=8
!      ny=8

      select case(nen)
!
      case(3)
!
!     linear - triangles
!
      hx=xl/nx
      hy=yl/ny
!
      nnodes= (nx+1)*(ny+1)
      nnx=nx+1
      nny=nx+1
      nelem = 2*nx*ny
      nedges= 3*nx*ny+nx+ny
!
!     allocating arrays
!
      allocate(coord(nnodes,nsd))
      allocate(el(nelem))
      allocate(ed(nedges))
!
!     coordinates
!
      do j=1,ny+1
      do i=1,nx+1
      l=i+(nx+1)*(j-1)
      coord( l,1) = (i-1)*hx
      coord( l,2) = (j-1)*hy
      end do ! i
      end do ! j
!
!      do i=1,nnodes
!      write(*,*) (coord(i,j),j=1,2)
!      end do
!      stop
!
!     conectivities: nodes
!
!      do j=1,ny
!      do i=1,nx
!      l=i+(nx)*(j-1)
!      el(l)%n(1)=l+j-1
!      el(l)%n(2)=l+j
!      el(l)%n(3)=l+j+nx+1
!      el(l)%n(4)=l+j+nx
!      end do ! i
!      end do ! j
      nk = 0
!
      do j = 1, ny
      ni = (j - 1)*nny
      do i = 1, nx
      k = nk + 1
!
      el(k)%n(1)= ni + i
      el(k)%n(2) = el(k)%n(1) + 1
      el(k)%n(3) = el(k)%n(2) + nx
!
      nk = k + 1
!
      el(nk)%n(1) = el(k)%n(3) + 1
      el(nk)%n(2) = el(k)%n(3)
      el(nk)%n(3) = el(k)%n(2)
      end do
      end do
!
!      do i=1,nelem
!      write(*,*) (el(i)%n(j),j=1,4)
!      end do
!      stop
!
      case(4)
!
!     bilinear
!
      hx=xl/nx
      hy=yl/ny
!
      nnodes= (nx+1)*(ny+1)
      nelem = nx*ny
      nedges= 2*nx*ny+nx+ny
!
!     allocating arrays
!
      allocate(coord(nnodes,nsd))
      allocate(el(nelem))
      allocate(ed(nedges))
!
!     coordinates
!
      do j=1,ny+1
      do i=1,nx+1
      l=i+(nx+1)*(j-1)
      coord( l,1) = (i-1)*hx
      coord( l,2) = (j-1)*hy
      end do ! i
      end do ! j
!
!      do i=1,nnodes
!      write(*,*) (coord(i,j),j=1,2)
!      end do
!      stop
!
!     conectivities: nodes
!
      do j=1,ny
      do i=1,nx
      l=i+(nx)*(j-1)
      el(l)%n(1)=l+j-1
      el(l)%n(2)=l+j
      el(l)%n(3)=l+j+nx+1
      el(l)%n(4)=l+j+nx
      ! s
!      do k = 1, nen
!         if (coord(el(l)%n(k),1)==0) then
!            el(l)%id(k) = 1
!         else if (coord(el(l)%n(k),2)==0) then
!            el(l)%id(k) = 1
!         else if (coord(el(l)%n(k),1)==xl) then
!            el(l)%id(k) = 1
!         else if (coord(el(l)%n(k),2)==yl) then
!            el(l)%id(k) = 1
!         else
!            el(l)%id(k) = 0
!         end if
!      end do ! k
      end do ! i
      end do ! j
!
!      do i=1,nelem
!      write(*,*) (el(i)%n(j),j=1,4)
!      end do
!      stop
!
      case(9)
!
!     biquadratic
      hx=xl/nx/2.d0
      hy=yl/ny/2.d0
!
      nnodes= (2*nx+1)*(2*ny+1)
      nelem = nx*ny
!      nedges= 2*nx*ny+nx+ny
!
!     allocating arrays
!
      allocate(coord(nnodes,nsd))
      allocate(el(nelem))
      allocate(ed(nedges))
!
!     coordinates
!
      do j=1,(2*ny+1)
      do i=1,(2*nx+1)
      l=i+(2*nx+1)*(j-1)
      coord( l,1) = (i-1)*hx
      coord( l,2) = (j-1)*hy
      end do ! i
      end do ! j
!
!      do i=1,nnodes
!      write(*,*) (coord(i,j),j=1,2)
!      end do
!      stop
!
!     conectivities: nodes
!
      do j=1,ny
      do i=1,nx
      l=i+(nx)*(j-1)
      el(l)%n(1)=2*(j-1)*(2*nx+1)+2*i-1
      el(l)%n(5)=el(l)%n(1)+1
      el(l)%n(2)=el(l)%n(1)+2

      el(l)%n(8)=el(l)%n(1)+(2*nx+1)
      el(l)%n(9)=el(l)%n(8)+1
      el(l)%n(6)=el(l)%n(8)+2


      el(l)%n(4)=el(l)%n(1)+2*(2*nx+1)
      el(l)%n(7)=el(l)%n(4)+1
      el(l)%n(3)=el(l)%n(4)+2
      ! contorno
!      do k = 1, nen
!         if ((dabs(coord(el(l)%n(k),1)))<1.e-8) then
!            el(l)%id(k) = 1
!         else if ((dabs(coord(el(l)%n(k),2)))<1.e-8) then
!            el(l)%id(k) = 1
!         else if ((dabs(coord(el(l)%n(k),1)-xl))<1.e-8) then
!            el(l)%id(k) = 2
!         else if ((dabs(coord(el(l)%n(k),2)-yl))<1.e-8) then
!            el(l)%id(k) = 2
!         else
!            el(l)%id(k) = 0
!         end if
!      end do ! k
      end do ! i
      end do ! j
!
!      do i=1,nelem
!      write(*,*) (el(i)%n(j),j=1,9)
!      end do
!      stop
!
      case(16)
!
!     bicubic
      hx=xl/nx/3.d0
      hy=yl/ny/3.d0
!
      nnodes= (3*nx+1)*(3*ny+1)
      nelem = nx*ny
!      nedges= 2*nx*ny+nx+ny
!
!     allocating arrays
!
      allocate(coord(nnodes,nsd))
      allocate(el(nelem))
      allocate(ed(nedges))
!
!     coordinates
!
      do j=1,(3*ny+1)
      do i=1,(3*nx+1)
      l=i+(3*nx+1)*(j-1)
      coord( l,1) = (i-1)*hx
      coord( l,2) = (j-1)*hy
      end do ! i
      end do ! j
!
!      do i=1,nnodes
!      write(*,*) (coord(i,j),j=1,2)
!      end do
!      stop
!
!     conectivities: nodes
!
      do j=1,ny
      do i=1,nx
      l=i+(nx)*(j-1)
      el(l)%n(1)=3*(j-1)*(3*nx+1)+3*i-2
      el(l)%n(5)=el(l)%n(1)+1
      el(l)%n(6)=el(l)%n(1)+2
      el(l)%n(2)=el(l)%n(1)+3

      el(l)%n(12)=el(l)%n(1)+(3*nx+1)
      el(l)%n(13)=el(l)%n(12)+1
      el(l)%n(14)=el(l)%n(12)+2
      el(l)%n(7)=el(l)%n(12)+3

      el(l)%n(11)=el(l)%n(1)+2*(3*nx+1)
      el(l)%n(16)=el(l)%n(11)+1
      el(l)%n(15)=el(l)%n(11)+2
      el(l)%n(8)=el(l)%n(11)+3

      el(l)%n(4)=el(l)%n(1)+3*(3*nx+1)
      el(l)%n(10)=el(l)%n(4)+1
      el(l)%n(9)=el(l)%n(4)+2
      el(l)%n(3)=el(l)%n(4)+3
      ! contorno
!      do k = 1, nen
!         if ((dabs(coord(el(l)%n(k),1)))<1.e-8) then
!            el(l)%id(k) = 1
!         else if ((dabs(coord(el(l)%n(k),2)))<1.e-8) then
!            el(l)%id(k) = 1
!         else if ((dabs(coord(el(l)%n(k),1)-xl))<1.e-8) then
!            el(l)%id(k) = 2
!         else if ((dabs(coord(el(l)%n(k),2)-yl))<1.e-8) then
!            el(l)%id(k) = 2
!         else
!            el(l)%id(k) = 0
!         end if
!      end do ! k
      end do ! i
      end do ! j
!
!      do i=1,nelem
!      write(*,*) (el(i)%n(j),j=1,16)
!      end do
!      stop
!
      case(8)
!
!     serendipity 8
!
      hx=xl/nx/2.d0
      hy=yl/ny/2.d0
!
      nnodes= (2*nx+1)*(2*ny+1) - nx*ny
      nelem = nx*ny
!      nedges= 2*nx*ny+nx+ny
!
!     allocating arrays
!
      allocate(coord(nnodes,nsd))
      allocate(el(nelem))
      allocate(ed(nedges))
!
!     coordinates
!
      l = 1
      do j=1,(ny)
      do i=1,(2*nx+1)

      coord( l,1) = (i-1)*hx
      coord( l,2) = (j-1)*2*hy
      l = l + 1
      end do ! i

      do i = 1,(nx+1)
      coord( l,1) = (i-1)*2*hx
      coord( l,2) = (2*j-1)*hy
      l = l + 1
      end do ! i

      end do ! j

      do i=1,(2*nx+1)

      coord( l,1) = (i-1)*hx
      coord( l,2) = 2*ny*hy
      l = l + 1
      end do ! i
!
!      do i=1,nnodes
!      write(*,*) i, (coord(i,j),j=1,2)
!      end do
!      stop
!
!     conectivities: nodes
!
      do j=1,ny
      do i=1,nx
      l=i+(nx)*(j-1)
      el(l)%n(1)=(j-1)*((nx+1)+2*nx+1)+2*i-1
      el(l)%n(5)=el(l)%n(1)+1
      el(l)%n(2)=el(l)%n(1)+2

!      el(l)%n(12)=el(l)%n(1)+(3*nx+1-(i-1)*2)
!      el(l)%n(7)=el(l)%n(12)+1

      el(l)%n(8)=el(l)%n(1)+(2*nx+1-(i-1))
      el(l)%n(6)=el(l)%n(8)+1

      el(l)%n(4)=el(l)%n(1)+(2*nx+1+(nx+1))
      el(l)%n(7)=el(l)%n(4)+1
      el(l)%n(3)=el(l)%n(4)+2

      ! contorno
!      do k = 1, nen
!         if ((dabs(coord(el(l)%n(k),1)))<1.e-8) then
!            el(l)%id(k) = 1
!         else if ((dabs(coord(el(l)%n(k),2)))<1.e-8) then
!            el(l)%id(k) = 1
!         else if ((dabs(coord(el(l)%n(k),1)-xl))<1.e-8) then
!            el(l)%id(k) = 2
!         else if ((dabs(coord(el(l)%n(k),2)-yl))<1.e-8) then
!            el(l)%id(k) = 2
!         else
!            el(l)%id(k) = 0
!         end if
!      end do ! k
      end do ! i
      end do ! j

!      do i=1,nelem
!      write(*,'(10i3x)') (el(i)%n(j),j=1,8)
!      end do
!      stop
!
      case default
!
      write(*,*) 'element degree not defined'
      stop
!
      end select
!
!     conectivities: element-edges
!
!     vertical edges
!
      k=0
      do j=1,ny
      do i=1,nx
      l=i+(nx)*(j-1)
      k=k+1
      el(l)%e(1)=k
      k=k+1
      el(l)%e(3)=k
      k=k-1
      end do ! i
      k=k+1
      end do ! j
!
!     horizontal edges
!
      do i=1,nx
      do j=1,ny
      l=i+(nx)*(j-1)
      k=k+1
      el(l)%e(2)=k
      k=k+1
      el(l)%e(4)=k
      k=k-1
      end do ! i
      k=k+1
      end do ! j
!
!     neighbors
!
      do i=1,nx
      do j=1,ny
!
      l =i+(nx)*(j-1)
      l1=l-1
      l2=l-nx
      l3=l+1
      l4=l+nx
!
      if(i.eq.1 ) l1=0
      if(j.eq.1 ) l2=0
      if(i.eq.nx) l3=0
      if(j.eq.ny) l4=0
!
      el(l)%ngb(1)=l1
      el(l)%ngb(2)=l2
      el(l)%ngb(3)=l3
      el(l)%ngb(4)=l4
!
      end do ! i
      end do ! j
!
!     indice para graficos temporais
!
!     posicao dos pontos na malha
!
!      ----------------
!     |  7    8    9   |
!     |  4    5    6   |
!     |  1    2    3   |
!      ----------------
!
!      ind(2) = (nx+2)/2 + ny/4*(nx+1)
!      ind(1) = ind(2) - nx/4
!      ind(3) = ind(2) + nx/4
!
!      ind(4) = ind(1) + ny/4*(nx+1)
!      ind(5) = ind(2) + ny/4*(nx+1)
!      ind(6) = ind(3) + ny/4*(nx+1)
!
!      ind(7) = ind(4) + ny/4*(nx+1)
!      ind(8) = ind(5) + ny/4*(nx+1)
!      ind(9) = ind(6) + ny/4*(nx+1)
!!      ----------------
!!     |  1    2    3   |
!!     |  4    5    6   |
!!     |  7    8    9   |
!!      ----------------
!!
      ind(8) = (nx+2)/2 + ny/4*(nx+1)
      ind(7) = ind(8) - nx/4
      ind(9) = ind(8) + nx/4

      ind(4) = ind(7) + ny/4*(nx+1)
      ind(5) = ind(8) + ny/4*(nx+1)
      ind(6) = ind(9) + ny/4*(nx+1)

      ind(1) = ind(4) + ny/4*(nx+1)
      ind(2) = ind(5) + ny/4*(nx+1)
      ind(3) = ind(6) + ny/4*(nx+1)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine readmesh
!
      use mgeometry
      implicit none
!
      integer :: i,j,l,k,l1,l2,l3,l4
      integer :: nnx,nny,nk,ni
      real(8) :: dist, dmax = 0.d0
!
      dt = tf/nt
!
!      nen = 4
!
!     bilinear
!
      open(unit=101,file='meshdata.dat',status='unknown')
      open(unit=102,file='meshcoord.dat',status='unknown')
      open(unit=103,file='meshconect.dat',status='unknown')
      open(unit=104,file='meshinindex.dat',status='unknown')
!
      read(101,*) nnodes
      read(101,*) nelem
      read(101,*) nedges
!
!     initial condition indexes
!
      ! do i=1,9
      ! read(104,*) ind(i)
      ! end do !i
! 
      ! write(*,*) nnodes, nelem, nedges
      ! stop
!
!     allocating arrays
!
      allocate(coord(nnodes,nsd))
      allocate(el(nelem))
      allocate(ed(nedges))
!
!     coordinates
!
      do i=1,nnodes
      read(102,*) coord(i,1), coord(i,2)
      end do !i
!
     ! do i=1,nnodes
     ! write(*,*) (coord(i,j),j=1,2)
     ! end do
     ! stop
!
!     conectivities: nodes
!
      do i=1,nelem
      ! read(103,*) el(i)%n(1),el(i)%n(2),el(i)%n(3),el(i)%n(4)
      read(103,*) el(i)%n(1),el(i)%n(2),el(i)%n(3)
      end do ! i
!
     ! do i=1,nelem
     ! write(*,*) (el(i)%n(j),j=1,4)
     ! end do
     ! stop
!
!     max element sizes
!
      do k=1,nelem
      do i=1,2!3
!
      j = el(k)%n(i)
      l = el(k)%n(i+1)
!
      dist = dsqrt((coord(j,1)-coord(l,1))**2+(coord(j,2)-coord(l,2))**2)
!
      if (dist > dmax) then
         dmax = dist
      end if
!
      end do !i
!
      j = el(k)%n(1)
      l = el(k)%n(3)!el(k)%n(4)
!
      dist = dsqrt((coord(j,1)-coord(l,1))**2+(coord(j,2)-coord(l,2))**2)
!
      if (dist > dmax) then
         dmax = dist
      end if
!
      end do !k
!
      hx = dmax
!
      ! write(*,*) dmax
      ! stop
!
!     conectivities: element-edges
!
!     neighbors
!
!     indexes
!
      ! ind(8) = (nx+2)/2 + ny/4*(nx+1)
      ! ind(7) = ind(8) - nx/4
      ! ind(9) = ind(8) + nx/4

      ! ind(4) = ind(7) + ny/4*(nx+1)
      ! ind(5) = ind(8) + ny/4*(nx+1)
      ! ind(6) = ind(9) + ny/4*(nx+1)

      ! ind(1) = ind(4) + ny/4*(nx+1)
      ! ind(2) = ind(5) + ny/4*(nx+1)
      ! ind(3) = ind(6) + ny/4*(nx+1)

!
      close(unit=101)
      close(unit=102)
      close(unit=103)
      close(unit=104)
!
      end subroutine
!
!----------------------------------------------------------------------
!
!       subroutine genmeshnonh
! !
!       use mgeometry
!       implicit none
! !
!       integer :: i,j,l,k,l1,l2,l3,l4
!       real(8) :: xl,yl
! !
!       xl=1.d0
!       yl=1.d0
! !
! !      nx=8
! !      ny=8

!       select case(nen)
! !
!       case(4)
! !
! !     bilinear
! !
!       hx=xl/nx
!       hy=yl/ny
! !
!       nnodes= (nx+1)*(ny+1)
!       nelem = nx*ny
!       nedges= 2*nx*ny+nx+ny
! !
! !     allocating arrays
! !
!       allocate(coord(nnodes,nsd))
!       allocate(el(nelem))
!       allocate(ed(nedges))
! !
! !     coordinates
! !
!       do j=1,ny+1
!       do i=1,nx+1
!       l=i+(nx+1)*(j-1)
!       coord( l,1) = (i-1)*hx
!       coord( l,2) = (j-1)*hy
!       end do ! i
!       end do ! j
! !
! !     coordinates-nonhomogneous
! !
!       do j=2,ny,2
!       do i=2,nx,2
!       l=i+(nx+1)*(j-1)
!       coord( l,2) = coord( l,2)+2.d0/3.d0*hy!0.25d0*hy
!       end do ! i
!       end do ! j
! !
!       do j=2,ny,2
!       do i=1,nx+1,2
!       l=i+(nx+1)*(j-1)
!       coord( l,2) = coord( l,2)-2.d0/3.d0*hy!0.25d0*hy
!       end do ! i
!       end do ! j
! !
! !      do i=1,nnodes
! !      write(*,*) (coord(i,j),j=1,2)
! !      end do
! !      stop
! !
! !     conectivities: nodes
! !
!       do j=1,ny
!       do i=1,nx
!       l=i+(nx)*(j-1)
!       el(l)%n(1)=l+j-1
!       el(l)%n(2)=l+j
!       el(l)%n(3)=l+j+nx+1
!       el(l)%n(4)=l+j+nx
!       ! contorno
!       do k = 1, nen
!          if (coord(el(l)%n(k),1)==0) then
!             el(l)%id(k) = 2
!          else if (coord(el(l)%n(k),2)==0) then
!             el(l)%id(k) = 2
!          else if (coord(el(l)%n(k),1)==xl) then
!             el(l)%id(k) = 2
!          else if (coord(el(l)%n(k),2)==yl) then
!             el(l)%id(k) = 2
!          else
!             el(l)%id(k) = 0
!          end if
!       end do ! k
!       end do ! i
!       end do ! j
! !
! !      do i=1,nelem
! !      write(*,*) (el(i)%n(j),j=1,4)
! !      end do
! !      stop
! !
!       case(9)
! !
! !     biquadratic
!       hx=xl/nx/2.d0
!       hy=yl/ny/2.d0
! !
!       nnodes= (2*nx+1)*(2*ny+1)
!       nelem = nx*ny
! !      nedges= 2*nx*ny+nx+ny
! !
! !     allocating arrays
! !
!       allocate(coord(nnodes,nsd))
!       allocate(el(nelem))
!       allocate(ed(nedges))
! !
! !     coordinates
! !
!       do j=1,(2*ny+1)
!       do i=1,(2*nx+1)
!       l=i+(2*nx+1)*(j-1)
!       coord( l,1) = (i-1)*hx
!       coord( l,2) = (j-1)*hy
!       end do ! i
!       end do ! j
! !
! !     coordinates-nonhomogneous
! !
!       do j=2,2*ny,4
!       do i=2,2*nx,4
!       l=i+1+(2*nx+1)*j
!       coord( l,2) = coord( l,2)+2.d0/3.d0*hy!0.5d0*hy

!       coord( l-(2*nx+1),2) = coord( l-(2*nx+1),2)+1.d0/3.d0*hy!0.5d0*hy 8
!       coord( l+(2*nx+1),2) = coord( l+(2*nx+1),2)+1.d0/3.d0*hy!0.5d0*hy 18

!       end do ! i
!       end do ! j
! !
!       do j=2,2*ny,4
!       do i=1,2*nx+1,4
!       l=i+(2*nx+1)*j
!       coord( l,2) = coord( l,2)-2.d0/3.d0*hy!0.5d0*hy

!       coord( l-(2*nx+1),2) = coord( l-(2*nx+1),2)-1.d0/3.d0*hy!0.5d0*hy 6 10
!       coord( l+(2*nx+1),2) = coord( l+(2*nx+1),2)-1.d0/3.d0*hy!0.5d0*hy 16 20
!       end do ! i
!       end do ! j
! !
! !      do i=1,nnodes
! !      write(*,*) (coord(i,j),j=1,2)
! !      end do
! !      stop
! !
! !     conectivities: nodes
! !
!       do j=1,ny
!       do i=1,nx
!       l=i+(nx)*(j-1)
!       el(l)%n(1)=2*(j-1)*(2*nx+1)+2*i-1
!       el(l)%n(5)=el(l)%n(1)+1
!       el(l)%n(2)=el(l)%n(1)+2

!       el(l)%n(8)=el(l)%n(1)+(2*nx+1)
!       el(l)%n(9)=el(l)%n(8)+1
!       el(l)%n(6)=el(l)%n(8)+2


!       el(l)%n(4)=el(l)%n(1)+2*(2*nx+1)
!       el(l)%n(7)=el(l)%n(4)+1
!       el(l)%n(3)=el(l)%n(4)+2
!       ! contorno
!       do k = 1, nen
!          if ((dabs(coord(el(l)%n(k),1)))<1.e-8) then
!             el(l)%id(k) = 2
!          else if ((dabs(coord(el(l)%n(k),2)))<1.e-8) then
!             el(l)%id(k) = 2
!          else if ((dabs(coord(el(l)%n(k),1)-xl))<1.e-8) then
!             el(l)%id(k) = 2
!          else if ((dabs(coord(el(l)%n(k),2)-yl))<1.e-8) then
!             el(l)%id(k) = 2
!          else
!             el(l)%id(k) = 0
!          end if
!       end do ! k
!       end do ! i
!       end do ! j
! !
! !      do i=1,nelem
! !      write(*,'(10i3x)') (el(i)%n(j),j=1,9)
! !      end do
! !      stop
! !
!       case(8)
! !
! !     serendipity 8
! !
!       hx=xl/nx/2.d0
!       hy=yl/ny/2.d0
! !
!       nnodes= (2*nx+1)*(2*ny+1) - nx*ny
!       nelem = nx*ny
! !      nedges= 2*nx*ny+nx+ny
! !
! !     allocating arrays
! !
!       allocate(coord(nnodes,nsd))
!       allocate(el(nelem))
!       allocate(ed(nedges))
! !
! !     coordinates
! !
!       l = 1
!       do j=1,(ny)
!       do i=1,(2*nx+1)

!       coord( l,1) = (i-1)*hx
!       coord( l,2) = (j-1)*2*hy
!       l = l + 1
!       end do ! i

!       do i = 1,(nx+1)
!       coord( l,1) = (i-1)*2*hx
!       coord( l,2) = (2*j-1)*hy
!       l = l + 1
!       end do ! i

!       end do ! j

!       do i=1,(2*nx+1)
!       coord( l,1) = (i-1)*hx
!       coord( l,2) = 2*ny*hy
!       l = l + 1
!       end do ! i
! !
! !     coordinates-nonhomogneous
! !
!       do j=1,ny/2
!       do i=4,2*nx,4
!       l=i-1+(3*nx+2)*(2*j-1)
!       coord( l,2) = coord( l,2)+2.d0/3.d0*hy!0.5d0*hy

!       coord( l+(nx-(i/2-1))*2+i/2,2) = coord( l+(nx-(i/2-1))*2+i/2,2)+1.d0/3.d0*hy!0.5d0*hy
!       coord( l-((i/2-1)*2+nx-(i/2-2)),2) = coord( l-((i/2-1)*2+nx-(i/2-2)),2)+1.d0/3.d0*hy!0.5d0*hy

!       end do ! i
!       end do ! j
! !
!       do j=1,ny/2
!       do i=1,2*nx+1,4
!       l=i+(3*nx+2)*(2*j-1)
!       coord( l,2) = coord( l,2)-2.d0/3.d0*hy!0.5d0*hy

!       coord( l+2*nx-(i-1)+(i-1)/2+1,2) = coord( l+2*nx-(i-1)+(i-1)/2+1,2)-1.d0/3.d0*hy!0.5d0*hy
!       coord( l-((i-1)+nx+1-(i-1)/2),2) = coord( l-((i-1)+nx+1-(i-1)/2),2)-1.d0/3.d0*hy!0.5d0*hy
!       end do ! i
!       end do ! j
! !
! !      do i=1,nnodes
! !      write(*,*) i, (coord(i,j),j=1,2)
! !      end do
! !      stop
! !
! !     conectivities: nodes
! !
!       do j=1,ny
!       do i=1,nx
!       l=i+(nx)*(j-1)
!       el(l)%n(1)=(j-1)*((nx+1)+2*nx+1)+2*i-1
!       el(l)%n(5)=el(l)%n(1)+1
!       el(l)%n(2)=el(l)%n(1)+2

!       el(l)%n(8)=el(l)%n(1)+(2*nx+1-(i-1))
!       el(l)%n(6)=el(l)%n(8)+1

!       el(l)%n(4)=el(l)%n(1)+(2*nx+1+(nx+1))
!       el(l)%n(7)=el(l)%n(4)+1
!       el(l)%n(3)=el(l)%n(4)+2

!       ! contorno
!       do k = 1, nen
!          if ((dabs(coord(el(l)%n(k),1)))<1.e-8) then
!             el(l)%id(k) = 2
!          else if ((dabs(coord(el(l)%n(k),2)))<1.e-8) then
!             el(l)%id(k) = 2
!          else if ((dabs(coord(el(l)%n(k),1)-xl))<1.e-8) then
!             el(l)%id(k) = 2
!          else if ((dabs(coord(el(l)%n(k),2)-yl))<1.e-8) then
!             el(l)%id(k) = 2
!          else
!             el(l)%id(k) = 0
!          end if
!       end do ! k
!       end do ! i
!       end do ! j

! !      do i=1,nelem
! !      write(*,'(10i3x)') (el(i)%n(j),j=1,8)
! !      end do
! !      stop
! !
!       case default
! !
!       write(*,*) 'element degree not defined'
!       stop
! !
!       end select
! !
! !     conectivities: element-edges
! !
! !     vertical edges
! !
!       k=0
!       do j=1,ny
!       do i=1,nx
!       l=i+(nx)*(j-1)
!       k=k+1
!       el(l)%e(1)=k
!       k=k+1
!       el(l)%e(3)=k
!       k=k-1
!       end do ! i
!       k=k+1
!       end do ! j
! !
! !     horizontal edges
! !
!       do i=1,nx
!       do j=1,ny
!       l=i+(nx)*(j-1)
!       k=k+1
!       el(l)%e(2)=k
!       k=k+1
!       el(l)%e(4)=k
!       k=k-1
!       end do ! i
!       k=k+1
!       end do ! j
! !
! !     neighbors
! !
!       do i=1,nx
!       do j=1,ny
! !
!       l =i+(nx)*(j-1)
!       l1=l-1
!       l2=l-nx
!       l3=l+1
!       l4=l+nx
! !
!       if(i.eq.1 ) l1=0
!       if(j.eq.1 ) l2=0
!       if(i.eq.nx) l3=0
!       if(j.eq.ny) l4=0
! !
!       el(l)%ngb(1)=l1
!       el(l)%ngb(2)=l2
!       el(l)%ngb(3)=l3
!       el(l)%ngb(4)=l4
! !
!       end do ! i
!       end do ! j

! !      do i=1,nelem
! !      write(*,'(10i3x)') (el(i)%ngb(j),j=1,4)
! !      end do
! !      stop
! !
!       end subroutine
!
!----------------------------------------------------------------------
!
      subroutine gausstable
!
!     store data for numerical integration
!
      use mgeometry
      use mgauss
!
      implicit none
!
      integer :: i,j,k,l,nint1d,nint2d
      integer, parameter :: ngt=5
!      real(8), dimension(ngt,ngt) :: xig,wgl
!
!     unidimensional
!
!     one point
!
      xig(1,1) = 0.d0
!
      wgl(1,1) = 2.d0
!
!     two points
!
      xig(1,2) = -dsqrt(3.d0)/3.d0
      xig(2,2) =  dsqrt(3.d0)/3.d0
!
      wgl(1,2) =  1.d0
      wgl(2,2) =  1.d0
!
!     three points
!
      xig(1,3) =  0.d0
      xig(2,3) = -dsqrt(3.d0/5.d0)
      xig(3,3) =  dsqrt(3.d0/5.d0)
!
      wgl(1,3) =  8.d0/9.d0
      wgl(2,3) =  5.d0/9.d0
      wgl(3,3) =  5.d0/9.d0
!
!     four points
!
      xig(1,4) = -dsqrt((3.d0-2.d0*dsqrt(6.d0/5.d0))/7.d0)
      xig(2,4) =  dsqrt((3.d0-2.d0*dsqrt(6.d0/5.d0))/7.d0)
      xig(3,4) = -dsqrt((3.d0+2.d0*dsqrt(6.d0/5.d0))/7.d0)
      xig(4,4) =  dsqrt((3.d0+2.d0*dsqrt(6.d0/5.d0))/7.d0)
!
      wgl(1,4) =  (18.d0+dsqrt(30.d0))/36.d0
      wgl(2,4) =  (18.d0+dsqrt(30.d0))/36.d0
      wgl(3,4) =  (18.d0-dsqrt(30.d0))/36.d0
      wgl(4,4) =  (18.d0-dsqrt(30.d0))/36.d0
!
!     five points
!
      xig(1,5) =  0.d0
      xig(2,5) = -dsqrt(5.d0-2.d0*dsqrt(10.d0/7.d0))/3.d0
      xig(3,5) =  dsqrt(5.d0-2.d0*dsqrt(10.d0/7.d0))/3.d0
      xig(4,5) = -dsqrt(5.d0+2.d0*dsqrt(10.d0/7.d0))/3.d0
      xig(5,5) =  dsqrt(5.d0+2.d0*dsqrt(10.d0/7.d0))/3.d0
!
      wgl(1,5) =  128.d0/225.d0
      wgl(2,5) =  (322.d0+13.d0*dsqrt(70.d0))/900.d0
      wgl(3,5) =  (322.d0+13.d0*dsqrt(70.d0))/900.d0
      wgl(4,5) =  (322.d0-13.d0*dsqrt(70.d0))/900.d0
      wgl(5,5) =  (322.d0-13.d0*dsqrt(70.d0))/900.d0
!
!     bidimensional
!
      do nint1d=1,ngt
!
      nint2d=nint1d**nsd
!
      l=0
!
      do j=1,nint1d
!
      do i=1,nint1d
!
      l=l+1
!
      xi(l,nint2d,1) = xig(i,nint1d)
      xi(l,nint2d,2) = xig(j,nint1d)
!
      wg(l,nint2d)=wgl(i,nint1d)*wgl(j,nint1d)
!
      end do ! i
!
      end do ! j
!
      end do ! nint1d
!
!
!     for triangular meshes
!
      xi(1,3,1) = 1.d0/6.d0
      xi(1,3,2) = 1.d0/6.d0
      xi(2,3,1) = 2.d0/3.d0
      xi(2,3,2) = 1.d0/6.d0
      xi(3,3,1) = 1.d0/6.d0
      xi(3,3,2) = 2.d0/3.d0    
!
      wg(1,3) = 1.d0/6.d0
      wg(2,3) = wg(1,3)
      wg(3,3) = wg(1,3)
!
!
!      do i=1,4
!      write(*,*) (xi(i,4,j),j=1,nsd)
!      end do ! i

!     numero de pontos de integracao
!      de acordo com o grau dos polinomios

      select case(nen)
!
      case(3)
!
      nint=3
!
      case(4)
!
      nint=4
!
      case(9)
!
      nint=16
!
      case(8)
!
      nint=16
!
      case(16)
!
      nint=25
!
      case default
!
      write(*,*) 'element degree not defined'
      stop
!
      end select
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine assemble
!
      use mgeometry
      use mcoeficientes
      use mgauss
      use marrays
!
      implicit none
!
      integer :: nel,i,j,k,l,nk,n1,n2,ii
      real(8), dimension(nsd,nsd) :: jf,invjf,tinvjf

      real(8), dimension(nsd)     :: xil,xxl,xx
      real(8), dimension(nen,nsd) :: xe
      real(8), dimension(nen,nen) :: meshe,meihe,merhe,mesve,meive
      real(8), dimension(nen,nen) :: mdshe,mdihe,mdrhe,mdsve,mdive
      real(8), dimension(nen,  1) :: fe, bte
      real(8), dimension(2) :: gradi, gradj
      real(8), dimension(4)       :: ee
      real(8), dimension(4,2)     :: normal
      real(8) :: detjf
      ! real(8) :: phi2d,dphi2dx,dphi2dy,f
      real(8) :: f
      real(8) :: ft,ftcont,dftcont
      real(8) :: xil1d
      integer :: ig,jg,ibanda,jbanda,ib,jb
!
      real(8), dimension(nen) :: phi
      real(8), dimension(nen,2) :: grad
      real(8) :: phi2,grad2
      real(8) :: a1dt,a2dt,a3dt,a4dt,a5dt
      real(8) :: ldt,sdmdt,mudt,lvdt,muvdt
!
!     constantes auxiliares
!
      a1dt = alpha(1)*dt
      a2dt = alpha(2)*dt
      a3dt = alpha(3)*dt
      a4dt = alpha(4)*dt
      a5dt = alpha(5)*dt
      ldt = lamb*dt
      sdmdt = (sig+delta+mui)*dt
      mudt = mu*dt
      lvdt = lambv*dt
      muvdt = muv*dt
!
      mesh=0.d0
      meih=0.d0
      merh=0.d0
      mesv=0.d0
      meiv=0.d0
      mdsh=0.d0
      mdih=0.d0
      mdrh=0.d0
      mdsv=0.d0
      mdiv=0.d0
!
      npos=0
      ilin=0
      jcol=0
!
      do nel=1,nelem
!
      meshe=0.d0
      meihe=0.d0
      merhe=0.d0
      mesve=0.d0
      meive=0.d0
      mdshe=0.d0
      mdihe=0.d0
      mdrhe=0.d0
      mdsve=0.d0
      mdive=0.d0
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(nel)%n(j),k)
      end do ! k
      end do ! j
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     real coordinates of the integration point
!
      xxl=0.d0
      do j=1,nen
      do k=1,nsd
      xxl(k)=xxl(k)+phi2d(j,xil)*xe(j,k)
      end do ! k
      end do ! j
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1)= jf(2,2)
      invjf(1,2)=-jf(1,2)
      invjf(2,1)=-jf(2,1)
      invjf(2,2)= jf(1,1)
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
!     funcoes de base e gradientes
!
      do i=1,nen
      phi(i) = phi2d(i,xil)
      grad(i,1) = tinvjf(1,1)*dphi2dx(i,xil)+tinvjf(1,2)*dphi2dy(i,xil)
      grad(i,2) = tinvjf(2,1)*dphi2dx(i,xil)+tinvjf(2,2)*dphi2dy(i,xil)
      end do !i
!
!     matrizes
!
      do i=1,nen
      do j=1,nen
!
      phi2 = phi(i)*phi(j)*detjf*wg(l,nint)
      grad2 = (grad(j,1)*grad(i,1) + grad(j,2)*grad(i,2))*detjf*wg(l,nint)
!
      meshe(i,j)=meshe(i,j) + (1.d0-ldt+mudt)*phi2 + a1dt*grad2
      meihe(i,j)=meihe(i,j) + (1.d0+sdmdt)*phi2 + a2dt*grad2
      merhe(i,j)=merhe(i,j) + (1.d0+mudt)*phi2 + a3dt*grad2
      ! mesve(i,j)=mesve(i,j) + (1.d0-lvdt+muvdt)*phi2 + a4dt*grad2
      mesve(i,j)=mesve(i,j) + (1.d0-lvdt)*phi2 + a4dt*grad2
      meive(i,j)=meive(i,j) + (1.d0+muvdt)*phi2 + a5dt*grad2
!
      mdshe(i,j)=mdshe(i,j) + phi2
      mdihe(i,j)=mdihe(i,j) + phi2
      mdrhe(i,j)=mdrhe(i,j) + phi2
      mdsve(i,j)=mdsve(i,j) + phi2
      mdive(i,j)=mdive(i,j) + phi2
!
      end do ! j
      end do ! i
!
      end do ! l
!
!     espalha a local na global
!
      do i=1,nen
      ig=el(nel)%n(i)
      do j=1,nen
      jg=el(nel)%n(j)
!
      ! ib=ibanda(ig,jg)
      ! jb=jbanda(ig,jg)
!
!     looking for index
!
      ii = 0
!
      do while (ii.le.npos)
!
      ii = ii + 1
!
      if ( (ilin(ii).eq.ig).and.(jcol(ii).eq.jg) ) then
         go to 10
      end if
!
      end do !while
!
      npos = ii
      ilin(ii) = ig
      jcol(ii) = jg
!
10    mesh(ii,1)=mesh(ii,1)+meshe(i,j)
      meih(ii,1)=meih(ii,1)+meihe(i,j)
      merh(ii,1)=merh(ii,1)+merhe(i,j)
      mesv(ii,1)=mesv(ii,1)+mesve(i,j)
      meiv(ii,1)=meiv(ii,1)+meive(i,j)
!
      mdsh(ii,1)=mdsh(ii,1)+mdshe(i,j)
      mdih(ii,1)=mdih(ii,1)+mdihe(i,j)
      mdrh(ii,1)=mdrh(ii,1)+mdrhe(i,j)
      mdsv(ii,1)=mdsv(ii,1)+mdsve(i,j)
      mdiv(ii,1)=mdiv(ii,1)+mdive(i,j)
!
! !     para o lapack
! !
!       mesh(ib,jb)=mesh(ib,jb)+meshe(i,j)
!       meih(ib,jb)=meih(ib,jb)+meihe(i,j)
!       merh(ib,jb)=merh(ib,jb)+merhe(i,j)
!       mesv(ib,jb)=mesv(ib,jb)+mesve(i,j)
!       meiv(ib,jb)=meiv(ib,jb)+meive(i,j)
! !
!       mdsh(ku+1-jg+ig,jb)=mdsh(ku+1-jg+ig,jb)+mdshe(i,j)
!       mdih(ku+1-jg+ig,jb)=mdih(ku+1-jg+ig,jb)+mdihe(i,j)
!       mdrh(ku+1-jg+ig,jb)=mdrh(ku+1-jg+ig,jb)+mdrhe(i,j)
!       mdsv(ku+1-jg+ig,jb)=mdsv(ku+1-jg+ig,jb)+mdsve(i,j)
!       mdiv(ku+1-jg+ig,jb)=mdiv(ku+1-jg+ig,jb)+mdive(i,j)
!
      end do ! j
      end do ! i
!
      end do ! nel
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine assemblenl
!
      use mgeometry
      use mcoeficientes
      use mgauss
      use marrays
!
      implicit none
!
      integer :: nel,i,j,k,l,nk,n1,n2,ii
      real(8), dimension(nsd,nsd) :: jf,invjf,tinvjf

      real(8), dimension(nsd)     :: xil,xxl,xx
      real(8), dimension(nen,nsd) :: xe
      real(8), dimension(nen,nen) :: meshnle,meihnle,mesvnle
      real(8), dimension(nen,1) :: fshe,fihe,frhe,fsve,five
      real(8), dimension(2) :: gradi, gradj
      real(8), dimension(nsd,  2) :: grad
      real(8), dimension(4)       :: ee
      real(8), dimension(4,2)     :: normal
      real(8) :: detjf
      ! real(8) :: phi2d,dphi2dx,dphi2dy,f
      real(8) :: f
      real(8) :: ft,ftcont,dftcont
      real(8) :: xil1d
      integer :: ig,jg,ibanda,jbanda,ib,jb,kg
!
      real(8), dimension(nen) :: phi
      real(8) :: phi1,phi2,phi3
      real(8) :: shk,ihk,rhk,svk,ivk,shj,ihj,rhj,svj,ivj
      real(8) :: ldtk,bdt,l2dtk,lsdt,ldt,udt,ddt,bsdt
      real(8) :: lvdtkv,bvdt,l2vdtkv,lvdt,ukdt
!
      real(8), dimension(5) :: fsource
!
!     constantes auxiliares
!
      ldtk = lamb*dt/kp*ikp
      bdt = beta*dt
      l2dtk = 2.d0*lamb*dt/kp*ikp
      lsdt = (lamb+sig)*dt
      ldt = lamb*dt
      ! udt = u*dt
      ddt = delta*dt
      lvdtkv = lambv*dt/kpv
      bvdt = betav*dt
      l2vdtkv = 2.d0*lambv*dt/kpv
      lvdt = lambv*dt
      bsdt = betas*dt
!
      meshnl=0.d0
      meihnl=0.d0
!      merhnl=0.d0
      mesvnl=0.d0
!      meivnl=0.d0
      fsh=0.d0
      fih=0.d0
      frh=0.d0
      fsv=0.d0
      fiv=0.d0
!
      do nel=1,nelem
!
      meshnle=0.d0
      meihnle=0.d0
      mesvnle=0.d0
      fshe=0.d0
      fihe=0.d0
      frhe=0.d0
      fsve=0.d0
      five=0.d0
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(nel)%n(j),k)
      end do ! k
      end do ! j
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     real coordinates of the integration point
!
      xxl=0.d0
      do j=1,nen
      do k=1,nsd
      xxl(k)=xxl(k)+phi2d(j,xil)*xe(j,k)
      end do ! k
      end do ! j
!
      call jacobf(xe,xil,jf,detjf)
!
      ! call source(fsource,xxl,tt+dt)
!
!     funcoes de base
!
      do i=1,nen
      phi(i) = phi2d(i,xil)
      end do !i
!
!     control at integration points
!
      ukdt = 0.d0
!
      do k=1,nen
!
      kg=el(nel)%n(k)
!
      ! ukdt = ukdt + dt*u(kg,it+1)*phi(k)
      ukdt = ukdt + dt*a*u(kg,it+1)*phi(k)
!
      end do !k
!
!     matrizes
!
      do i=1,nen
!
      phi1 = phi(i)*detjf*wg(l,nint)
!
      do j=1,nen
!
      jg=el(nel)%n(j)
!
      phi2 = phi(j)*phi1
!
      shj = suh(jg,it+1)
      ihj = inh(jg,it+1)
      rhj = reh(jg,it+1)
      svj = suv(jg,it+1)
      ivj = inv(jg,it+1)
!
      do k=1,nen
!
      kg=el(nel)%n(k)
!
      phi3 = phi(k)*phi2
!
      shk = suh(kg,it+1)
      ihk = inh(kg,it+1)
      rhk = reh(kg,it+1)
      svk = suv(kg,it+1)
      ivk = inv(kg,it+1)
!
      meshnle(i,j)=meshnle(i,j)+(ldtk*shk+bdt*ivk+bsdt*ihk+l2dtk*(ihk+rhk))*phi3
      meihnle(i,j)=meihnle(i,j)-bsdt*shk*phi3
      mesvnle(i,j)=mesvnle(i,j)+(lvdtkv*svk+bvdt*ihk+l2vdtkv*ivk)*phi3
      fshe(i,1)=fshe(i,1)+(-ldtk*ihj*ihk-l2dtk*ihj*rhk-ldtk*rhj*rhk)*phi3
      fihe(i,1)=fihe(i,1)+(bdt*shj*ivk)*phi3
      fsve(i,1)=fsve(i,1)+(-lvdtkv*ivj*ivk)*phi3
      five(i,1)=five(i,1)+(bvdt*svj*ihk)*phi3
!
      end do ! k
!
      meshnle(i,j)=meshnle(i,j)+ukdt*phi2
      fshe(i,1)=fshe(i,1) + (lsdt*ihj+ldt*rhj)*phi2
      frhe(i,1)=frhe(i,1) + (ukdt*shj+ddt*ihj)*phi2
      fsve(i,1)=fsve(i,1) + (lvdt*ivj)*phi2
!
      end do ! j
!
      ! fshe(i,1) = fshe(i,1) + dt*fsource(1)*phi1
      ! fihe(i,1) = fihe(i,1) + dt*fsource(2)*phi1
      ! frhe(i,1) = frhe(i,1) + dt*fsource(3)*phi1
      ! fsve(i,1) = fsve(i,1) + dt*fsource(4)*phi1
      ! five(i,1) = five(i,1) + dt*fsource(5)*phi1
!
      end do ! i
!
      end do ! l
!
!     espalha a local na global
!
      do i=1,nen
      ig=el(nel)%n(i)
      fsh(ig,1)=fsh(ig,1)+fshe(i,1)
      fih(ig,1)=fih(ig,1)+fihe(i,1)
      frh(ig,1)=frh(ig,1)+frhe(i,1)
      fsv(ig,1)=fsv(ig,1)+fsve(i,1)
      fiv(ig,1)=fiv(ig,1)+five(i,1)

      do j=1,nen
      jg=el(nel)%n(j)
!
!     looking for index
!
      ii = 0
!
      do while (ii.le.npos)
!
      ii = ii + 1
!
      if ( (ilin(ii).eq.ig).and.(jcol(ii).eq.jg) ) then
         go to 10
      end if
!
      end do !while
!
      npos = ii
      ilin(ii) = ig
      jcol(ii) = jg
!
10    meshnl(ii,1)=meshnl(ii,1)+meshnle(i,j)
      meihnl(ii,1)=meihnl(ii,1)+meihnle(i,j)
      mesvnl(ii,1)=mesvnl(ii,1)+mesvnle(i,j)
!
      ! ib=ibanda(ig,jg)
      ! jb=jbanda(ig,jg)
! !
! !     para o lapack
! !
!       meshnl(ib,jb)=meshnl(ib,jb)+meshnle(i,j)
!       mesvnl(ib,jb)=mesvnl(ib,jb)+mesvnle(i,j)
!
      end do ! j
      end do ! i
!
      end do ! nel
!
!     adding new cases in cities
!
!       if ((tt>=14.d0-dt).and.(tt<=21.d0-dt)) then!8.d0-dt)) then
! !
!       do i=1,3
! !
!       ig = el(156)%n(i)
! ! 
!       fih(ig,1)=fih(ig,1)+dt*i02*1.d-3!2.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
!       fsh(ig,1)=fsh(ig,1)-dt*i02*1.d-3!2.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
! !
!       end do !i
! !
!       end if

!       if ((tt>=28.d0-dt).and.(tt<=35.d0-dt)) then
! !
!       do i=1,3
! !
!       ig = el(156)%n(i)
! ! 
!       fih(ig,1)=fih(ig,1)+dt*3.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
!       fsh(ig,1)=fsh(ig,1)-dt*3.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
! !
!       end do !i
! !
!       end if
! !
!       if ((tt>=35.d0-dt).and.(tt<=42.d0-dt)) then
! !
!       do i=1,3
! !
!       ig = el(156)%n(i)
! ! 
!       fih(ig,1)=fih(ig,1)+dt*7.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
!       fsh(ig,1)=fsh(ig,1)-dt*7.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
! !
!       end do !i
! !
!       end if

!       if ((tt>=42.d0-dt).and.(tt<=49.d0-dt)) then
! !
!       do i=1,3
! !
!       ig = el(156)%n(i)
! ! 
!       fih(ig,1)=fih(ig,1)+dt*6.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
!       fsh(ig,1)=fsh(ig,1)-dt*6.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
! !
!       end do !i
! !
!       end if

!       if ((tt>=35.d0-dt).and.(tt<=42.d0-dt)) then
! !
!       do i=1,3
! !
!       ig = el(247)%n(i)
! ! 
!       fih(ig,1)=fih(ig,1)+dt*i03*1.d-3!2.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/4.8986d-3/10195921.63d0/7.d0/0.2d0
!       fsh(ig,1)=fsh(ig,1)-dt*i03*1.d-3!2.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/4.8986d-3/10195921.63d0/7.d0/0.2d0
! !
!       end do !i

!       end if
!
!       if ((tt>=49.d0-dt).and.(tt<=63.d0-dt)) then
! !
!       do i=1,3
! !
!       ig = el(247)%n(i)
! ! 
!       fih(ig,1)=fih(ig,1)+dt*1.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/4.8986d-3/10195921.63d0/7.d0/0.2d0
!       fsh(ig,1)=fsh(ig,1)-dt*1.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/4.8986d-3/10195921.63d0/7.d0/0.2d0
! !
!       end do !i
! !
!       end if

!       if ((tt>=63.d0-dt).and.(tt<=70.d0-dt)) then
! !
!       do i=1,3
! !
!       ig = el(247)%n(i)
! ! 
!       fih(ig,1)=fih(ig,1)+dt*87.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/4.8986d-3/10195921.63d0/7.d0/0.2d0
!       fsh(ig,1)=fsh(ig,1)-dt*87.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/4.8986d-3/10195921.63d0/7.d0/0.2d0
! !
!       end do !i
! !
      ! end if



      if ((tt>=7.d0-dt).and.(tt<=14.d0-dt)) then!8.d0-dt)) then
!
      do i=1,3
!
      ig = el(156)%n(i)
! 
      fih(ig,1)=fih(ig,1)+dt*i02*1.d-4!2.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
      fsh(ig,1)=fsh(ig,1)-dt*i02*1.d-4!2.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/5.5118d-3/10195921.63d0/7.d0/0.2d0
      ! fiv(ig,1)=fiv(ig,1)+dt*i02*1.d-4*5.d0
      ! fsv(ig,1)=fsv(ig,1)-dt*i02*1.d-4*5.d0
!
      end do !i
!
      end if

      if ((tt>=28.d0-dt).and.(tt<=35.d0-dt)) then
!
      do i=1,3
!
      ig = el(247)%n(i)
! 
      fih(ig,1)=fih(ig,1)+dt*i03*1.d-4!2.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/4.8986d-3/10195921.63d0/7.d0/0.2d0
      fsh(ig,1)=fsh(ig,1)-dt*i03*1.d-4!2.d0/10195921.63d0/7.d0/3.d0/0.2d0!2.d0/4.8986d-3/10195921.63d0/7.d0/0.2d0
      ! fiv(ig,1)=fiv(ig,1)+dt*i03*1.d-4*5.d0
      ! fsv(ig,1)=fsv(ig,1)-dt*i03*1.d-4*5.d0
!
      end do !i

      end if
!
      end subroutine
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      subroutine assemblead
!
      use mgeometry
      use mcoeficientes
      use mgauss
      use marrays
!
      implicit none
!
      integer :: nel,i,j,k,l,nk,n1,n2,ii
      real(8), dimension(nsd,nsd) :: jf,invjf,tinvjf

      real(8), dimension(nsd)     :: xil,xxl,xx
      real(8), dimension(nen,nsd) :: xe
      real(8), dimension(nen,nen) :: mel1e,mel2e,mel3e,mel4e,mel5e
      real(8), dimension(nen,nen) :: mdl1e,mdl2e,mdl3e,mdl4e,mdl5e
      real(8), dimension(2) :: gradi, gradj
      real(8), dimension(4)       :: ee
      real(8), dimension(4,2)     :: normal
      real(8) :: detjf
      ! real(8) :: phi2d,dphi2dx,dphi2dy,f
      real(8) :: f
      real(8) :: xil1d
      ! integer :: ig,jg,kg,ibanda,jbanda,ib,jb
      integer :: ig,jg,kg,ib,jb
!
      real(8), dimension(nen) :: phi
      real(8), dimension(nen,2) :: grad
      real(8) :: phi2,grad2,phi3
      real(8) :: a1dt,a2dt,a3dt,a4dt,a5dt
      real(8) :: bdt,ldt,l2dtk,sdmdt,mudt,bvdt,lvdt,lv2dtkv,muvdt,udt
      real(8) :: shk,ihk,rhk,svk,ivk,ukdt,bsdt
!
!     auxiliary constants
!
      a1dt = alpha(1)*dt
      a2dt = alpha(2)*dt
      a3dt = alpha(3)*dt
      a4dt = alpha(4)*dt
      a5dt = alpha(5)*dt
!
      bdt = beta*dt
      ldt = lamb*dt
      l2dtk = 2.d0*lamb*dt/kp*ikp
      ! udt = u*dt
      sdmdt = (sig+delta+mui)*dt
      mudt = mu*dt
      bvdt = betav*dt
      lvdt = lambv*dt
      lv2dtkv = 2.d0*lambv*dt/kpv
      muvdt = muv*dt
      bsdt = betas*dt
!
      mel1=0.d0
      mel2=0.d0
      mel3=0.d0
      mel4=0.d0
      mel5=0.d0
      mdl1=0.d0
      mdl2=0.d0
      mdl3=0.d0
      mdl4=0.d0
      mdl5=0.d0
!
      do nel=1,nelem
!
      mel1e=0.d0
      mel2e=0.d0
      mel3e=0.d0
      mel4e=0.d0
      mel5e=0.d0
      mdl1e=0.d0
      mdl2e=0.d0
      mdl3e=0.d0
      mdl4e=0.d0
      mdl5e=0.d0
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(nel)%n(j),k)
      end do ! k
      end do ! j
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     real coordinates of the integration point
!
      xxl=0.d0
      do j=1,nen
      do k=1,nsd
      xxl(k)=xxl(k)+phi2d(j,xil)*xe(j,k)
      end do ! k
      end do ! j
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1)= jf(2,2)
      invjf(1,2)=-jf(1,2)
      invjf(2,1)=-jf(2,1)
      invjf(2,2)= jf(1,1)
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
!     basis functions and gradients
!
      do i=1,nen
      phi(i) = phi2d(i,xil)
      grad(i,1) = tinvjf(1,1)*dphi2dx(i,xil)+tinvjf(1,2)*dphi2dy(i,xil)
      grad(i,2) = tinvjf(2,1)*dphi2dx(i,xil)+tinvjf(2,2)*dphi2dy(i,xil)
      end do !i
!
!     state solutions and control at integration points
!
      shk = 0.d0
      ihk = 0.d0
      rhk = 0.d0
      svk = 0.d0
      ivk = 0.d0
      ukdt = 0.d0
!
      do k=1,nen
!
      kg=el(nel)%n(k)
!
      shk = shk + suh(kg,itt-1)*phi(k)
      ihk = ihk + inh(kg,itt-1)*phi(k)
      rhk = rhk + reh(kg,itt-1)*phi(k)
      svk = svk + suv(kg,itt-1)*phi(k)
      ivk = ivk + inv(kg,itt-1)*phi(k)
      ukdt = ukdt + dt*u(kg,itt-1)*phi(k)
      ! ukdt = ukdt + dt*a*u(kg,itt-1)*phi(k)
!
      end do !k
!
!     matrices
!
      do i=1,nen
      do j=1,nen
!
      phi2 = phi(i)*phi(j)*detjf*wg(l,nint)
      grad2 = (grad(j,1)*grad(i,1) + grad(j,2)*grad(i,2))*detjf*wg(l,nint)
!
      mel1e(i,j)=mel1e(i,j) + a1dt*grad2 + (1.d0+bdt*ivk+bsdt*ihk-ldt+mudt &
         +l2dtk*(shk+ihk+rhk)+ukdt)*phi2
      mel2e(i,j)=mel2e(i,j) + a2dt*grad2 + (1.d0+sdmdt-bsdt*ihk)*phi2
      mel3e(i,j)=mel3e(i,j) + a3dt*grad2 + (1.d0+mudt)*phi2
      mel4e(i,j)=mel4e(i,j) + a4dt*grad2 + (1.d0+bvdt*ihk-lvdt+muvdt &
         +lv2dtkv*(svk+ivk))*phi2
      mel5e(i,j)=mel5e(i,j) + a5dt*grad2 + (1.d0+muvdt)*phi2
!
      mdl1e(i,j)=mdl1e(i,j) + phi2
      mdl2e(i,j)=mdl2e(i,j) + phi2
      mdl3e(i,j)=mdl3e(i,j) + phi2
      mdl4e(i,j)=mdl4e(i,j) + phi2
      mdl5e(i,j)=mdl5e(i,j) + phi2
!
      end do ! j
      end do ! i
!
      end do ! l
!
!     scatter local in global matrices
!
      do i=1,nen
      ig=el(nel)%n(i)
      do j=1,nen
      jg=el(nel)%n(j)
!
!     looking for index
!
      ii = 0
!
      do while (ii.le.npos)
!
      ii = ii + 1
!
      if ( (ilin(ii).eq.ig).and.(jcol(ii).eq.jg) ) then
         go to 10
      end if
!
      end do !while
!
      npos = ii
      ilin(ii) = ig
      jcol(ii) = jg
!
10    mel1(ii,1)=mel1(ii,1)+mel1e(i,j)
      mel2(ii,1)=mel2(ii,1)+mel2e(i,j)
      mel3(ii,1)=mel3(ii,1)+mel3e(i,j)
      mel4(ii,1)=mel4(ii,1)+mel4e(i,j)
      mel5(ii,1)=mel5(ii,1)+mel5e(i,j)
!
      mdl1(ii,1)=mdl1(ii,1)+mdl1e(i,j)
      mdl2(ii,1)=mdl2(ii,1)+mdl2e(i,j)
      mdl3(ii,1)=mdl3(ii,1)+mdl3e(i,j)
      mdl4(ii,1)=mdl4(ii,1)+mdl4e(i,j)
      mdl5(ii,1)=mdl5(ii,1)+mdl5e(i,j)
!
      ! ib=ibanda(ig,jg)
      ! jb=jbanda(ig,jg)
!
!     for lapack
!
!       mel1(ib,jb)=mel1(ib,jb)+mel1e(i,j)
!       mel2(ib,jb)=mel2(ib,jb)+mel2e(i,j)
!       mel3(ib,jb)=mel3(ib,jb)+mel3e(i,j)
!       mel4(ib,jb)=mel4(ib,jb)+mel4e(i,j)
!       mel5(ib,jb)=mel5(ib,jb)+mel5e(i,j)
! !
!       mdl1(ku+1-jg+ig,jb)=mdl1(ku+1-jg+ig,jb)+mdl1e(i,j)
!       mdl2(ku+1-jg+ig,jb)=mdl2(ku+1-jg+ig,jb)+mdl2e(i,j)
!       mdl3(ku+1-jg+ig,jb)=mdl3(ku+1-jg+ig,jb)+mdl3e(i,j)
!       mdl4(ku+1-jg+ig,jb)=mdl4(ku+1-jg+ig,jb)+mdl4e(i,j)
!       mdl5(ku+1-jg+ig,jb)=mdl5(ku+1-jg+ig,jb)+mdl5e(i,j)
!
      end do ! j
      end do ! i
!
      end do ! nel
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine assemblenlad
!
      use mgeometry
      use mcoeficientes
      use mgauss
      use marrays
!
      implicit none
!
      integer :: nel,i,j,k,l,nk,n1,n2
      real(8), dimension(nsd,nsd) :: jf,invjf,tinvjf

      real(8), dimension(nsd)     :: xil,xxl,xx
      real(8), dimension(nen,nsd) :: xe
      real(8), dimension(nen,1) :: fl1e,fl2e,fl3e,fl4e,fl5e
      real(8), dimension(2) :: gradi, gradj
      real(8), dimension(nsd,  2) :: grad
      real(8), dimension(4)       :: ee
      real(8), dimension(4,2)     :: normal
      real(8) :: detjf
      ! real(8) :: phi2d,dphi2dx,dphi2dy,f
      real(8) :: f
      real(8) :: ft,ftcont,dftcont
      real(8) :: xil1d
      integer :: ig,jg,ibanda,jbanda,ib,jb,kg
!
      real(8), dimension(nen) :: phi
      real(8) :: phi1,phi2,phi3
      real(8) :: l1j,l2j,l3j,l4j,l5j
      real(8) :: bdt,udt,c1udt,ldt,lsdt,l2dtk,ddt,bvdt,lvdt,lv2dtkv
      real(8) :: shk,ihk,rhk,svk,ivk,ukdt,bsdt
!
      real(8), dimension(5) :: fsource
!
!     auxiliary constants
!
      bdt = beta*dt
      ! udt = u*dt
      ! c1udt = c1*u*dt
      ldt = lamb*dt
      lsdt = (lamb+sig)*dt
      l2dtk = 2.d0*lamb*dt/kp*ikp
      ddt = delta*dt
      bvdt = betav*dt
      lvdt = lambv*dt
      lv2dtkv = 2.d0*lambv*dt/kpv
      bsdt = betas*dt
!
      fl1=0.d0
      fl2=0.d0
      fl3=0.d0
      fl4=0.d0
      fl5=0.d0
!
      do nel=1,nelem
!
      fl1e=0.d0
      fl2e=0.d0
      fl3e=0.d0
      fl4e=0.d0
      fl5e=0.d0
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(nel)%n(j),k)
      end do ! k
      end do ! j
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     real coordinates of the integration point
!
      xxl=0.d0
      do j=1,nen
      do k=1,nsd
      xxl(k)=xxl(k)+phi2d(j,xil)*xe(j,k)
      end do ! k
      end do ! j
!
      call jacobf(xe,xil,jf,detjf)
!
      ! call sourcead(fsource,xxl,tt-dt)
!
!     basis functions
!
      do i=1,nen
      phi(i) = phi2d(i,xil)
      end do !i
!
!     state solutions and control at integration points
!
      shk = 0.d0
      ihk = 0.d0
      rhk = 0.d0
      svk = 0.d0
      ivk = 0.d0
      ukdt = 0.d0
!
      do k=1,nen
!
      kg=el(nel)%n(k)
!
      shk = shk + suh(kg,itt-1)*phi(k)
      ihk = ihk + inh(kg,itt-1)*phi(k)
      rhk = rhk + reh(kg,itt-1)*phi(k)
      svk = svk + suv(kg,itt-1)*phi(k)
      ivk = ivk + inv(kg,itt-1)*phi(k)
      ukdt = ukdt + dt*u(kg,itt-1)*phi(k)
      ! ukdt = ukdt + dt*a*u(kg,itt-1)*phi(k)
!
      end do !k
!
!     matrices
!
      do i=1,nen
!
      phi1 = phi(i)*detjf*wg(l,nint)
!
      do j=1,nen
!
      jg=el(nel)%n(j)
!
      phi2 = phi(j)*phi1
!
      l1j = l1(jg,itt-1)
      l2j = l2(jg,itt-1)
      l3j = l3(jg,itt-1)
      l4j = l4(jg,itt-1)
      l5j = l5(jg,itt-1)
!
      fl1e(i,1) = fl1e(i,1) + ((bdt*ivk+bsdt*ihk)*l2j+ukdt*l3j)*phi2
      fl2e(i,1) = fl2e(i,1) + ((lsdt-l2dtk*(shk+ihk+rhk)-bsdt*shk)*l1j&
         +ddt*l3j-bvdt*svk*l4j+bvdt*svk*l5j)*phi2
      fl3e(i,1) = fl3e(i,1) + ((ldt-l2dtk*(shk+ihk+rhk))*l1j)*phi2
      fl4e(i,1) = fl4e(i,1) + (bvdt*ihk*l5j)*phi2
      fl5e(i,1) = fl5e(i,1) + (-bdt*shk*l1j+bdt*shk*l2j&
         +(lvdt-lv2dtkv*(svk+ivk))*l4j)*phi2
!
      end do ! j
!
      fl1e(i,1) = fl1e(i,1) + c2*ukdt*phi1
      ! fl1e(i,1) = fl1e(i,1) + c2*ukdt*phi1/a
      fl2e(i,1) = fl2e(i,1) + c1*dt*phi1
!
      ! fl1e(i,1) = fl1e(i,1) - dt*fsource(1)*phi1
      ! fl2e(i,1) = fl2e(i,1) - dt*fsource(2)*phi1
      ! fl3e(i,1) = fl3e(i,1) - dt*fsource(3)*phi1
      ! fl4e(i,1) = fl4e(i,1) - dt*fsource(4)*phi1
      ! fl5e(i,1) = fl5e(i,1) - dt*fsource(5)*phi1
!
      end do ! i
!
      end do ! l
!
!     scatter local in global matrices
!
      do i=1,nen
!
      ig=el(nel)%n(i)
!
      fl1(ig,1)=fl1(ig,1)+fl1e(i,1)
      fl2(ig,1)=fl2(ig,1)+fl2e(i,1)
      fl3(ig,1)=fl3(ig,1)+fl3e(i,1)
      fl4(ig,1)=fl4(ig,1)+fl4e(i,1)
      fl5(ig,1)=fl5(ig,1)+fl5e(i,1)
!
      end do ! i
!
      end do ! nel
!
      end subroutine
!
!-----------------------------------------------------------------------
!
      function phi(i,nen,xi)
!
!     one dimensional lagrangian functions
!
      implicit none
!
      real(8) :: phi,xi
      integer :: i,nen
!
      select case(nen)
!
      case(1)
!
!     constant
!
      phi=1.d0
!
      case(2)
!
!     linear
!
      if(i.eq.1) then
      phi=0.5d0*(1.d0-xi)
      else
      if(i.eq.2) then
      phi=0.5d0*(1.d0+xi)
      else
      write(*,*) 'inexistent node'
      stop
      end if
      end if
!
      case(3)
!
!     quadratic
!
      if(i.eq.1) then
      phi= 0.5d0*xi*(xi-1.d0)
      else
      if(i.eq.2) then
      phi=(1.d0-xi)*(1.d0+xi)
      else
      if(i.eq.3) then
      phi= 0.5d0*xi*(xi+1.d0)
      else
      write(*,*) 'inexistent node'
      stop
      end if
      end if
      end if
!
      case(4)
!
!     cubic
!
      if(i.eq.1) then
      phi= (xi+1.d0/3.d0)*(xi-1.d0/3.d0)*(xi-1.d0)*(-9.d0/16.d0)
      else
      if(i.eq.2) then
      phi= (xi+1.d0)*(xi-1.d0/3.d0)*(xi-1.d0)*(27.d0/16.d0)
      else
      if(i.eq.3) then
      phi= (xi+1.d0)*(xi+1.d0/3.d0)*(xi-1.d0)*(-27.d0/16.d0)
      else
      if(i.eq.4) then
      phi= (xi+1.d0)*(xi+1.d0/3.d0)*(xi-1.d0/3.d0)*(9.d0/16.d0)
      else
      write(*,*) 'inexistent node'
      stop
      end if
      end if
      end if
      end if
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function phi
!
!-----------------------------------------------------------------------
!
      function dphi(i,nen,xi)
!
!     derivative of the one dimensional lagrangian functions
!
      implicit none
!
      real(8) :: dphi,xi
      integer :: i,nen
!
      select case(nen)
!
      case(1)
!
!     constant
!
      dphi=0.d0
!
      case(2)
!
!     linear
!
      if(i.eq.1) then
      dphi=-0.5d0
      else
      dphi= 0.5d0
      end if
!
      case(3)
!
!     quadratic
!
      if(i.eq.1) then
      dphi= xi-0.5d0
      else
      if(i.eq.2) then
      dphi=-2.d0*xi
      else
      dphi= xi+0.5d0
      end if
      end if
!
      case(4)
!
!     cubic
!
      if(i.eq.1) then
      dphi= (-27.d0*xi**2+18.d0*xi+1.d0)/16.d0
      else
      if(i.eq.2) then
      dphi= (3.d0*xi**2-2.d0/3.d0*xi-1.d0)*(27.d0/16.d0)
      else
      if(i.eq.3) then
      dphi= (3.d0*xi**2+2.d0/3.d0*xi-1.d0)*(-27.d0/16.d0)
      else
      dphi= (27.d0*xi**2+18.d0*xi-1.d0)/16.d0
      end if
      end if
      end if
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function dphi

!
!-----------------------------------------------------------------------
!
      function phi2d(i,xi)
      use mgeometry, only : nsd,nen
      implicit none
!
      real(8) :: phi2d
      real(8) :: phi
      real(8), dimension(nsd) :: xi
      integer :: i,n1d
!
      select case(nen)
!
!........................
!
      case(1)
!
!     constant
!
      phi2d=1.d0
!
!........................
!
      case(3)
!
!     linear - triangles
!
      n1d=2
!
      select case(i)
!
      case(1)
      phi2d=1.d0-xi(1)-xi(2)
!
      case(2)
      phi2d=xi(1)
!
      case(3)
      phi2d=xi(2)
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(4)
!
!     bilinear
!
      n1d=2
!
      select case(i)
!
      case(1)
      phi2d=phi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      phi2d=phi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      phi2d=phi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(4)
      phi2d=phi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!

      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(9)
!
!     biquadratic
!
      n1d=3
!
      select case(i)
!
      case(1)
      phi2d=phi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      phi2d=phi(3,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      phi2d=phi(3,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(4)
      phi2d=phi(1,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(5)
      phi2d=phi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(6)
      phi2d=phi(3,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(7)
      phi2d=phi(2,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(8)
      phi2d=phi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(9)
      phi2d=phi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(16)
!
!     bicubic
!
      n1d=4
!
      select case(i)
!
      case(1)
      phi2d=phi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      phi2d=phi(4,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      phi2d=phi(4,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(4)
      phi2d=phi(1,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(5)
      phi2d=phi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(6)
      phi2d=phi(3,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(7)
      phi2d=phi(4,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(8)
      phi2d=phi(4,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(9)
      phi2d=phi(3,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(10)
      phi2d=phi(2,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(11)
      phi2d=phi(1,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(12)
      phi2d=phi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(13)
      phi2d=phi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(14)
      phi2d=phi(3,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(15)
      phi2d=phi(3,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(16)
      phi2d=phi(2,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(8)
!
!     serendipity 8
!
      n1d=4
!
      select case(i)
!
      case(1)
      phi2d=-1.d0/4.d0*(1.d0-xi(1))*(1.d0-xi(2))*(1.d0+xi(1)+xi(2))
!
      case(2)
      phi2d=-1.d0/4.d0*(1.d0+xi(1))*(1.d0-xi(2))*(1.d0-xi(1)+xi(2))
!
      case(3)
      phi2d=-1.d0/4.d0*(1.d0+xi(1))*(1.d0+xi(2))*(1.d0-xi(1)-xi(2))
!
      case(4)
      phi2d=-1.d0/4.d0*(1.d0-xi(1))*(1.d0+xi(2))*(1.d0+xi(1)-xi(2))
!
      case(5)
      phi2d=1.d0/2.d0*(1.d0-xi(1))*(1.d0+xi(1))*(1.d0-xi(2))
!
      case(6)
      phi2d=1.d0/2.d0*(1.d0+xi(1))*(1.d0+xi(2))*(1.d0-xi(2))
!
      case(7)
      phi2d=1.d0/2.d0*(1.d0-xi(1))*(1.d0+xi(1))*(1.d0+xi(2))
!
      case(8)
      phi2d=1.d0/2.d0*(1.d0-xi(1))*(1.d0+xi(2))*(1.d0-xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function phi2d
!
!
!-----------------------------------------------------------------------
!
      function dphi2dx(i,xi)
      use mgeometry, only : nsd,nen
      implicit none
!
      real(8) :: dphi2dx
      real(8) :: phi,dphi
      real(8), dimension(nsd) :: xi
      integer :: i,n1d
!
      select case(nen)
!
!........................
!
      case(1)
!
!     constant
!
      dphi2dx=0.d0
!
!........................
!
      case(3)
!
!     linear - triangles
!
      n1d=2
!
      select case(i)
!
      case(1)
      dphi2dx=-1.d0
!
      case(2)
      dphi2dx=1.d0
!
      case(3)
      dphi2dx=0.d0
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(4)
!
!     bilinear
!
      n1d=2
!
      select case(i)
!
      case(1)
      dphi2dx=dphi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      dphi2dx=dphi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      dphi2dx=dphi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(4)
      dphi2dx=dphi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(9)
!
!     biquadratic
!
      n1d=3
!
      select case(i)
!
      case(1)
      dphi2dx=dphi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      dphi2dx=dphi(3,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      dphi2dx=dphi(3,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(4)
      dphi2dx=dphi(1,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(5)
      dphi2dx=dphi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(6)
      dphi2dx=dphi(3,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(7)
      dphi2dx=dphi(2,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(8)
      dphi2dx=dphi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(9)
      dphi2dx=dphi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(16)
!
!     bicubic
!
      n1d=4
!
      select case(i)
!
      case(1)
      dphi2dx=dphi(1,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(2)
      dphi2dx=dphi(4,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(3)
      dphi2dx=dphi(4,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(4)
      dphi2dx=dphi(1,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(5)
      dphi2dx=dphi(2,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(6)
      dphi2dx=dphi(3,n1d,xi(1))*phi(1,n1d,xi(2))
!
      case(7)
      dphi2dx=dphi(4,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(8)
      dphi2dx=dphi(4,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(9)
      dphi2dx=dphi(3,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(10)
      dphi2dx=dphi(2,n1d,xi(1))*phi(4,n1d,xi(2))
!
      case(11)
      dphi2dx=dphi(1,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(12)
      dphi2dx=dphi(1,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(13)
      dphi2dx=dphi(2,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(14)
      dphi2dx=dphi(3,n1d,xi(1))*phi(2,n1d,xi(2))
!
      case(15)
      dphi2dx=dphi(3,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case(16)
      dphi2dx=dphi(2,n1d,xi(1))*phi(3,n1d,xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(8)
!
!     serendipity 8
!
      n1d=4
!
      select case(i)
!
      case(1)
      dphi2dx=-1.d0/4.d0*(-1.d0+xi(2))*(2.d0*xi(1)+xi(2))
!
      case(2)
      dphi2dx=1.d0/4.d0*(-1.d0+xi(2))*(xi(2)-2.d0*xi(1))
!
      case(3)
      dphi2dx=1.d0/4.d0*(1.d0+xi(2))*(2.d0*xi(1)+xi(2))
!
      case(4)
      dphi2dx=-1.d0/4.d0*(1.d0+xi(2))*(xi(2)-2.d0*xi(1))
!
      case(5)
      dphi2dx=xi(1)*(-1.d0+xi(2))
!
      case(6)
      dphi2dx=-1.d0/2.d0*(1.d0+xi(2))*(-1.d0+xi(2))
!
      case(7)
      dphi2dx=-xi(1)*(1.d0+xi(2))
!
      case(8)
      dphi2dx=1.d0/2.d0*(1.d0+xi(2))*(-1.d0+xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function dphi2dx
!
!-----------------------------------------------------------------------
!
      function dphi2dy(i,xi)
      use mgeometry, only : nsd,nen
      implicit none
!
      real(8) :: dphi2dy
      real(8) :: phi,dphi
      real(8), dimension(nsd) :: xi
      integer :: i,n1d
!
      select case(nen)
!
!........................
!
      case(1)
!
!     constant
!
      dphi2dy=0.d0
!
!........................
!
      case(3)
!
!     bilinear
!
      n1d=2
!
      select case(i)
!
      case(1)
      dphi2dy=-1.d0
!
      case(2)
      dphi2dy=0.d0
!
      case(3)
      dphi2dy=1.d0
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(4)
!
!     bilinear
!
      n1d=2
!
      select case(i)
!
      case(1)
      dphi2dy=phi(1,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(2)
      dphi2dy=phi(2,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(3)
      dphi2dy=phi(2,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(4)
      dphi2dy=phi(1,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(9)
!
!     biquadratic
!
      n1d=3
!
      select case(i)
!
      case(1)
      dphi2dy=phi(1,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(2)
      dphi2dy=phi(3,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(3)
      dphi2dy=phi(3,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(4)
      dphi2dy=phi(1,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(5)
      dphi2dy=phi(2,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(6)
      dphi2dy=phi(3,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(7)
      dphi2dy=phi(2,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(8)
      dphi2dy=phi(1,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(9)
      dphi2dy=phi(2,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(16)
!
!     bicubic
!
      n1d=4
!
      select case(i)
!
      case(1)
      dphi2dy=phi(1,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(2)
      dphi2dy=phi(4,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(3)
      dphi2dy=phi(4,n1d,xi(1))*dphi(4,n1d,xi(2))
!
      case(4)
      dphi2dy=phi(1,n1d,xi(1))*dphi(4,n1d,xi(2))
!
      case(5)
      dphi2dy=phi(2,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(6)
      dphi2dy=phi(3,n1d,xi(1))*dphi(1,n1d,xi(2))
!
      case(7)
      dphi2dy=phi(4,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(8)
      dphi2dy=phi(4,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(9)
      dphi2dy=phi(3,n1d,xi(1))*dphi(4,n1d,xi(2))
!
      case(10)
      dphi2dy=phi(2,n1d,xi(1))*dphi(4,n1d,xi(2))
!
      case(11)
      dphi2dy=phi(1,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(12)
      dphi2dy=phi(1,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(13)
      dphi2dy=phi(2,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(14)
      dphi2dy=phi(3,n1d,xi(1))*dphi(2,n1d,xi(2))
!
      case(15)
      dphi2dy=phi(3,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case(16)
      dphi2dy=phi(2,n1d,xi(1))*dphi(3,n1d,xi(2))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
!........................
!
      case(8)
!
!     serendipity 8
!
      n1d=4
!
      select case(i)
!
      case(1)
      dphi2dy=-1.d0/4.d0*(-1.d0+xi(1))*(2.d0*xi(2)+xi(1))
!
      case(2)
      dphi2dy=1.d0/4.d0*(1.d0+xi(1))*(2.d0*xi(2)-xi(1))
!
      case(3)
      dphi2dy=1.d0/4.d0*(1.d0+xi(1))*(2.d0*xi(2)+xi(1))
!
      case(4)
      dphi2dy=-1.d0/4.d0*(-1.d0+xi(1))*(2.d0*xi(2)-xi(1))
!
      case(5)
      dphi2dy=1.d0/2.d0*(1.d0+xi(1))*(-1+xi(1))
!
      case(6)
      dphi2dy=-xi(2)*(1.d0+xi(1))
!
      case(7)
      dphi2dy=-1.d0/2.d0*(1.d0+xi(1))*(-1.d0+xi(1))
!
      case(8)
      dphi2dy=xi(2)*(-1.d0+xi(1))
!
      case default
!
      write(*,*) 'node not defined'
      stop
!
      end select
!
      case default
!
      write(*,*) 'degree not defined'
      stop
!
      end select
!
      end function dphi2dy
!
!-----------------------------------------------------------------------
!
      subroutine jacobf(xe,xi,jf,detjf)
!
!     jacobian of the bilinear transformation
!
      use mgeometry
      implicit none
      real(8), dimension(nsd,nsd) :: jf
      real(8), dimension(nen,nsd) :: xe
      real(8), dimension(nsd)     :: xi
      real(8) :: detjf,eps=1.d-6
!
      ! real(8) :: dphi2dx,dphi2dy
      integer :: i
!
      jf=0.d0
!
      do i=1,nen
!
      jf(1,1)=jf(1,1)+xe(i,1)*dphi2dx(i,xi)
      jf(1,2)=jf(1,2)+xe(i,1)*dphi2dy(i,xi)
!
      jf(2,1)=jf(2,1)+xe(i,2)*dphi2dx(i,xi)
      jf(2,2)=jf(2,2)+xe(i,2)*dphi2dy(i,xi)
!
!
!      write(*,*) i
!
!      write(*,*) dphi2dx(i,xi),dphi2dy(i,xi)
!
      end do ! i
!
!      write(*,"(5e13.4)") jf(1,1),jf(1,2),jf(2,1),jf(2,2)

      detjf=jf(1,1)*jf(2,2)-jf(1,2)*jf(2,1)
!      write(*,*) 'detjf=',detjf
!
!      detjf=dabs(detjf)
!
      if(dabs(detjf).le.eps) then
      write(*,*) 'degenerated element'
      stop
      end if
!
      end subroutine
!
!----------------------------------------------------------------------
!!
!      subroutine prtgnuplot
!!
!      use mgeometry
!      use marrays
!!
!      implicit none
!!
!      integer :: i,j,k
!      real(8) :: zero=0.d0
!!
!      open(unit = 10,file = 'elements.dat',status = 'unknown')
!      open(unit = 20,file = 'solution.dat',status = 'unknown')
!!
!      do i=1,nelem
!!
!      do j=1,4!nen
!!
!      write(10,"(10(f25.15,2x))") coord(el(i)%n(j),1), coord(el(i)%n(j),2), zero
!      write(20,"(10(f25.15,2x))") coord(el(i)%n(j),1), coord(el(i)%n(j),2), p(i,j)
!!
!      end do ! j
!!
!      write(10,"(10(f25.15,2x))") coord(el(i)%n(1),1), coord(el(i)%n(1),2), zero
!      write(20,"(10(f25.15,2x))") coord(el(i)%n(1),1), coord(el(i)%n(1),2), p(i,1)
!!
!      write(10,*)
!      write(10,*)
!      write(20,*)
!      write(20,*)
!!
!      end do ! i
!!
!      close(unit = 10)
!      close(unit = 20)
!!
!      open(unit = 30,file = 'solution2.dat',status = 'unknown')
!      do i = 1, nnodes
!         write(30,"(10(f25.15,2x))") bb(i,1)
!      end do
!      close(unit = 30)
!!
!      end subroutine
!!----------------------------------------------------------------------
!!
!      subroutine postprocessing
!!
!      use mgeometry
!      use marrays
!!
!      implicit none
!!
!      real(8), dimension(nsd)     :: xil,xx
!      real(8), dimension(nen,nsd) :: xe,xer
!!
!      real(8) :: dphi2d,phi2d
!      integer :: i,j,k,n,ii,jj
!!
!!      write(*,*) ndof
!!      do i = 1,ndof
!!      write(*,*) (b(i,j),j=1,ndof)
!!      end do
!!
!!     coordinates of the nodes in the reference element
!!
!      select case(nen)
!!
!      case(4)
!!
!!     bilinear
!!
!      xer(1,1) = -1.d0
!      xer(1,2) = -1.d0
!      xer(2,1) =  1.d0
!      xer(2,2) = -1.d0
!      xer(3,1) =  1.d0
!      xer(3,2) =  1.d0
!      xer(4,1) = -1.d0
!      xer(4,2) =  1.d0
!!
!      case(9)
!!
!!     biquadratic
!!
!      xer(1,1) = -1.d0
!      xer(1,2) = -1.d0
!      xer(2,1) =  1.d0
!      xer(2,2) = -1.d0
!      xer(3,1) =  1.d0
!      xer(3,2) =  1.d0
!      xer(4,1) = -1.d0
!      xer(4,2) =  1.d0
!!
!      xer(5,1) =  0.d0
!      xer(5,2) = -1.d0
!      xer(6,1) =  1.d0
!      xer(6,2) =  0.d0
!      xer(7,1) =  0.d0
!      xer(7,2) =  1.d0
!      xer(8,1) = -1.d0
!      xer(8,2) =  0.d0
!      xer(9,1) =  0.d0
!      xer(9,2) =  0.d0
!!
!      case(16)
!!
!!     bicubic
!!
!      xer(1,1) = -1.d0
!      xer(1,2) = -1.d0
!      xer(2,1) =  1.d0
!      xer(2,2) = -1.d0
!      xer(3,1) =  1.d0
!      xer(3,2) =  1.d0
!      xer(4,1) = -1.d0
!      xer(4,2) =  1.d0
!!
!      xer(5,1) = -1.d0/3.d0
!      xer(5,2) = -1.d0
!      xer(6,1) =  1.d0/3.d0
!      xer(6,2) = -1.d0
!      xer(7,1) =  1.d0
!      xer(7,2) = -1.d0/3.d0
!      xer(8,1) =  1.d0
!      xer(8,2) =  1.d0/3.d0
!      xer(9,1) =  1.d0/3.d0
!      xer(9,2) =  1.d0
!      xer(10,1) = -1.d0/3.d0
!      xer(10,2) =  1.d0
!      xer(11,1) = -1.d0
!      xer(11,2) =  1.d0/3.d0
!      xer(12,1) = -1.d0
!      xer(12,2) = -1.d0/3.d0
!      xer(13,1) = -1.d0/3.d0
!      xer(13,2) = -1.d0/3.d0
!      xer(14,1) =  1.d0/3.d0
!      xer(14,2) = -1.d0/3.d0
!      xer(15,1) =  1.d0/3.d0
!      xer(15,2) =  1.d0/3.d0
!      xer(16,1) = -1.d0/3.d0
!      xer(16,2) =  1.d0/3.d0
!!
!      case(8)
!!
!!     serendipity 8
!!
!      xer(1,1) = -1.d0
!      xer(1,2) = -1.d0
!      xer(2,1) =  1.d0
!      xer(2,2) = -1.d0
!      xer(3,1) =  1.d0
!      xer(3,2) =  1.d0
!      xer(4,1) = -1.d0
!      xer(4,2) =  1.d0
!!
!      xer(5,1) =  0.d0
!      xer(5,2) = -1.d0
!      xer(6,1) =  1.d0
!      xer(6,2) =  0.d0
!      xer(7,1) =  0.d0
!      xer(7,2) =  1.d0
!      xer(8,1) = -1.d0
!      xer(8,2) =  0.d0
!!
!      case default
!!
!      write(*,*) 'element degree not defined'
!      stop
!!
!      end select
!
!
!!
!      p   =0.d0
!!
!!     loop on elements
!!
!      do i=1,nelem
!!
!!     coordinates of the element nodes
!!
!      do j=1,nen
!      do k=1,nsd
!      xe(j,k)=coord(el(i)%n(j),k)
!      end do ! k
!      end do ! j
!!
!!     solution in the nodes of the element
!!
!      do j=1,nen
!!
!      do k=1,nsd
!      xx(k)=xer(j,k)
!      end do ! k
!!
!!     construct of solution in each element:
!!
!      do k = 1,nen
!      p(i,j) = p(i,j) + phi2d(k,xx)*bb(el(i)%n(k),1)
!!
!      end do ! k
!      end do ! j
!!
!      end do ! i
!!
!!      do i = 1,nelem
!!         do j = 1,nelem
!!            p(i,j) = ph((i-1)*nelem + j)
!!         end do
!!      end do
!!
!!
!      end subroutine postprocessing
!!
!----------------------------------------------------------------------
!
      subroutine norms
!
      use mgeometry
      use mgauss
      use marrays
      use merro
!
      implicit none
!
      real(8), dimension(nsd)         :: xil,xx
      real(8), dimension(nen,nsd)     :: xe
      real(8), dimension(nsd  ,nsd  ) :: jf, invjf, tinvjf
      ! real(8) :: phi2d, dphi2dx, dphi2dy
      real(8) :: detjf
      integer :: i,j,k,l,ni,nj
      real(8), dimension(2) :: grad
!
      real(8) :: elsh,elih,elrh,elsv,eliv
      real(8) :: el2sh,el2ih,el2rh,el2sv,el2iv
      real(8) :: she,ihe,rhe,sve,ive
      real(8) :: shh,ihh,rhh,svh,ivh
      real(8) :: el2dsh,el2dih,el2drh,el2dsv,el2div
      real(8), dimension(2) :: dshe,dihe,drhe,dsve,dive
      real(8), dimension(2) :: dshh,dihh,drhh,dsvh,divh
      real(8), dimension(2) :: eldsh,eldih,eldrh,eldsv,eldiv
!
!      integer :: nint
!
      ! nint=(5)**nsd
!
!     integrating
!
      el2sh = 0.d0
      el2ih = 0.d0
      el2rh = 0.d0
      el2sv = 0.d0
      el2iv = 0.d0
      el2dsh = 0.d0
      el2dih = 0.d0
      el2drh = 0.d0
      el2dsv = 0.d0
      el2div = 0.d0
!
!     loop on elements
!
      do i=1,nelem
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(i)%n(j),k)
      end do ! k
      end do ! j
!
!     loop on integration points
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     Jacobian of the transformation
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1)= jf(2,2)
      invjf(1,2)=-jf(1,2)
      invjf(2,1)=-jf(2,1)
      invjf(2,2)= jf(1,1)
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
!     real coordinates of the integration point
!
      xx=0.d0
      do ni=1,nen
      do nj=1,nsd
      xx(nj)=xx(nj)+phi2d(ni,xil)*xe(ni,nj)
      end do ! nj
      end do ! ni
!
!     exact solutions at the integration points
!
      call exacts(she,ihe,rhe,sve,ive,xx,tt+dt)
      call exactd(dshe,dihe,drhe,dsve,dive,xx,tt+dt)
!
!     approximated solutions at the integration points
!
!     construct of solution in each element:
!
      shh = 0.d0
      ihh = 0.d0
      rhh = 0.d0
      svh = 0.d0
      ivh = 0.d0
      dshh = 0.d0
      dihh = 0.d0
      drhh = 0.d0
      dsvh = 0.d0
      divh = 0.d0
!
      do k = 1,nen
!
      shh = shh + phi2d(k,xil)*suh(el(i)%n(k),nt+1)
      ihh = ihh + phi2d(k,xil)*inh(el(i)%n(k),nt+1)
      rhh = rhh + phi2d(k,xil)*reh(el(i)%n(k),nt+1)
      svh = svh + phi2d(k,xil)*suv(el(i)%n(k),nt+1)
      ivh = ivh + phi2d(k,xil)*inv(el(i)%n(k),nt+1)
!
      grad(1) = tinvjf(1,1)*dphi2dx(k,xil)+tinvjf(1,2)*dphi2dy(k,xil)
      grad(2) = tinvjf(2,1)*dphi2dx(k,xil)+tinvjf(2,2)*dphi2dy(k,xil)
!
      dshh(1) = dshh(1) + grad(1)*suh(el(i)%n(k),nt+1)
      dihh(1) = dihh(1) + grad(1)*inh(el(i)%n(k),nt+1)
      drhh(1) = drhh(1) + grad(1)*reh(el(i)%n(k),nt+1)
      dsvh(1) = dsvh(1) + grad(1)*suv(el(i)%n(k),nt+1)
      divh(1) = divh(1) + grad(1)*inv(el(i)%n(k),nt+1)
!
      dshh(2) = dshh(2) + grad(2)*suh(el(i)%n(k),nt+1)
      dihh(2) = dihh(2) + grad(2)*inh(el(i)%n(k),nt+1)
      drhh(2) = drhh(2) + grad(2)*reh(el(i)%n(k),nt+1)
      dsvh(2) = dsvh(2) + grad(2)*suv(el(i)%n(k),nt+1)
      divh(2) = divh(2) + grad(2)*inv(el(i)%n(k),nt+1)
!
      end do ! k
!
      elsh  = dabs(shh-she)
      elih  = dabs(ihh-ihe)
      elrh  = dabs(rhh-rhe)
      elsv  = dabs(svh-sve)
      eliv  = dabs(ivh-ive)

      eldsh = dabs(dshh-dshe)
      eldih = dabs(dihh-dihe)
      eldrh = dabs(drhh-drhe)
      eldsv = dabs(dsvh-dsve)
      eldiv = dabs(divh-dive)
!
!     l2 norm
!
      el2sh  = el2sh + (elsh**2)*wg(l,nint)*detjf
      el2ih  = el2ih + (elih**2)*wg(l,nint)*detjf
      el2rh  = el2rh + (elrh**2)*wg(l,nint)*detjf
      el2sv  = el2sv + (elsv**2)*wg(l,nint)*detjf
      el2iv  = el2iv + (eliv**2)*wg(l,nint)*detjf
!
!     grad norm
!
      el2dsh = el2dsh + (eldsh(1)**2 + eldsh(2)**2)*wg(l,nint)*detjf
      el2dih = el2dih + (eldih(1)**2 + eldih(2)**2)*wg(l,nint)*detjf
      el2drh = el2drh + (eldrh(1)**2 + eldrh(2)**2)*wg(l,nint)*detjf
      el2dsv = el2dsv + (eldsv(1)**2 + eldsv(2)**2)*wg(l,nint)*detjf
      el2div = el2div + (eldiv(1)**2 + eldiv(2)**2)*wg(l,nint)*detjf
!
      end do ! l
!
      end do ! i
!
      el2sh  = dsqrt(el2sh)
      el2ih  = dsqrt(el2ih)
      el2rh  = dsqrt(el2rh)
      el2sv  = dsqrt(el2sv)
      el2iv  = dsqrt(el2iv)
!
      el2dsh = dsqrt(el2dsh)
      el2dih = dsqrt(el2dih)
      el2drh = dsqrt(el2drh)
      el2dsv = dsqrt(el2dsv)
      el2div = dsqrt(el2div)
!
      open(unit = 111,file = 'norms.dat',status = 'unknown')
!
      errol2(expnum,1) = el2sh
      errol2(expnum,2) = el2ih
      errol2(expnum,3) = el2rh
      errol2(expnum,4) = el2sv
      errol2(expnum,5) = el2iv
      errol2d(expnum,1) = el2dsh
      errol2d(expnum,2) = el2dih
      errol2d(expnum,3) = el2drh
      errol2d(expnum,4) = el2dsv
      errol2d(expnum,5) = el2div
!
      write(111,"(15(e15.7))") errol2(expnum,1),errol2(expnum,2),&
         errol2(expnum,3),errol2(expnum,4),errol2(expnum,5), &
         errol2d(expnum,1),errol2d(expnum,2),&
         errol2d(expnum,3),errol2d(expnum,4),errol2d(expnum,5)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine normsad
!
      use mgeometry
      use mgauss
      use marrays
      use merro
!
      implicit none
!
      real(8), dimension(nsd)         :: xil,xx
      real(8), dimension(nen,nsd)     :: xe
      real(8), dimension(nsd  ,nsd  ) :: jf, invjf, tinvjf
      ! real(8) :: phi2d, dphi2dx, dphi2dy
      real(8) :: detjf
      integer :: i,j,k,l,ni,nj
      real(8), dimension(2) :: grad
!
      real(8) :: elsh,elih,elrh,elsv,eliv
      real(8) :: el2sh,el2ih,el2rh,el2sv,el2iv
      real(8) :: she,ihe,rhe,sve,ive
      real(8) :: shh,ihh,rhh,svh,ivh
      real(8) :: el2dsh,el2dih,el2drh,el2dsv,el2div
      real(8), dimension(2) :: dshe,dihe,drhe,dsve,dive
      real(8), dimension(2) :: dshh,dihh,drhh,dsvh,divh
      real(8), dimension(2) :: eldsh,eldih,eldrh,eldsv,eldiv
!
!      integer :: nint
!
      ! nint=(5)**nsd
!
!     integrating
!
      el2sh = 0.d0
      el2ih = 0.d0
      el2rh = 0.d0
      el2sv = 0.d0
      el2iv = 0.d0
      el2dsh = 0.d0
      el2dih = 0.d0
      el2drh = 0.d0
      el2dsv = 0.d0
      el2div = 0.d0
!
!     loop on elements
!
      do i=1,nelem
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(i)%n(j),k)
      end do ! k
      end do ! j
!
!     loop on integration points
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     Jacobian of the transformation
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1)= jf(2,2)
      invjf(1,2)=-jf(1,2)
      invjf(2,1)=-jf(2,1)
      invjf(2,2)= jf(1,1)
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
!     real coordinates of the integration point
!
      xx=0.d0
      do ni=1,nen
      do nj=1,nsd
      xx(nj)=xx(nj)+phi2d(ni,xil)*xe(ni,nj)
      end do ! nj
      end do ! ni
!
!     exact solutions at the integration points
!
      call exacts(she,ihe,rhe,sve,ive,xx,tt-dt)
      call exactd(dshe,dihe,drhe,dsve,dive,xx,tt-dt)
!
!     approximated solutions at the integration points
!
!     construct of solution in each element:
!
      shh = 0.d0
      ihh = 0.d0
      rhh = 0.d0
      svh = 0.d0
      ivh = 0.d0
      dshh = 0.d0
      dihh = 0.d0
      drhh = 0.d0
      dsvh = 0.d0
      divh = 0.d0
!
      do k = 1,nen
!
      shh = shh + phi2d(k,xil)*l1(el(i)%n(k),1)
      ihh = ihh + phi2d(k,xil)*l2(el(i)%n(k),1)
      rhh = rhh + phi2d(k,xil)*l3(el(i)%n(k),1)
      svh = svh + phi2d(k,xil)*l4(el(i)%n(k),1)
      ivh = ivh + phi2d(k,xil)*l5(el(i)%n(k),1)
!
      grad(1) = tinvjf(1,1)*dphi2dx(k,xil)+tinvjf(1,2)*dphi2dy(k,xil)
      grad(2) = tinvjf(2,1)*dphi2dx(k,xil)+tinvjf(2,2)*dphi2dy(k,xil)
!
      dshh(1) = dshh(1) + grad(1)*l1(el(i)%n(k),1)
      dihh(1) = dihh(1) + grad(1)*l2(el(i)%n(k),1)
      drhh(1) = drhh(1) + grad(1)*l3(el(i)%n(k),1)
      dsvh(1) = dsvh(1) + grad(1)*l4(el(i)%n(k),1)
      divh(1) = divh(1) + grad(1)*l5(el(i)%n(k),1)
!
      dshh(2) = dshh(2) + grad(2)*l1(el(i)%n(k),1)
      dihh(2) = dihh(2) + grad(2)*l2(el(i)%n(k),1)
      drhh(2) = drhh(2) + grad(2)*l3(el(i)%n(k),1)
      dsvh(2) = dsvh(2) + grad(2)*l4(el(i)%n(k),1)
      divh(2) = divh(2) + grad(2)*l5(el(i)%n(k),1)
!
      end do ! k
!
      elsh  = dabs(shh-she)
      elih  = dabs(ihh-ihe)
      elrh  = dabs(rhh-rhe)
      elsv  = dabs(svh-sve)
      eliv  = dabs(ivh-ive)

      eldsh = dabs(dshh-dshe)
      eldih = dabs(dihh-dihe)
      eldrh = dabs(drhh-drhe)
      eldsv = dabs(dsvh-dsve)
      eldiv = dabs(divh-dive)
!
!     l2 norm
!
      el2sh  = el2sh + (elsh**2)*wg(l,nint)*detjf
      el2ih  = el2ih + (elih**2)*wg(l,nint)*detjf
      el2rh  = el2rh + (elrh**2)*wg(l,nint)*detjf
      el2sv  = el2sv + (elsv**2)*wg(l,nint)*detjf
      el2iv  = el2iv + (eliv**2)*wg(l,nint)*detjf
!
!     grad norm
!
      el2dsh = el2dsh + (eldsh(1)**2 + eldsh(2)**2)*wg(l,nint)*detjf
      el2dih = el2dih + (eldih(1)**2 + eldih(2)**2)*wg(l,nint)*detjf
      el2drh = el2drh + (eldrh(1)**2 + eldrh(2)**2)*wg(l,nint)*detjf
      el2dsv = el2dsv + (eldsv(1)**2 + eldsv(2)**2)*wg(l,nint)*detjf
      el2div = el2div + (eldiv(1)**2 + eldiv(2)**2)*wg(l,nint)*detjf
!
      end do ! l
!
      end do ! i
!
      el2sh  = dsqrt(el2sh)
      el2ih  = dsqrt(el2ih)
      el2rh  = dsqrt(el2rh)
      el2sv  = dsqrt(el2sv)
      el2iv  = dsqrt(el2iv)
!
      el2dsh = dsqrt(el2dsh)
      el2dih = dsqrt(el2dih)
      el2drh = dsqrt(el2drh)
      el2dsv = dsqrt(el2dsv)
      el2div = dsqrt(el2div)
!
      open(unit = 111,file = 'norms.dat',status = 'unknown')
!
      errol2(expnum,1) = el2sh
      errol2(expnum,2) = el2ih
      errol2(expnum,3) = el2rh
      errol2(expnum,4) = el2sv
      errol2(expnum,5) = el2iv
      errol2d(expnum,1) = el2dsh
      errol2d(expnum,2) = el2dih
      errol2d(expnum,3) = el2drh
      errol2d(expnum,4) = el2dsv
      errol2d(expnum,5) = el2div
!
      write(111,"(15(e15.7))") errol2(expnum,1),errol2(expnum,2),&
         errol2(expnum,3),errol2(expnum,4),errol2(expnum,5), &
         errol2d(expnum,1),errol2d(expnum,2),&
         errol2d(expnum,3),errol2d(expnum,4),errol2d(expnum,5)
!
      end subroutine
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      subroutine exacts(she,ihe,rhe,sve,ive,xx,tt)
!
      use mgeometry, only : nsd
!
      implicit none
!
      real(8) :: she,ihe,rhe,sve,ive
      real(8) :: tt, pi, aux
      real(8), dimension(nsd) :: xx
!
      pi = 4.d0*datan(1.d0)
      ! aux = dcos(pi*xx(1))*dcos(pi*xx(2))*dcos(2.d0*pi*tt)!dexp(-tt/10.d0)
      ! aux = dcos(pi*xx(1))*dcos(pi*xx(2))+dcos(2.d0*pi*tt)!dexp(-tt/10.d0)
      ! aux = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(xx(2)**3/3.d0-xx(2)**2/2.d0)&
      !    *dexp(-tt/10.d0)
      ! aux = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(xx(2)**3/3.d0-xx(2)**2/2.d0)&
      !    *dsin(pi*tt)
      aux = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(xx(2)**3/3.d0-xx(2)**2/2.d0)&
         *dcos(pi*tt)
!
      she = aux + 10.d0
      ihe = aux + 5.d0
      rhe = aux + 3.d0
      sve = aux + 2.d0
      ive = aux + 1.d0
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine exactd(dshe,dihe,drhe,dsve,dive,xx,tt)
!
      use mgeometry, only : nsd
!
      implicit none
!
      real(8), dimension(2) :: dshe,dihe,drhe,dsve,dive
      real(8) :: tt, pi, auxx, auxy
      real(8), dimension(nsd) :: xx
!
      pi = 4.d0*datan(1.d0)
      ! auxx = -pi*dsin(pi*xx(1))*dcos(pi*xx(2))*dcos(2.d0*pi*tt)!dexp(-tt/10.d0)
      ! auxy = -pi*dcos(pi*xx(1))*dsin(pi*xx(2))*dcos(2.d0*pi*tt)!dexp(-tt/10.d0)
      ! auxx = -pi*dsin(pi*xx(1))*dcos(pi*xx(2))
      ! auxy = -pi*dcos(pi*xx(1))*dsin(pi*xx(2))
      ! auxx = (xx(1)**2-xx(1))*(xx(2)**3/3.d0-xx(2)**2/2.d0)*dexp(-tt/10.d0)
      ! auxy = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(xx(2)**2-xx(2))*dexp(-tt/10.d0)
      ! auxx = (xx(1)**2-xx(1))*(xx(2)**3/3.d0-xx(2)**2/2.d0)*dsin(pi*tt)
      ! auxy = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(xx(2)**2-xx(2))*dsin(pi*tt)
      auxx = (xx(1)**2-xx(1))*(xx(2)**3/3.d0-xx(2)**2/2.d0)*dcos(pi*tt)
      auxy = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(xx(2)**2-xx(2))*dcos(pi*tt)
!
      dshe(1) = auxx
      dihe(1) = auxx
      drhe(1) = auxx
      dsve(1) = auxx
      dive(1) = auxx
   !
      dshe(2) = auxy
      dihe(2) = auxy
      drhe(2) = auxy
      dsve(2) = auxy
      dive(2) = auxy
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine exactd2(d2she,d2ihe,d2rhe,d2sve,d2ive,xx,tt)
!
      use mgeometry, only : nsd
!
      implicit none
!
      real(8), dimension(2) :: d2she,d2ihe,d2rhe,d2sve,d2ive
      real(8) :: tt, pi, pi2, aux, auxx, auxy
      real(8), dimension(nsd) :: xx
!
      pi = 4.d0*datan(1.d0)
      pi2 = pi*pi
      ! aux = -pi2*dcos(pi*xx(1))*dcos(pi*xx(2))*dcos(2.d0*pi*tt)!dexp(-tt/10.d0)
      ! aux = -pi2*dcos(pi*xx(1))*dcos(pi*xx(2))
      ! auxx = (2.d0*xx(1)-1.d0)*(xx(2)**3/3.d0-xx(2)**2/2.d0)*dexp(-tt/10.d0)
      ! auxy = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(2.d0*xx(2)-1.d0)*dexp(-tt/10.d0)
      ! auxx = (2.d0*xx(1)-1.d0)*(xx(2)**3/3.d0-xx(2)**2/2.d0)*dsin(pi*tt)
      ! auxy = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(2.d0*xx(2)-1.d0)*dsin(pi*tt)
      auxx = (2.d0*xx(1)-1.d0)*(xx(2)**3/3.d0-xx(2)**2/2.d0)*dcos(pi*tt)
      auxy = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(2.d0*xx(2)-1.d0)*dcos(pi*tt)
!
      d2she(1) = auxx
      d2ihe(1) = auxx
      d2rhe(1) = auxx
      d2sve(1) = auxx
      d2ive(1) = auxx
!
      d2she(2) = auxy
      d2ihe(2) = auxy
      d2rhe(2) = auxy
      d2sve(2) = auxy
      d2ive(2) = auxy
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine exactdt(dtshe,dtihe,dtrhe,dtsve,dtive,xx,tt)
!
      use mgeometry, only : nsd
!
      implicit none
!
      real(8) :: dtshe,dtihe,dtrhe,dtsve,dtive
      real(8) :: tt, pi, aux
      real(8), dimension(nsd) :: xx
!
      pi = 4.d0*datan(1.d0)
      ! aux = -dcos(pi*xx(1))*dcos(pi*xx(2))*dexp(-tt/10.d0)/10.d0
      ! aux = -2.d0*pi*dcos(pi*xx(1))*dcos(pi*xx(2))*dsin(2.d0*pi*tt)
      ! aux = -2.d0*pi*dsin(2.d0*pi*tt)
      ! aux = (xx(1)**3/3.d0-xx(1)**2/2.d0)*(xx(2)**3/3.d0-xx(2)**2/2.d0)&
      !    *dexp(-tt/10.d0)*(-1/10.d0)
      aux = -(xx(1)**3/3.d0-xx(1)**2/2.d0)*(xx(2)**3/3.d0-xx(2)**2/2.d0)&
         *dsin(pi*tt)*pi
!
      dtshe = aux
      dtihe = aux
      dtrhe = aux
      dtsve = aux
      dtive = aux
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine source(fsource,xx,tt)
!
      use mgeometry, only : nsd
      use mcoeficientes
!
      implicit none
!
      real(8), dimension(5) :: fsource
      real(8) :: tt
      real(8), dimension(nsd) :: xx
      real(8) :: she,ihe,rhe,sve,ive,u
      real(8) :: dtshe,dtihe,dtrhe,dtsve,dtive
      real(8), dimension(2) :: dshe,dihe,drhe,dsve,dive
      real(8), dimension(2) :: d2she,d2ihe,d2rhe,d2sve,d2ive
!
      call exacts(she,ihe,rhe,sve,ive,xx,tt)
      call exactd2(d2she,d2ihe,d2rhe,d2sve,d2ive,xx,tt)
      call exactdt(dtshe,dtihe,dtrhe,dtsve,dtive,xx,tt)
!
      u = min(umax,max((she-rhe-c1)*she/(2.d0*c2),0.d0))
!
      fsource(1) = dtshe - alpha(1)*(d2she(1)+d2she(2)) &
         +beta*she*ive-lamb*(she+ihe+rhe)*(1.d0-(she+ihe+rhe)/kp*ikp) &
         +mu*she-sig*ihe+u*she
      fsource(2) = dtihe - alpha(2)*(d2ihe(1)+d2ihe(2)) &
         -beta*she*ive+(sig+delta+mui)*ihe
      fsource(3) = dtrhe - alpha(3)*(d2rhe(1)+d2rhe(2)) &
         +mu*rhe-u*she-delta*ihe
      fsource(4) = dtsve - alpha(4)*(d2sve(1)+d2sve(2)) &
         +betav*sve*ihe-lambv*(sve+ive)*(1.d0-(sve+ive)/kpv) !&
         !+muv*sve
      fsource(5) = dtive - alpha(5)*(d2ive(1)+d2ive(2)) &
         -betav*sve*ihe+muv*ive
!
      end subroutine
!
!----------------------------------------------------------------------
!
!
      subroutine sourcead(fsource,xx,tt)
!
      use mgeometry, only : nsd
      use mcoeficientes
!
      implicit none
!
      real(8), dimension(5) :: fsource
      real(8) :: tt
      real(8), dimension(nsd) :: xx
      real(8) :: l1e,l2e,l3e,l4e,l5e
      real(8) :: dtl1e,dtl2e,dtl3e,dtl4e,dtl5e
      real(8), dimension(2) :: d2l1e,d2l2e,d2l3e,d2l4e,d2l5e
      real(8) :: shj,ihj,rhj,svj,ivj,u
!
      call exacts(shj,ihj,rhj,svj,ivj,xx,tt)
!
      call exacts(l1e,l2e,l3e,l4e,l5e,xx,tt)
      call exactd2(d2l1e,d2l2e,d2l3e,d2l4e,d2l5e,xx,tt)
      call exactdt(dtl1e,dtl2e,dtl3e,dtl4e,dtl5e,xx,tt)
!
      u = min(umax,max((l1e-l3e-c1)*shj/(2.d0*c2),0.d0))
!
      fsource(1) = dtl1e + alpha(1)*(d2l1e(1)+d2l1e(2)) &
         - (beta*ivj-lamb+mu+2.d0*lamb*(shj+ihj+rhj)/kp*ikp+u)*l1e &
         + beta*ivj*l2e + u*l3e + c1*u
      fsource(2) = dtl2e + alpha(2)*(d2l2e(1)+d2l2e(2)) &
         - (-lamb+2.d0*lamb*(shj+ihj+rhj)/kp*ikp-sig)*l1e &
         -(sig+mui+delta)*l2e + delta*l3e - betav*svj*l4e + betav*svj*l5e + 1.d0
      fsource(3) = dtl3e + alpha(3)*(d2l3e(1)+d2l3e(2)) &
         - (-lamb+2.d0*lamb*(shj+ihj+rhj)/kp*ikp)*l1e - mu*l3e
      fsource(4) = dtl4e + alpha(4)*(d2l4e(1)+d2l4e(2)) &
         ! - (betav*ihj-lambv+muv+2.d0*lambv*(svj+ivj)/kpv)*l4e &
         - (betav*ihj-lambv+2.d0*lambv*(svj+ivj)/kpv)*l4e &
         +betav*ihj*l5e
      fsource(5) = dtl5e + alpha(5)*(d2l5e(1)+d2l5e(2)) &
         - beta*shj*l1e+beta*shj*l2e &
         -(-lambv+2.d0*lambv*(svj+ivj)/kpv)*l4e-muv*l5e
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine initial
!
      use mgeometry
      use marrays
      use mcoeficientes
!
      implicit none
!
      integer :: i,j
      real(8) :: she,ihe,rhe,sve,ive
      real(8), dimension(2) :: x
      ! integer, dimension(23,1) :: elem
      integer, dimension(13,1) :: elem
!
      suh = 0.d0
      inh = 0.d0
      reh = 0.d0
      suv = 0.d0
      inv = 0.d0
! 
      suh(:,1) = s0
      reh(:,1) = r0
      suv(:,1) = sv0
      inv(:,1) = iv0
! 
      ! open(unit = 400,file = 'elemcities.dat',status = 'old')
      ! elem = 0
      ! do i = 1,14
      !    read(400,*) (elem(i,j),j=1,7)
      ! end do
      ! close(unit=400)
      ! elem(:,1) = (/299,353,426,29,70,156,69,48,69,247,60,169,406,386/) !6 weeks
     elem(:,1) = (/132,74,247,66,299,426,156,20,406,24,386,316,69/)
 ! elem(:,1) = 1
! 
      ! do i=1,6
      do j=1,3
! 
      inh(el(elem(5,1))%n(j),1) = i01!33.d0/4.3646d-3/10195921.63d0/0.2d0
      ! inv(el(elem(5,1))%n(j),1) = 5.d0*inh(el(elem(5,1))%n(j),1)
! 
!       inh(el(elem(3,1))%n(j),1) = i02!33.d0/4.3646d-3/10195921.63d0/0.2d0
!       inv(el(elem(3,1))%n(j),1) = 5.d0*inh(el(elem(3,1))%n(j),1)
! ! 
!       inh(el(elem(7,1))%n(j),1) = i03!33.d0/4.3646d-3/10195921.63d0/0.2d0
!       inv(el(elem(7,1))%n(j),1) = 5.d0*inh(el(elem(7,1))%n(j),1)
! 
      end do !j
      ! end do !i
! 
      suh(:,1) = suh(:,1) - inh(:,1)
      suv(:,1) = suv(:,1) - inv(:,1)
!
!     testing matlab comparison
     ! inh = 0.d0
      ! elem(:,1) = (/77,84,103,131,132,137,142,176,186,189,193,194,197,202,222,230,245,246,247,248,255,260,264/)

      ! do i=1,23
      ! inh(elem(i,1),1) = i01
      ! end do !i
!
!......................................
!       ! open(unit=21,file='solutions/noise0/sol001suh.dat',status='old')
!       ! open(unit=22,file='solutions/noise0/sol001inh.dat',status='old')
!       ! open(unit=23,file='solutions/noise0/sol001reh.dat',status='old')
!       ! open(unit=24,file='solutions/noise0/sol001suv.dat',status='old')
!       ! open(unit=25,file='solutions/noise0/sol001inv.dat',status='old')
!       open(unit=21,file='solutions/noise0t70/sol001suh.dat',status='old')
!       open(unit=22,file='solutions/noise0t70/sol001inh.dat',status='old')
!       open(unit=23,file='solutions/noise0t70/sol001reh.dat',status='old')
!       open(unit=24,file='solutions/noise0t70/sol001suv.dat',status='old')
!       open(unit=25,file='solutions/noise0t70/sol001inv.dat',status='old')
! !
!       do i = 1,ndof
!          read(21,*) (suh(i,j),j=1,(nt+1))
!          read(22,*) (inh(i,j),j=1,(nt+1))
!          read(23,*) (reh(i,j),j=1,(nt+1))
!          read(24,*) (suv(i,j),j=1,(nt+1))
!          read(25,*) (inv(i,j),j=1,(nt+1))
!       end do !i
! !
!       close(unit=21)
!       close(unit=22)
!       close(unit=23)
!       close(unit=24)
!       close(unit=25)
! !
!       suh(:,1) = suh(:,nt+1)
!       inh(:,1) = inh(:,nt+1)
!       reh(:,1) = reh(:,nt+1)
!       suv(:,1) = suv(:,nt+1)
!       inv(:,1) = inv(:,nt+1)
! ! 
!       suh(:,2:nt+1) = 0.d0
!       inh(:,2:nt+1) = 0.d0
!       reh(:,2:nt+1) = 0.d0
!       suv(:,2:nt+1) = 0.d0
!       inv(:,2:nt+1) = 0.d0
!......................................
!
      suh(:,2) = suh(:,1)
      inh(:,2) = inh(:,1)
      reh(:,2) = reh(:,1)
      suv(:,2) = suv(:,1)
      inv(:,2) = inv(:,1)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine initialad
!
      use mgeometry
      use marrays
      use mcoeficientes
!
      implicit none
!
      integer :: i,j
      real(8) :: she,ihe,rhe,sve,ive
      real(8), dimension(2) :: x
!
      l1 = 0.d0
      l2 = 0.d0
      l3 = 0.d0
      l4 = 0.d0
      l5 = 0.d0
!
!       do i=1,nelem
!       do j=1,nen
!       x(1) = coord(el(i)%n(j),1)
!       x(2) = coord(el(i)%n(j),2)
! !
!       call exacts(she,ihe,rhe,sve,ive,x,0.d0)
! !
!       l1(el(i)%n(j),nt+1) = she
!       l2(el(i)%n(j),nt+1) = ihe
!       l3(el(i)%n(j),nt+1) = rhe
!       l4(el(i)%n(j),nt+1) = sve
!       l5(el(i)%n(j),nt+1) = ive
! !
!       end do !j
!       end do !i
!
      l1(:,nt) = l1(:,nt+1)
      l2(:,nt) = l2(:,nt+1)
      l3(:,nt) = l3(:,nt+1)
      l4(:,nt) = l4(:,nt+1)
      l5(:,nt) = l5(:,nt+1)
!
      ! l1 = 0.d0
      ! l2 = 0.d0
      ! l3 = 0.d0
      ! l4 = 0.d0
      ! l5 = 0.d0
      ! l1st = 0.d0
      ! l2st = 0.d0
      ! l3st = 0.d0
      ! l4st = 0.d0
      ! l5st = 0.d0
      ! l1t = 0.d0
      ! l2t = 0.d0
      ! l3t = 0.d0
      ! l4t = 0.d0
      ! l5t = 0.d0
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine plotsols
!
      use mgeometry
      use marrays
!
      implicit none
!
      integer :: i, j
      character(len = 128) :: namesh,nameih,namerh,namesv,nameiv,namect,exp
!
      write(exp,'(i3.3)') expnum
!
      namesh = 'solutions/sol'//trim(exp)//trim('suh')//'.dat'
      nameih = 'solutions/sol'//trim(exp)//trim('inh')//'.dat'
      namerh = 'solutions/sol'//trim(exp)//trim('reh')//'.dat'
      namesv = 'solutions/sol'//trim(exp)//trim('suv')//'.dat'
      nameiv = 'solutions/sol'//trim(exp)//trim('inv')//'.dat'
      namect = 'solutions/sol'//trim(exp)//trim('ctl')//'.dat'
!
      open(unit = 201, file = namesh, status = 'unknown')
      open(unit = 202, file = nameih, status = 'unknown')
      open(unit = 203, file = namerh, status = 'unknown')
      open(unit = 204, file = namesv, status = 'unknown')
      open(unit = 205, file = nameiv, status = 'unknown')
      open(unit = 206, file = namect, status = 'unknown')
!
      do i=1,ndof
!
      write(201,*) (suh(i,j), j=1,(nt+1))
      write(202,*) (inh(i,j), j=1,(nt+1))
      write(203,*) (reh(i,j), j=1,(nt+1))
      write(204,*) (suv(i,j), j=1,(nt+1))
      write(205,*) (inv(i,j), j=1,(nt+1))
      write(206,*) (u(i,j), j=1,(nt+1))
!
      end do !i
!
      close(unit = 201)
      close(unit = 202)
      close(unit = 203)
      close(unit = 204)
      close(unit = 205)
      close(unit = 206)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine plotints
!
      use mgeometry
      use marrays
!
      implicit none
!
      character(len = 128) :: namesh,nameih,namerh,namesv,nameiv,namect,namej,exp
!
      write(exp,'(i3.3)') expnum
!
      namesh = 'solutions/int'//trim(exp)//trim('suh')//'.dat'
      nameih = 'solutions/int'//trim(exp)//trim('inh')//'.dat'
      namerh = 'solutions/int'//trim(exp)//trim('reh')//'.dat'
      namesv = 'solutions/int'//trim(exp)//trim('suv')//'.dat'
      nameiv = 'solutions/int'//trim(exp)//trim('inv')//'.dat'
      namect = 'solutions/int'//trim(exp)//trim('ctl')//'.dat'
      namej = 'solutions/int'//trim(exp)//trim('j')//'.dat'
!
      open(unit = 101, file = namesh, status = 'unknown')
      open(unit = 102, file = nameih, status = 'unknown')
      open(unit = 103, file = namerh, status = 'unknown')
      open(unit = 104, file = namesv, status = 'unknown')
      open(unit = 105, file = nameiv, status = 'unknown')
      open(unit = 106, file = namect, status = 'unknown')
      open(unit = 107, file = namej, status = 'unknown')
!
      write(101,*) intsuh
      write(102,*) intinh
      write(103,*) intreh
      write(104,*) intsuv
      write(105,*) intinv
      write(106,*) intctl
      write(107,*) intj
!
      close(unit = 101)
      close(unit = 102)
      close(unit = 103)
      close(unit = 104)
      close(unit = 105)
      close(unit = 106)
      close(unit = 107)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine plotintslhs
!
      use mgeometry
      use marrays
!
      implicit none
!
      character(len = 128) :: namej,namen,exp
!
      write(exp,'(i3.3)') expnum
!
      namej = 'solutions/lhs/int'//trim(exp)//trim('j')//'.dat'
      namen = 'solutions/lhs/int'//trim(exp)//trim('nnc')//'.dat'
!
      open(unit = 107, file = namej, status = 'unknown')
      open(unit = 108, file = namen, status = 'unknown')
!
      write(107,*) intj*10195921.6313475d0*66.67d0
      write(108,*) intnnc*10195921.6313475d0
!
      close(unit = 107)
      close(unit = 108)
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine grafsurf(u,pop)
!----------------------------------------------------------------------
!     chamada dos modulos
!--------------------------------------------------------------------------
      use mgeometry
!--------------------------------------------------------------------------
!     declaracao das variaveis
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: i, aux
      real(8), dimension(nnodes,1) :: u
!      real(8) :: aux
      character(len = 128) :: filename, itaux, expnumaux
      character(len = 3) :: pop
!
!--------------------------------------------------------------------------
!     salvando em arquivo
!
      aux = (100*it)/nt
!      aux = it/8
      write (itaux,'(i3.3)') aux
!      aux = it/8
!      write (itaux,'(i3.3)') aux
!
      write(expnumaux,'(i3.3)') expnum
!
      filename = 'solutions/sol'//trim(expnumaux)//trim(pop)//trim(itaux)//'.dat'
!
      open(unit = 10, file = filename, status = 'unknown')
!
!     write(10,*), '#        u'
!     write(*,*), '#        u'
!
      do i = 1, nnodes
         write(10,200) u(i,1)
         !write(*,200) u(i)
200      format(e22.15)
      end do
!
      close(unit = 10)
!
      end subroutine
!----------------------------------------------------------------------
!
      subroutine graftemp(u,uaux,pop,intu)

!     chamada dos modulos
!--------------------------------------------------------------------------
      use mgeometry
!--------------------------------------------------------------------------
!     declaracao das variaveis
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: i, j
      real(8) :: intu
      real(8), dimension(nnodes,1) :: u, aux
      real(8), dimension(nt+1,10) :: uaux
      character(len = 128) :: filename, itaux, expnumaux
      character(len = 3) :: pop
!
!--------------------------------------------------------------------------
!     guarda os valores para os graficos temporais
!
      do i = 1, 9
         uaux(it+1,i) = u(ind(i),1)
      end do
      uaux(it+1,10) = intu
!      do i = 1, 9
!         uaux(it+1,i) = 0.0d0
!      end do
!      uaux(it+1,1) = u(2*(nx+1)+3,1)
!      uaux(it+1,2) = u(nnodes-2*(nx+1)-2,1)
!
!     salvando em arquivo
!
      if (it == nt) then
      write(expnumaux,'(i3.3)') expnum
      filename = 'solutions/sol'//trim(expnumaux)//'temp'//trim(pop)//'.dat'
!
      open(unit = 10, file = filename, status = 'unknown')
!
      do i = 1, nt+1
         write(10,200) (uaux(i,j),j=1,10)
200      format(10(e21.15,2x))
      end do
!
      close(unit = 10)

      end if
!--------------------------------------------------------------------------
      end subroutine
!--------------------------------------------------------------------------
!
      subroutine integral
!
      use mgeometry
      use mgauss
      use marrays
      use mcoeficientes
!
      implicit none
!
      real(8), dimension(nsd)         :: xil,xx
      real(8), dimension(nen,nsd)     :: xe
      real(8), dimension(nsd  ,nsd  ) :: jf, invjf, tinvjf
      ! real(8), dimension(nt+1,1) :: intjmax
!
      ! real(8) :: phi2d
      real(8) :: suhh,inhh,rehh,suvh,invh,ctlh,jh,jmaxh
      real(8) :: nnch,vach
      real(8) :: detjf
      integer :: i,j,k,l,ni,nj,nk
!
      real(8) :: jfunc, nnc, vac
!
      ! nint=(5)**nsd
!
      intsuh = 0.d0
      intinh = 0.d0
      intreh = 0.d0
      intsuv = 0.d0
      intinv = 0.d0
      intctl = 0.d0
      ! intjmax = 0.d0
      intj = 0.d0
      intnnc = 0.d0
      intvac = 0.d0
!
      do it=1,(nt+1)
!
!     integrating
!
!     loop on elements
!
      do i=1,nelem
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(i)%n(j),k)
      end do ! k
      end do ! j
!
!     loop on integration points
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     Jacobian of the transformation
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1)= jf(2,2)
      invjf(1,2)=-jf(1,2)
      invjf(2,1)=-jf(2,1)
      invjf(2,2)= jf(1,1)
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
!     real coordinates of the integration point
!
      xx=0.d0
      do ni=1,nen
      do nj=1,nsd
      xx(nj)=xx(nj)+phi2d(ni,xil)*xe(ni,nj)
      end do ! nj
      end do ! ni
!
!     approximated solutions at the integration points
!
!      ph=ph+phip(k,xil)*el(i)%p(k,1)!
!     construct of solution in each element:
!
      suhh = 0.d0
      inhh = 0.d0
      rehh = 0.d0
      suvh = 0.d0
      invh = 0.d0
      ctlh = 0.d0
      jh = 0.d0
      jmaxh= 0.d0
      nnch = 0.d0
      vach = 0.d0
! 
      do k = 1,nen
      suhh = suhh + phi2d(k,xil)*suh(el(i)%n(k),it)
      inhh = inhh + phi2d(k,xil)*inh(el(i)%n(k),it)
      rehh = rehh + phi2d(k,xil)*reh(el(i)%n(k),it)
      suvh = suvh + phi2d(k,xil)*suv(el(i)%n(k),it)
      invh = invh + phi2d(k,xil)*inv(el(i)%n(k),it)
      ctlh = ctlh + phi2d(k,xil)*u(el(i)%n(k),it)
      ! jh = jh + phi2d(k,xil)*inh(el(i)%n(k),it)*10195921.6313475d0*66.67d0
      jh = jh + c1*phi2d(k,xil)*inh(el(i)%n(k),it)
      ! jmaxh = jmaxh + c3*phi2d(k,xil)*(inh(el(i)%n(k),it)+ &
      !       (c1)*umax*suh(el(i)%n(k),it)+(c2)*umax*umax)
      do nk = 1,nen
      ! jh = jh + c3*phi2d(k,xil)*(phi2d(nk,xil)*&
      !    ((c1)*u(el(i)%n(nk),it)*suh(el(i)%n(k),it)*10195921.6313475d0*66.67d0 + &
      !    (c2)*u(el(i)%n(nk),it)*u(el(i)%n(k),it)*10195921.6313475d0*66.67d0/398.2514/398.2514))
      jh = jh + phi2d(k,xil)*(phi2d(nk,xil)*&
         (c2*u(el(i)%n(nk),it)*suh(el(i)%n(k),it) + &
         c3*u(el(i)%n(nk),it)*u(el(i)%n(k),it)))
      !
      nnch = nnch + (beta*suh(el(i)%n(nk),it)*inv(el(i)%n(k),it) +&
         betas*suh(el(i)%n(nk),it)*inh(el(i)%n(k),it))*&
         phi2d(k,xil)*phi2d(nk,xil)
      vach = vach + u(el(i)%n(nk),it)*suh(el(i)%n(k),it)*&
         phi2d(k,xil)*phi2d(nk,xil)
      !
      end do ! nk
      end do ! k
!
!     integral
!
      intsuh(it,1) = intsuh(it,1) + suhh*wg(l,nint)*detjf
      intinh(it,1) = intinh(it,1) + inhh*wg(l,nint)*detjf
      intreh(it,1) = intreh(it,1) + rehh*wg(l,nint)*detjf
      intsuv(it,1) = intsuv(it,1) + suvh*wg(l,nint)*detjf
      intinv(it,1) = intinv(it,1) + invh*wg(l,nint)*detjf
      intctl(it,1) = intctl(it,1) + ctlh*wg(l,nint)*detjf
      intj(it,1) = intj(it,1) + jh*wg(l,nint)*detjf
      ! intjmax(it,1) = intjmax(it,1) + jmaxh*wg(l,nint)*detjf
      intnnc(it,1) = intnnc(it,1) + nnch*wg(l,nint)*detjf
      intvac(it,1) = intvac(it,1) + vach*wg(l,nint)*detjf
!
      end do ! l
!
      end do ! i
!
      end do !it
!
      ! call simpson(intjmax,0.d0,tf,nt+1,jfunc)
      ! write(*,*) '   J at upper bound is', jfunc*10195921.6313475d0*66.67d0
!---------------------------
      write(*,*) '    infected at t = T is', intinh(nt+1,1)*10195921.6313475d0
      write(*,*) '     removed at t = T is', intreh(nt+1,1)*10195921.6313475d0
!
      call simpson(intj,0.d0,tf,nt+1,jfunc)
!
      write(*,*) '   J functional value is', jfunc*10195921.6313475d0*66.67d0
      ! write(123,*) jfunc*10195921.6313475d0*66.67d0
!
      call simpson(intinh,0.d0,tf,nt+1,jfunc)
      write(*,*) '  total # of infected is', jfunc*10195921.6313475d0
!
      call simpson(intnnc,0.d0,tf,nt+1,nnc)
      write(*,*) ' total # of new cases is', nnc*10195921.6313475d0
      ! write(124,*) nnc*10195921.6313475d0
! 
      call simpson(intvac,0.d0,tf,nt+1,vac)
      write(*,*) 'total # of vaccinated is', vac*10195921.6313475d0
!---------------------------- 
      ! write(123,*) jfunc*10195921.6313475d0 , nnc*10195921.6313475d0, vac*10195921.6313475d0
!      
     ! open(unit = 100,file = 'integrais.dat',status = 'unknown')
!
!      write(*,*) 'integrais de suscetiveis e infectados:'
!      write(*,"(15(e15.7))") intsus, intinf
!      write(100,"(15(e15.7))") intsus, intinf
!
      end subroutine
!
!--------------------------------------------------------------------------
!
      subroutine integralad
!
      use mgeometry
      use mgauss
      use marrays
!
      implicit none
!
      real(8), dimension(nsd)         :: xil,xx
      real(8), dimension(nen,nsd)     :: xe
      real(8), dimension(nsd  ,nsd  ) :: jf, invjf, tinvjf
!
      ! real(8) :: phi2d,isuh,iinh,ireh,isuv,iinv
      real(8) :: isuh,iinh,ireh,isuv,iinv
      real(8) :: suhh,inhh,rehh,suvh,invh
      real(8) :: detjf
      integer :: i,j,k,l,ni,nj
!
      ! nint=(5)**nsd
!
      intl1 = 0.d0
      intl2 = 0.d0
      intl3 = 0.d0
      intl4 = 0.d0
      intl5 = 0.d0
!
!     integrating
!
      do it=1,(nt+1)
!
!     loop on elements
!
      do i=1,nelem
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(i)%n(j),k)
      end do ! k
      end do ! j
!
!     loop on integration points
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     Jacobian of the transformation
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1)= jf(2,2)
      invjf(1,2)=-jf(1,2)
      invjf(2,1)=-jf(2,1)
      invjf(2,2)= jf(1,1)
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
!     real coordinates of the integration point
!
      xx=0.d0
      do ni=1,nen
      do nj=1,nsd
      xx(nj)=xx(nj)+phi2d(ni,xil)*xe(ni,nj)
      end do ! nj
      end do ! ni
!
!     approximated solutions at the integration points
!
!      ph=ph+phip(k,xil)*el(i)%p(k,1)!
!     construct of solution in each element:
!
      suhh = 0.d0
      inhh = 0.d0
      rehh = 0.d0
      suvh = 0.d0
      invh = 0.d0
      do k = 1,nen
      suhh = suhh + phi2d(k,xil)*l1(el(i)%n(k),it)
      inhh = inhh + phi2d(k,xil)*l2(el(i)%n(k),it)
      rehh = rehh + phi2d(k,xil)*l3(el(i)%n(k),it)
      suvh = suvh + phi2d(k,xil)*l4(el(i)%n(k),it)
      invh = invh + phi2d(k,xil)*l5(el(i)%n(k),it)
      end do ! k
!
!     integral
!
      intl1(it,1) = intl1(it,1) + suhh*wg(l,nint)*detjf
      intl2(it,1) = intl2(it,1) + inhh*wg(l,nint)*detjf
      intl3(it,1) = intl3(it,1) + rehh*wg(l,nint)*detjf
      intl4(it,1) = intl4(it,1) + suvh*wg(l,nint)*detjf
      intl5(it,1) = intl5(it,1) + invh*wg(l,nint)*detjf
!
      end do ! l
!
      end do ! i
!
      end do ! it
!
!      open(unit = 100,file = 'integrais.dat',status = 'unknown')
!
!      write(*,*) 'integrais de suscetiveis e infectados:'
!      write(*,"(15(e15.7))") intsus, intinf
!      write(100,"(15(e15.7))") intsus, intinf
!
      end subroutine
!
!--------------------------------------------------------------------------
!
      subroutine simpson(f,a,b,n,s)
!
!     f: vetor com os valores da funcao
!     [a,b]: intervalo
!     n: numero de pontos
!
!
      implicit none
!
      integer :: i, j, n
      real(8) :: a, b, h, s
      real(8), dimension(*) :: f
!
      s = f(1) + f(n)
      h = (b-a)/dble(n-1)
!
      do i=1,(n-1)/2
         s = s + 4.d0*f(2*i)
      end do
!
      do i=1,(n-1)/2-1
         s = s + 2.d0*f(2*i+1)
      end do
!
      s = s*h/3.d0
!
!
      end subroutine
!
!--------------------------------------------------------------------------
!
      subroutine elemint(integ,elem,itini,itend)
!
!     calculate integral over one element
!
      use mgeometry
      use marrays
      use mgauss
      use mcoeficientes
!
      implicit none
!
      real(8), dimension(nsd)         :: xil,xx
      real(8), dimension(nen,nsd)     :: xe
      real(8), dimension(nsd  ,nsd  ) :: jf, invjf, tinvjf
      real(8), dimension(int(7/dt)+1,2) :: integ
      integer :: itini, itend
      integer, dimension(*) :: elem
!
      ! real(8) :: phi2d
      real(8) :: inth, inths
      real(8) :: detjf
      integer :: i,j,k,l,ni,nj,nk,iit,ii
!
      ! nint=(5)**nsd
!
      integ = 0.d0
!
!     integrating
!
      do ii=1,1!7
!
!     integration element
!
      i = elem(ii)
! 
      if (i>0) then
!
!     coordinates of the element nodes
!
      do j=1,nen
      do k=1,nsd
      xe(j,k)=coord(el(i)%n(j),k)
      end do ! k
      end do ! j
!
!     loop on integration points
!
      do l=1,nint
!
!     integration points
!
      do k=1,nsd
      xil(k)=xi(l,nint,k)
      end do ! k
!
!     Jacobian of the transformation
!
      call jacobf(xe,xil,jf,detjf)
!
!     Inverse of the Jacobian
!
      invjf(1,1)= jf(2,2)
      invjf(1,2)=-jf(1,2)
      invjf(2,1)=-jf(2,1)
      invjf(2,2)= jf(1,1)
!
      invjf = invjf/detjf
!
!     Transpose
!
      tinvjf = transpose(invjf)
!
!     real coordinates of the integration point
!
      xx=0.d0
      do ni=1,nen
      do nj=1,nsd
      xx(nj)=xx(nj)+phi2d(ni,xil)*xe(ni,nj)
      end do ! nj
      end do ! ni
!
      do it=1,int(7/dt)+1!itini,itend
!
      iit = it+itini-1
!
!     construct of solution in each element:
!
      inth = 0.d0
      inths = 0.d0
      do k = 1,nen
      inth = inth + phi2d(k,xil)*beta*inv(el(i)%n(k),iit)*suh(el(i)%n(k),iit)
      inth = inth + phi2d(k,xil)*betas*inh(el(i)%n(k),iit)*suh(el(i)%n(k),iit)
      end do ! k
!
!     integral
!
      integ(it,1) = integ(it,1) + inth*wg(l,nint)*detjf
      ! integ(it,2) = integ(it,2) + inths*wg(l,nint)*detjf
!
      end do !it
!
      end do ! l
!
      end if 
!
      end do !ii
!
      end subroutine
!
!--------------------------------------------------------------------------
!
      subroutine mosqpdes
!
      use mgeometry
      use marrays
      use mcoeficientes
      use merro
      use mgauss      
!
      implicit none
!
      real(8), dimension(5) :: est
!
!     stability check
!
      est(1) = alpha(1)*dt/(hx*hx)
      est(2) = alpha(2)*dt/(hx*hx)
      est(3) = alpha(3)*dt/(hx*hx)
      est(4) = alpha(4)*dt/(hx*hx)
      est(5) = alpha(5)*dt/(hx*hx)
!
      if (maxval(est)>=0.5d0) then
         write(*,*) maxval(est)
         stop 'stability violated'
      end if
!
!     assemble constant state matrices
!
      call assemble
!
      ! call initial
!
!     solving state equations
!
      call mosquitoes
!
!     plotting solutions
!
      ! call integral
      ! call plotsols
      ! call plotints
!
      end subroutine
!
!--------------------------------------------------------------------------
!
      subroutine contpdes
!
      use mgeometry
      use marrays
      use mcoeficientes
      use merro
      use mgauss      
!
      implicit none
!
      real(8) :: cflag
!
!     assemble constant state matrices
!
      call assemble
!
      ! call initial
      ! call initialad
!
!     sweeping loop
!
      cflag = -1.d0
      iswp = 0
      u = 0.d0
!
      ! u(:,nt1:nt2+1) = umax ! constant control
!
      do while ((cflag<0.d0).and.(iswp<50))
!
      iswp = iswp + 1
!
!     saving old solutions
!
      shold = suh
      ihold = inh
      rhold = reh
      svold = suv
      ivold = inv
      l1old = l1
      l2old = l2
      l3old = l3
      l4old = l4
      l5old = l5
      uold = u
!
!     solving state equations
!
      call mosquitoes
! exit ! constant control
!
!     solving adjoint equations
!
      call adjoint
!
!     updating control
!
      call updtctrl
!
!     testing convergence
!
      call converg(cflag)
!
      end do !while
!
!     plotting solutions
!
      ! call integral
      ! call plotsols
      ! call plotints
!
      end subroutine
!
!--------------------------------------------------------------------------
!
      end module
!
!--------------------------------------------------------------------------
!
      module moptim
!
      implicit none
!
      integer, dimension(:), allocatable :: nst
      real(8), dimension(:), allocatable :: fst
      real(8), dimension(:,:), allocatable :: xst
      integer :: ii, npar
      integer, dimension(20) :: indp
!
      contains
!
!--------------------------------------------------------------------------
!
      subroutine newseed

      implicit none

      integer, dimension(8) :: time
      integer, dimension(:), allocatable :: seed
      real(8) :: r
      integer :: k
! 
      call date_and_time(values=time)
! 
      call random_seed(size=k)
! 
      allocate(seed(1:k))
!
      seed(:) = time(8)
! 
      call random_seed(put=seed)
!
      end subroutine
!
!--------------------------------------------------------------------------
!
      subroutine residual(q,r)
!
!     q: parameters
!     r: residual
!
      use mgeometry
      use mcoeficientes
      use marrays
      use mpdes

      implicit none
!
      integer :: i, j, k, l
      real(8) :: r, aux, auxs, rdata
      real(8), dimension(*) :: q
      real(8), dimension(20) :: paramfit
      ! real(8), dimension(nnodes,nt) :: f
      ! real(8), dimension(nt) :: ut
      real(8), dimension(int(7/dt)+1,2) :: integ
      integer, dimension(21) :: itind
      ! integer, dimension(14,7) :: elem
      integer, dimension(13,1) :: elem=0, indcit=0
      real(8), dimension(13,20) :: obstest = 0.d0
!
!     defining some indexes
!
      do i=1,21!20
      itind(i) = int((i-1)*7/dt+1)
      end do !i
!
      ! open(unit = 400,file = 'elemcities.dat',status = 'old')
      ! elem = 0
      ! do i = 1,14
      !    read(400,*) (elem(i,j),j=1,7)
      ! end do
      ! close(unit=400)
      elem(:,1) = (/132,74,247,66,299,426,156,20,406,24,386,316,69/) !20 weeks
      indcit(1:13,1) = (/1,2,3,4,5,6,7,8,9,10,11,12,13/)
      ! indcit(1:9,1) = (/4,5,6,8,9,10,11,12,13/)
      ! indcit(1:3,1) = (/1,2,7/)
      ! indcit(1:4,1) = (/1,2,3,7/)
!
!     defining parameters
!
!       do i=1,npar
! !
! !     updating params according to indp
!       paramfit(indp(i)) = q(i)
! !
!       end do !i
!
      ! paramfit = params
! 
      ! paramfit(4) = q(1)
      ! paramfit(9) = q(2)
      ! i01 = q(3)
      beta  = q(1) !!
      betav = q(2) !!
      i01   = q(3) !!

      i02 =   q(4) !q(1) !
      i03 =   q(5) !q(2) !
      betas = 0.d0 ! q(6) !
!
      ! call updtparam(paramfit)
!
! beta = 0.4409d0   !0.4409d+0
! betav= 0.003428d0 !0.3428d-2
! i01  = 0.8934d-3  !0.8934d-3
! i02  = 0.5629d-1  !0.d0 !0.5629d-1
! i03  = 0.4576d-1! 0.d0 !0.4576d-1
! betas= 0.7615d-3  !0.d0 !0.7615d-3
!     solving PDEs
!
      nt1 = 1
      nt2 = nt!int(15*7/dt)!
      call initial
      call mosqpdes
!
      ! r = 0.d0
      ! do i = 1,nnodes
      ! do j = 2,(nt+1)
      !    r = r + (inh(i,j)-obsinh(i,j))*(inh(i,j)-obsinh(i,j))
      ! end do !j
      ! end do !i
!
! l=0
      r = 0.d0
      rdata = 0.d0
      do k=1,13!6!14 !cities loop
! 
      i = indcit(k,1)
! 
      do j=1,20!time loop
! 
      if ((obsinh(i,j)>=1.d-14).and.(i>0)) then
      ! if (i>0) then
! l=l+1
! 
      call elemint(integ,elem(i,:),itind(j),itind(j+1))
      ! write(*,*) (integ(k,1),k=1,10)
      call simpson(integ(:,1),(itind(j)-1.d0)*dt,(itind(j+1)-1.d0)*dt,int(7/dt+1),aux)
      ! call simpson(integ(:,2),(itind(j)-1.d0)*dt,(itind(j+1)-1.d0)*dt,int(7/dt+1),auxs)
      ! print *, 'aux', aux
!
      ! r = r + (obsinh(i,j)-aux)*(obsinh(i,j)-aux)*1.0396d+14
      r = r + (obsinh(i,j)-aux)*(obsinh(i,j)-aux)!*dble(j)*7.d0!*1.0396d+14
      rdata = rdata + obsinh(i,j)*obsinh(i,j)
! 
      obstest(i,j) = aux!*dsqrt(1.0396d+14)
      ! obstest(i,j) = auxs*dsqrt(1.0396d+14)
! 
      end if
!
      end do !j
      end do !k
! 
!----------------------------
!       call integral

!       r = 0.d0
!       rdata = 0.d0
!       obstest = 0.d0
! !
!       do j=1,20!time loop
!       call simpson(intnnc(itind(j):itind(j+1),1), &
!             (itind(j)-1.d0)*dt,(itind(j+1)-1.d0)*dt,int(7/dt+1),aux)
! !
!       r = r + (obsinh(1,j)-aux)*(obsinh(1,j)-aux)!*dble(j)*7.d0!*1.0396d+14
!       rdata = rdata + obsinh(1,j)*obsinh(1,j)
! ! 
!       obstest(1,j) = aux!*dsqrt(1.0396d+14)
! ! 
!       end do !j
!----------------------------
! 
      r = dsqrt(r/rdata)
! 
      if ((r > 1.d6).or.(isnan(r)).or.(r /= r)) then
      r = 1.d6
      end if
! 
!     penalization
      ! do i=1,5
      ! r = r + 0.001d0*q(i)*q(i)
      ! end do !i
      ! r = dsqrt(r)
!
! print*, r
!       do i=1,13
!          write(12,*) (obstest(i,j),j=1,20)
!       end do
!       stop
!
      end subroutine
!
!--------------------------------------------------------------------------
!
      subroutine loaddata
!
      use mgeometry
      use marrays
      use mpdes
!
      implicit none
!
      integer :: i, j
!
      character(len = 128) :: filename, expnumaux
!
      write(expnumaux,'(i3.3)') expnum
      ! filename = 'dados/noise01/sol'//trim(expnumaux)//'inh.dat'
      ! filename = 'dados/rndata/rndata2.dat'
      filename = 'dados/noise0/data.dat'
      ! filename = 'dados/noise/noise100.dat'
      open(unit = 400,file = filename,status = 'old')
!
      obsinh = 0.d0
      do i = 1,13!14 !6
         read(400,*) (obsinh(i,j),j=1,20)!19
      end do
      close(unit=400)
!
      ! obsinh = obsinh/10195921.63d0/0.2d0
! 
      end subroutine
!
!--------------------------------------------------------------------------
!
!--------------------------------------------------------------------------
!
      end module
!
!--------------------------------------------------------------------------
!
      module mbobyqa
!
      contains
!
!--------------------------------------------------------------------------
!
      SUBROUTINE ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT, &
        KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,GLAG,HCOL,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(NPT,*),XOPT(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*), &
        SU(*),XNEW(*),XALT(*),GLAG(*),HCOL(*),W(*)
!
!     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have
!       the same meanings as the corresponding arguments of BOBYQB.
!     KOPT is the index of the optimal interpolation point.
!     KNEW is the index of the interpolation point that is going to be moved.
!     ADELT is the current trust region bound.
!     XNEW will be set to a suitable new position for the interpolation point
!       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region
!       bounds and it should provide a large denominator in the next call of
!       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the
!       straight lines through XOPT and another interpolation point.
!     XALT also provides a large value of the modulus of the KNEW-th Lagrange
!       function subject to the constraints that have been mentioned, its main
!       difference from XNEW being that XALT-XOPT is a constrained version of
!       the Cauchy step within the trust region. An exception is that XALT is
!       not calculated if all components of GLAG (see below) are zero.
!     ALPHA will be set to the KNEW-th diagonal element of the H matrix.
!     CAUCHY will be set to the square of the KNEW-th Lagrange function at
!       the step XALT-XOPT from XOPT for the vector XALT that is returned,
!       except that CAUCHY is set to zero if XALT is not calculated.
!     GLAG is a working space vector of length N for the gradient of the
!       KNEW-th Lagrange function at XOPT.
!     HCOL is a working space vector of length NPT for the second derivative
!       coefficients of the KNEW-th Lagrange function.
!     W is a working space vector of length 2N that is going to hold the
!       constrained Cauchy step from XOPT of the Lagrange function, followed
!       by the downhill version of XALT when the uphill step is calculated.
!
!     Set the first NPT components of W to the leading elements of the
!     KNEW-th column of the H matrix.
!
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      CONST=ONE+DSQRT(2.0D0)
      DO 10 K=1,NPT
   10 HCOL(K)=ZERO
      DO 20 J=1,NPT-N-1
      TEMP=ZMAT(KNEW,J)
      DO 20 K=1,NPT
   20 HCOL(K)=HCOL(K)+TEMP*ZMAT(K,J)
      ALPHA=HCOL(KNEW)
      HA=HALF*ALPHA
!
!     Calculate the gradient of the KNEW-th Lagrange function at XOPT.
!
      DO 30 I=1,N
   30 GLAG(I)=BMAT(KNEW,I)
      DO 50 K=1,NPT
      TEMP=ZERO
      DO 40 J=1,N
   40 TEMP=TEMP+XPT(K,J)*XOPT(J)
      TEMP=HCOL(K)*TEMP
      DO 50 I=1,N
   50 GLAG(I)=GLAG(I)+TEMP*XPT(K,I)
!
!     Search for a large denominator along the straight lines through XOPT
!     and another interpolation point. SLBD and SUBD will be lower and upper
!     bounds on the step along each of these lines in turn. PREDSQ will be
!     set to the square of the predicted denominator for each line. PRESAV
!     will be set to the largest admissible value of PREDSQ that occurs.
!
      PRESAV=ZERO
      DO 80 K=1,NPT
      IF (K .EQ. KOPT) GOTO 80
      DDERIV=ZERO
      DISTSQ=ZERO
      DO 60 I=1,N
      TEMP=XPT(K,I)-XOPT(I)
      DDERIV=DDERIV+GLAG(I)*TEMP
   60 DISTSQ=DISTSQ+TEMP*TEMP
      SUBD=ADELT/DSQRT(DISTSQ)
      SLBD=-SUBD
      ILBD=0
      IUBD=0
      SUMIN=DMIN1(ONE,SUBD)
!
!     Revise SLBD and SUBD if necessary because of the bounds in SL and SU.
!
      DO 70 I=1,N
      TEMP=XPT(K,I)-XOPT(I)
      IF (TEMP .GT. ZERO) THEN
          IF (SLBD*TEMP .LT. SL(I)-XOPT(I)) THEN
              SLBD=(SL(I)-XOPT(I))/TEMP
              ILBD=-I
          END IF
          IF (SUBD*TEMP .GT. SU(I)-XOPT(I)) THEN
              SUBD=DMAX1(SUMIN,(SU(I)-XOPT(I))/TEMP)
              IUBD=I
          END IF
      ELSE IF (TEMP .LT. ZERO) THEN
          IF (SLBD*TEMP .GT. SU(I)-XOPT(I)) THEN
              SLBD=(SU(I)-XOPT(I))/TEMP
              ILBD=I
          END IF
          IF (SUBD*TEMP .LT. SL(I)-XOPT(I)) THEN
              SUBD=DMAX1(SUMIN,(SL(I)-XOPT(I))/TEMP)
              IUBD=-I
          END IF
      END IF
   70 CONTINUE
!
!     Seek a large modulus of the KNEW-th Lagrange function when the index
!     of the other interpolation point on the line through XOPT is KNEW.
!
      IF (K .EQ. KNEW) THEN
          DIFF=DDERIV-ONE
          STEP=SLBD
          VLAG=SLBD*(DDERIV-SLBD*DIFF)
          ISBD=ILBD
          TEMP=SUBD*(DDERIV-SUBD*DIFF)
          IF (DABS(TEMP) .GT. DABS(VLAG)) THEN
              STEP=SUBD
              VLAG=TEMP
              ISBD=IUBD
          END IF
          TEMPD=HALF*DDERIV
          TEMPA=TEMPD-DIFF*SLBD
          TEMPB=TEMPD-DIFF*SUBD
          IF (TEMPA*TEMPB .LT. ZERO) THEN
              TEMP=TEMPD*TEMPD/DIFF
              IF (DABS(TEMP) .GT. DABS(VLAG)) THEN
                  STEP=TEMPD/DIFF
                  VLAG=TEMP
                  ISBD=0
              END IF
          END IF
!
!     Search along each of the other lines through XOPT and another point.
!
      ELSE
          STEP=SLBD
          VLAG=SLBD*(ONE-SLBD)
          ISBD=ILBD
          TEMP=SUBD*(ONE-SUBD)
          IF (DABS(TEMP) .GT. DABS(VLAG)) THEN
              STEP=SUBD
              VLAG=TEMP
              ISBD=IUBD
          END IF
          IF (SUBD .GT. HALF) THEN
              IF (DABS(VLAG) .LT. 0.25D0) THEN
                  STEP=HALF
                  VLAG=0.25D0
                  ISBD=0
              END IF
          END IF
          VLAG=VLAG*DDERIV
      END IF
!
!     Calculate PREDSQ for the current line search and maintain PRESAV.
!
      TEMP=STEP*(ONE-STEP)*DISTSQ
      PREDSQ=VLAG*VLAG*(VLAG*VLAG+HA*TEMP*TEMP)
      IF (PREDSQ .GT. PRESAV) THEN
          PRESAV=PREDSQ
          KSAV=K
          STPSAV=STEP
          IBDSAV=ISBD
      END IF
   80 CONTINUE
!
!     Construct XNEW in a way that satisfies the bound constraints exactly.
!
      DO 90 I=1,N
      TEMP=XOPT(I)+STPSAV*(XPT(KSAV,I)-XOPT(I))
   90 XNEW(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
      IF (IBDSAV .LT. 0) XNEW(-IBDSAV)=SL(-IBDSAV)
      IF (IBDSAV .GT. 0) XNEW(IBDSAV)=SU(IBDSAV)
!
!     Prepare for the iterative method that assembles the constrained Cauchy
!     step in W. The sum of squares of the fixed components of W is formed in
!     WFIXSQ, and the free components of W are set to BIGSTP.
!
      BIGSTP=ADELT+ADELT
      IFLAG=0
  100 WFIXSQ=ZERO
      GGFREE=ZERO
      DO 110 I=1,N
      W(I)=ZERO
      TEMPA=DMIN1(XOPT(I)-SL(I),GLAG(I))
      TEMPB=DMAX1(XOPT(I)-SU(I),GLAG(I))
      IF (TEMPA .GT. ZERO .OR. TEMPB .LT. ZERO) THEN
          W(I)=BIGSTP
          GGFREE=GGFREE+GLAG(I)**2
      END IF
  110 CONTINUE
      IF (GGFREE .EQ. ZERO) THEN
          CAUCHY=ZERO
          GOTO 200
      END IF
!
!     Investigate whether more components of W can be fixed.
!
  120 TEMP=ADELT*ADELT-WFIXSQ
      IF (TEMP .GT. ZERO) THEN
          WSQSAV=WFIXSQ
          STEP=DSQRT(TEMP/GGFREE)
          GGFREE=ZERO
          DO 130 I=1,N
          IF (W(I) .EQ. BIGSTP) THEN
              TEMP=XOPT(I)-STEP*GLAG(I)
              IF (TEMP .LE. SL(I)) THEN
                  W(I)=SL(I)-XOPT(I)
                  WFIXSQ=WFIXSQ+W(I)**2
              ELSE IF (TEMP .GE. SU(I)) THEN
                  W(I)=SU(I)-XOPT(I)
                  WFIXSQ=WFIXSQ+W(I)**2
              ELSE
                  GGFREE=GGFREE+GLAG(I)**2
              END IF
          END IF
  130     CONTINUE
          IF (WFIXSQ .GT. WSQSAV .AND. GGFREE .GT. ZERO) GOTO 120
      END IF
!
!     Set the remaining free components of W and all components of XALT,
!     except that W may be scaled later.
!
      GW=ZERO
      DO 140 I=1,N
      IF (W(I) .EQ. BIGSTP) THEN
          W(I)=-STEP*GLAG(I)
          XALT(I)=DMAX1(SL(I),DMIN1(SU(I),XOPT(I)+W(I)))
      ELSE IF (W(I) .EQ. ZERO) THEN
          XALT(I)=XOPT(I)
      ELSE IF (GLAG(I) .GT. ZERO) THEN
          XALT(I)=SL(I)
      ELSE
          XALT(I)=SU(I)
      END IF
  140 GW=GW+GLAG(I)*W(I)
!
!     Set CURV to the curvature of the KNEW-th Lagrange function along W.
!     Scale W by a factor less than one if that can reduce the modulus of
!     the Lagrange function at XOPT+W. Set CAUCHY to the final value of
!     the square of this function.
!
      CURV=ZERO
      DO 160 K=1,NPT
      TEMP=ZERO
      DO 150 J=1,N
  150 TEMP=TEMP+XPT(K,J)*W(J)
  160 CURV=CURV+HCOL(K)*TEMP*TEMP
      IF (IFLAG .EQ. 1) CURV=-CURV
      IF (CURV .GT. -GW .AND. CURV .LT. -CONST*GW) THEN
          SCALE=-GW/CURV
          DO 170 I=1,N
          TEMP=XOPT(I)+SCALE*W(I)
  170     XALT(I)=DMAX1(SL(I),DMIN1(SU(I),TEMP))
          CAUCHY=(HALF*GW*SCALE)**2
      ELSE
          CAUCHY=(GW+HALF*CURV)**2
      END IF
!
!     If IFLAG is zero, then XALT is calculated as before after reversing
!     the sign of GLAG. Thus two XALT vectors become available. The one that
!     is chosen is the one that gives the larger value of CAUCHY.
!
      IF (IFLAG .EQ. 0) THEN
          DO 180 I=1,N
          GLAG(I)=-GLAG(I)
  180     W(N+I)=XALT(I)
          CSAVE=CAUCHY
          IFLAG=1
          GOTO 100
      END IF
      IF (CSAVE .GT. CAUCHY) THEN
          DO 190 I=1,N
  190     XALT(I)=W(N+I)
          CAUCHY=CSAVE
      END IF
  200 RETURN
      END

!--------------------------------------------------------------------------

      SUBROUTINE BOBYQA (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),W(*)
!
!     This subroutine seeks the least value of a function of many variables,
!     by applying a trust region method that forms quadratic models by
!     interpolation. There is usually some freedom in the interpolation
!     conditions, which is taken up by minimizing the Frobenius norm of
!     the change to the second derivative of the model, beginning with the
!     zero matrix. The values of the variables are constrained by upper and
!     lower bounds. The arguments of the subroutine are as follows.
!
!     N must be set to the number of variables and must be at least two.
!     NPT is the number of interpolation conditions. Its value must be in
!       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not
!       recommended.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper
!       bounds, respectively, on X(I). The construction of quadratic models
!       requires XL(I) to be strictly less than XU(I) for each I. Further,
!       the contribution to a model from changes to the I-th variable is
!       damaged severely by rounding errors if XU(I)-XL(I) is too small.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND no greater than
!       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest
!       expected change to a variable, while RHOEND should indicate the
!       accuracy that is required in the final values of the variables. An
!       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N,
!       is less than 2*RHOBEG.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing. Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1. Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       (NPT+5)*(NPT+N)+3*N*(N+5)/2.
!
!     SUBROUTINE CALFUN (N,X,F) has to be provided by the user. It must set
!     F to the value of the objective function for the current values of the
!     variables X(1),X(2),...,X(N), which are generated automatically in a
!     way that satisfies the bounds given in XL and XU.
!
!     Return if the value of NPT is unacceptable.
!
      NP=N+1
      IF (NPT .LT. N+2 .OR. NPT .GT. ((N+2)*NP)/2) THEN
          PRINT 10
   10     FORMAT (/4X,'Return from BOBYQA because NPT is not in', &
            ' the required interval')
          GO TO 40
      END IF
!
!     Partition the working space array, so that different parts of it can
!     be treated separately during the calculation of BOBYQB. The partition
!     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the
!     space that is taken by the last array in the argument list of BOBYQB.
!
      NDIM=NPT+N
      IXB=1
      IXP=IXB+N
      IFV=IXP+N*NPT
      IXO=IFV+NPT
      IGO=IXO+N
      IHQ=IGO+N
      IPQ=IHQ+(N*NP)/2
      IBMAT=IPQ+NPT
      IZMAT=IBMAT+NDIM*N
      ISL=IZMAT+NPT*(NPT-NP)
      ISU=ISL+N
      IXN=ISU+N
      IXA=IXN+N
      ID=IXA+N
      IVL=ID+N
      IW=IVL+NDIM
!
!     Return if there is insufficient space between the bounds. Modify the
!     initial X if necessary in order to avoid conflicts between the bounds
!     and the construction of the first quadratic model. The lower and upper
!     bounds on moves from the updated X are set now, in the ISL and ISU
!     partitions of W, in order to provide useful and exact information about
!     components of X that become within distance RHOBEG from their bounds.
!
      ZERO=0.0D0
      DO 30 J=1,N
      TEMP=XU(J)-XL(J)
      IF (TEMP .LT. RHOBEG+RHOBEG) THEN
          PRINT 20
   20     FORMAT (/4X,'Return from BOBYQA because one of the', &
            ' differences XU(I)-XL(I)'/6X,' is less than 2*RHOBEG.')
          GO TO 40
      END IF
      JSL=ISL+J-1
      JSU=JSL+N
      W(JSL)=XL(J)-X(J)
      W(JSU)=XU(J)-X(J)
      IF (W(JSL) .GE. -RHOBEG) THEN
          IF (W(JSL) .GE. ZERO) THEN
              X(J)=XL(J)
              W(JSL)=ZERO
              W(JSU)=TEMP
          ELSE
              X(J)=XL(J)+RHOBEG
              W(JSL)=-RHOBEG
              W(JSU)=DMAX1(XU(J)-X(J),RHOBEG)
          END IF
      ELSE IF (W(JSU) .LE. RHOBEG) THEN
          IF (W(JSU) .LE. ZERO) THEN
              X(J)=XU(J)
              W(JSL)=-TEMP
              W(JSU)=ZERO
          ELSE
              X(J)=XU(J)-RHOBEG
              W(JSL)=DMIN1(XL(J)-X(J),-RHOBEG)
              W(JSU)=RHOBEG
          END IF
      END IF
   30 CONTINUE
!
!     Make the call of BOBYQB.
!
      CALL BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT,MAXFUN,W(IXB), &
        W(IXP),W(IFV),W(IXO),W(IGO),W(IHQ),W(IPQ),W(IBMAT),W(IZMAT), &
        NDIM,W(ISL),W(ISU),W(IXN),W(IXA),W(ID),W(IVL),W(IW))
   40 RETURN
      END

!--------------------------------------------------------------------------

      SUBROUTINE BOBYQB (N,NPT,X,XL,XU,RHOBEG,RHOEND,IPRINT, &
        MAXFUN,XBASE,XPT,FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM, &
        SL,SU,XNEW,XALT,D,VLAG,W)
!
      use moptim, only : nst,fst,xst,ii
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*), &
        XOPT(*),GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*), &
        SL(*),SU(*),XNEW(*),XALT(*),D(*),VLAG(*),W(*)
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN
!       are identical to the corresponding arguments in SUBROUTINE BOBYQA.
!     XBASE holds a shift of origin that should reduce the contributions
!       from rounding errors to values of the model and Lagrange functions.
!     XPT is a two-dimensional array that holds the coordinates of the
!       interpolation points relative to XBASE.
!     FVAL holds the values of F at the interpolation points.
!     XOPT is set to the displacement from XBASE of the trust region centre.
!     GOPT holds the gradient of the quadratic model at XBASE+XOPT.
!     HQ holds the explicit second derivatives of the quadratic model.
!     PQ contains the parameters of the implicit second derivatives of the
!       quadratic model.
!     BMAT holds the last N columns of H.
!     ZMAT holds the factorization of the leading NPT by NPT submatrix of H,
!       this factorization being ZMAT times ZMAT^T, which provides both the
!       correct rank and positive semi-definiteness.
!     NDIM is the first dimension of BMAT and has the value NPT+N.
!     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively.
!       All the components of every XOPT are going to satisfy the bounds
!       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when
!       XOPT is on a constraint boundary.
!     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the
!       vector of variables for the next call of CALFUN. XNEW also satisfies
!       the SL and SU constraints in the way that has just been mentioned.
!     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW
!       in order to increase the denominator in the updating of UPDATE.
!     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT.
!     VLAG contains the values of the Lagrange functions at a new point X.
!       They are part of a product that requires VLAG to be of length NDIM.
!     W is a one-dimensional array that is used for working space. Its length
!       must be at least 3*NDIM = 3*(NPT+N).
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      TEN=10.0D0
      TENTH=0.1D0
      TWO=2.0D0
      ZERO=0.0D0
      NP=N+1
      NPTM=NPT-NP
      NH=(N*NP)/2
!
!     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, with the corresponding values of
!     of NF and KOPT, which are the number of calls of CALFUN so far and the
!     index of the interpolation point at the trust region centre. Then the
!     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is
!     less than NPT. GOPT will be updated if KOPT is different from KBASE.
!
      CALL PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE,XPT, &
        FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
      XOPTSQ=ZERO
      DO 10 I=1,N
      XOPT(I)=XPT(KOPT,I)
   10 XOPTSQ=XOPTSQ+XOPT(I)**2
      FSAVE=FVAL(1)
      IF (NF .LT. NPT) THEN
          IF (IPRINT .GT. 0) PRINT 390
          GOTO 720
      END IF
      KBASE=1
!
!     Complete the settings that are required for the iterative procedure.
!
      RHO=RHOBEG
      DELTA=RHO
      NRESC=NF
      NTRITS=0
      DIFFA=ZERO
      DIFFB=ZERO
      ITEST=0
      NFSAV=NF
!
!     Update GOPT if necessary before the first iteration and after each
!     call of RESCUE that makes a call of CALFUN.
!
   20 IF (KOPT .NE. KBASE) THEN
          IH=0
          DO 30 J=1,N
          DO 30 I=1,J
          IH=IH+1
          IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*XOPT(I)
   30     GOPT(I)=GOPT(I)+HQ(IH)*XOPT(J)
          IF (NF .GT. NPT) THEN
              DO 50 K=1,NPT
              TEMP=ZERO
              DO 40 J=1,N
   40         TEMP=TEMP+XPT(K,J)*XOPT(J)
              TEMP=PQ(K)*TEMP
              DO 50 I=1,N
   50         GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
          END IF
      END IF
!
!     Generate the next point in the trust region that provides a small value
!     of the quadratic model subject to the constraints on the variables.
!     The integer NTRITS is set to the number "trust region" iterations that
!     have occurred since the last "alternative" iteration. If the length
!     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to
!     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW.
!
   60 CALL TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA,XNEW,D, &
        W,W(NP),W(NP+N),W(NP+2*N),W(NP+3*N),DSQ,CRVMIN)
      DNORM=DMIN1(DELTA,DSQRT(DSQ))
      IF (DNORM .LT. HALF*RHO) THEN
          NTRITS=-1
          DISTSQ=(TEN*RHO)**2
          IF (NF .LE. NFSAV+2) GOTO 650
!
!     The following choice between labels 650 and 680 depends on whether or
!     not our work with the current RHO seems to be complete. Either RHO is
!     decreased or termination occurs if the errors in the quadratic model at
!     the last three interpolation points compare favourably with predictions
!     of likely improvements to the model within distance HALF*RHO of XOPT.
!
          ERRBIG=DMAX1(DIFFA,DIFFB,DIFFC)
          FRHOSQ=0.125D0*RHO*RHO
          IF (CRVMIN .GT. ZERO .AND. ERRBIG .GT. FRHOSQ*CRVMIN) &
             GOTO 650
          BDTOL=ERRBIG/RHO
          DO 80 J=1,N
          BDTEST=BDTOL
          IF (XNEW(J) .EQ. SL(J)) BDTEST=W(J)
          IF (XNEW(J) .EQ. SU(J)) BDTEST=-W(J)
          IF (BDTEST .LT. BDTOL) THEN
              CURV=HQ((J+J*J)/2)
              DO 70 K=1,NPT
   70         CURV=CURV+PQ(K)*XPT(K,J)**2
              BDTEST=BDTEST+HALF*CURV*RHO
              IF (BDTEST .LT. BDTOL) GOTO 650
          END IF
   80     CONTINUE
          GOTO 680
      END IF
      NTRITS=NTRITS+1
!
!     Severe cancellation is likely to occur if XOPT is too far from XBASE.
!     If the following test holds, then XBASE is shifted so that XOPT becomes
!     zero. The appropriate changes are made to BMAT and to the second
!     derivatives of the current model, beginning with the changes to BMAT
!     that do not depend on ZMAT. VLAG is used temporarily for working space.
!
   90 IF (DSQ .LE. 1.0D-3*XOPTSQ) THEN
          FRACSQ=0.25D0*XOPTSQ
          SUMPQ=ZERO
          DO 110 K=1,NPT
          SUMPQ=SUMPQ+PQ(K)
          SUM=-HALF*XOPTSQ
          DO 100 I=1,N
  100     SUM=SUM+XPT(K,I)*XOPT(I)
          W(NPT+K)=SUM
          TEMP=FRACSQ-HALF*SUM
          DO 110 I=1,N
          W(I)=BMAT(K,I)
          VLAG(I)=SUM*XPT(K,I)+TEMP*XOPT(I)
          IP=NPT+I
          DO 110 J=1,I
  110     BMAT(IP,J)=BMAT(IP,J)+W(I)*VLAG(J)+VLAG(I)*W(J)
!
!     Then the revisions of BMAT that depend on ZMAT are calculated.
!
          DO 150 JJ=1,NPTM
          SUMZ=ZERO
          SUMW=ZERO
          DO 120 K=1,NPT
          SUMZ=SUMZ+ZMAT(K,JJ)
          VLAG(K)=W(NPT+K)*ZMAT(K,JJ)
  120     SUMW=SUMW+VLAG(K)
          DO 140 J=1,N
          SUM=(FRACSQ*SUMZ-HALF*SUMW)*XOPT(J)
          DO 130 K=1,NPT
  130     SUM=SUM+VLAG(K)*XPT(K,J)
          W(J)=SUM
          DO 140 K=1,NPT
  140     BMAT(K,J)=BMAT(K,J)+SUM*ZMAT(K,JJ)
          DO 150 I=1,N
          IP=I+NPT
          TEMP=W(I)
          DO 150 J=1,I
  150     BMAT(IP,J)=BMAT(IP,J)+TEMP*W(J)
!
!     The following instructions complete the shift, including the changes
!     to the second derivative parameters of the quadratic model.
!
          IH=0
          DO 170 J=1,N
          W(J)=-HALF*SUMPQ*XOPT(J)
          DO 160 K=1,NPT
          W(J)=W(J)+PQ(K)*XPT(K,J)
  160     XPT(K,J)=XPT(K,J)-XOPT(J)
          DO 170 I=1,J
          IH=IH+1
          HQ(IH)=HQ(IH)+W(I)*XOPT(J)+XOPT(I)*W(J)
  170     BMAT(NPT+I,J)=BMAT(NPT+J,I)
          DO 180 I=1,N
          XBASE(I)=XBASE(I)+XOPT(I)
          XNEW(I)=XNEW(I)-XOPT(I)
          SL(I)=SL(I)-XOPT(I)
          SU(I)=SU(I)-XOPT(I)
  180     XOPT(I)=ZERO
          XOPTSQ=ZERO
      END IF
      IF (NTRITS .EQ. 0) GOTO 210
      GOTO 230
!
!     XBASE is also moved to XOPT by a call of RESCUE. This calculation is
!     more expensive than the previous shift, because new matrices BMAT and
!     ZMAT are generated from scratch, which may include the replacement of
!     interpolation points whose positions seem to be causing near linear
!     dependence in the interpolation conditions. Therefore RESCUE is called
!     only if rounding errors have reduced by at least a factor of two the
!     denominator of the formula for updating the H matrix. It provides a
!     useful safeguard, but is not invoked in most applications of BOBYQA.
!
  190 NFSAV=NF
      KBASE=KOPT
      CALL RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT,FVAL, &
        XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA,KOPT, &
        VLAG,W,W(N+NP),W(NDIM+NP))
!
!     XOPT is updated now in case the branch below to label 720 is taken.
!     Any updating of GOPT occurs after the branch below to label 20, which
!     leads to a trust region iteration as does the branch to label 60.
!
      XOPTSQ=ZERO
      IF (KOPT .NE. KBASE) THEN
          DO 200 I=1,N
          XOPT(I)=XPT(KOPT,I)
  200     XOPTSQ=XOPTSQ+XOPT(I)**2
      END IF
      IF (NF .LT. 0) THEN
          NF=MAXFUN
          IF (IPRINT .GT. 0) PRINT 390
          GOTO 720
      END IF
      NRESC=NF
      IF (NFSAV .LT. NF) THEN
          NFSAV=NF
          GOTO 20
      END IF
      IF (NTRITS .GT. 0) GOTO 60
!
!     Pick two alternative vectors of variables, relative to XBASE, that
!     are suitable as new positions of the KNEW-th interpolation point.
!     Firstly, XNEW is set to the point on a line through XOPT and another
!     interpolation point that minimizes the predicted value of the next
!     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL
!     and SU bounds. Secondly, XALT is set to the best feasible point on
!     a constrained version of the Cauchy step of the KNEW-th Lagrange
!     function, the corresponding value of the square of this function
!     being returned in CAUCHY. The choice between these alternatives is
!     going to be made when the denominator is calculated.
!
  210 CALL ALTMOV (N,NPT,XPT,XOPT,BMAT,ZMAT,NDIM,SL,SU,KOPT, &
        KNEW,ADELT,XNEW,XALT,ALPHA,CAUCHY,W,W(NP),W(NDIM+1))
      DO 220 I=1,N
  220 D(I)=XNEW(I)-XOPT(I)
!
!     Calculate VLAG and BETA for the current choice of D. The scalar
!     product of D with XPT(K,.) is going to be held in W(NPT+K) for
!     use when VQUAD is calculated.
!
  230 DO 250 K=1,NPT
      SUMA=ZERO
      SUMB=ZERO
      SUM=ZERO
      DO 240 J=1,N
      SUMA=SUMA+XPT(K,J)*D(J)
      SUMB=SUMB+XPT(K,J)*XOPT(J)
  240 SUM=SUM+BMAT(K,J)*D(J)
      W(K)=SUMA*(HALF*SUMA+SUMB)
      VLAG(K)=SUM
  250 W(NPT+K)=SUMA
      BETA=ZERO
      DO 270 JJ=1,NPTM
      SUM=ZERO
      DO 260 K=1,NPT
  260 SUM=SUM+ZMAT(K,JJ)*W(K)
      BETA=BETA-SUM*SUM
      DO 270 K=1,NPT
  270 VLAG(K)=VLAG(K)+SUM*ZMAT(K,JJ)
      DSQ=ZERO
      BSUM=ZERO
      DX=ZERO
      DO 300 J=1,N
      DSQ=DSQ+D(J)**2
      SUM=ZERO
      DO 280 K=1,NPT
  280 SUM=SUM+W(K)*BMAT(K,J)
      BSUM=BSUM+SUM*D(J)
      JP=NPT+J
      DO 290 I=1,N
  290 SUM=SUM+BMAT(JP,I)*D(I)
      VLAG(JP)=SUM
      BSUM=BSUM+SUM*D(J)
  300 DX=DX+D(J)*XOPT(J)
      BETA=DX*DX+DSQ*(XOPTSQ+DX+DX+HALF*DSQ)+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
!
!     If NTRITS is zero, the denominator may be increased by replacing
!     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if
!     rounding errors have damaged the chosen denominator.
!
      IF (NTRITS .EQ. 0) THEN
          DENOM=VLAG(KNEW)**2+ALPHA*BETA
          IF (DENOM .LT. CAUCHY .AND. CAUCHY .GT. ZERO) THEN
              DO 310 I=1,N
              XNEW(I)=XALT(I)
  310         D(I)=XNEW(I)-XOPT(I)
              CAUCHY=ZERO
              GO TO 230
          END IF
          IF (DENOM .LE. HALF*VLAG(KNEW)**2) THEN
              IF (NF .GT. NRESC) GOTO 190
              IF (IPRINT .GT. 0) PRINT 320
  320         FORMAT (/5X,'Return from BOBYQA because of much', &
                ' cancellation in a denominator.')
              GOTO 720
          END IF
!
!     Alternatively, if NTRITS is positive, then set KNEW to the index of
!     the next interpolation point to be deleted to make room for a trust
!     region step. Again RESCUE may be called if rounding errors have damaged
!     the chosen denominator, which is the reason for attempting to select
!     KNEW before calculating the next value of the objective function.
!
      ELSE
          DELSQ=DELTA*DELTA
          SCADEN=ZERO
          BIGLSQ=ZERO
          KNEW=0
          DO 350 K=1,NPT
          IF (K .EQ. KOPT) GOTO 350
          HDIAG=ZERO
          DO 330 JJ=1,NPTM
  330     HDIAG=HDIAG+ZMAT(K,JJ)**2
          DEN=BETA*HDIAG+VLAG(K)**2
          DISTSQ=ZERO
          DO 340 J=1,N
  340     DISTSQ=DISTSQ+(XPT(K,J)-XOPT(J))**2
          TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
          IF (TEMP*DEN .GT. SCADEN) THEN
              SCADEN=TEMP*DEN
              KNEW=K
              DENOM=DEN
          END IF
          BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
  350     CONTINUE
          IF (SCADEN .LE. HALF*BIGLSQ) THEN
              IF (NF .GT. NRESC) GOTO 190
              IF (IPRINT .GT. 0) PRINT 320
              GOTO 720
          END IF
      END IF
!
!     Put the variables for the next calculation of the objective function
!       in XNEW, with any adjustments for the bounds.
!
!
!     Calculate the value of the objective function at XBASE+XNEW, unless
!       the limit on the number of calculations of F has been reached.
!
  360 DO 380 I=1,N
      X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XNEW(I)),XU(I))
      IF (XNEW(I) .EQ. SL(I)) X(I)=XL(I)
      IF (XNEW(I) .EQ. SU(I)) X(I)=XU(I)
  380 CONTINUE
      IF (NF .GE. MAXFUN) THEN
          IF (IPRINT .GT. 0) PRINT 390
  390     FORMAT (/4X,'Return from BOBYQA because CALFUN has been', &
            ' called MAXFUN times.')
          GOTO 720
      END IF
      NF=NF+1
      CALL CALFUN (N,X,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 400, NF,F,(X(I),I=1,N)
  400      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10, &
             '    The corresponding X is:'/(2X,5D15.6))
      END IF
      IF (NTRITS .EQ. -1) THEN
          FSAVE=F
          GOTO 720
      END IF
!
!     Use the quadratic model to predict the change in F due to the step D,
!       and set DIFF to the error of this prediction.
!
      FOPT=FVAL(KOPT)
      VQUAD=ZERO
      IH=0
      DO 410 J=1,N
      VQUAD=VQUAD+D(J)*GOPT(J)
      DO 410 I=1,J
      IH=IH+1
      TEMP=D(I)*D(J)
      IF (I .EQ. J) TEMP=HALF*TEMP
  410 VQUAD=VQUAD+HQ(IH)*TEMP
      DO 420 K=1,NPT
  420 VQUAD=VQUAD+HALF*PQ(K)*W(NPT+K)**2
      DIFF=F-FOPT-VQUAD
      DIFFC=DIFFB
      DIFFB=DIFFA
      DIFFA=DABS(DIFF)
      IF (DNORM .GT. RHO) NFSAV=NF
!
!     Pick the next value of DELTA after a trust region step.
!
      IF (NTRITS .GT. 0) THEN
          IF (VQUAD .GE. ZERO) THEN
              IF (IPRINT .GT. 0) PRINT 430
  430         FORMAT (/4X,'Return from BOBYQA because a trust', &
                ' region step has failed to reduce Q.')
              GOTO 720
          END IF
          RATIO=(F-FOPT)/VQUAD
          IF (RATIO .LE. TENTH) THEN
              DELTA=DMIN1(HALF*DELTA,DNORM)
          ELSE IF (RATIO .LE. 0.7D0) THEN
              DELTA=DMAX1(HALF*DELTA,DNORM)
          ELSE
              DELTA=DMAX1(HALF*DELTA,DNORM+DNORM)
          END IF
          IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
!
!     Recalculate KNEW and DENOM if the new F is less than FOPT.
!
          IF (F .LT. FOPT) THEN
              KSAV=KNEW
              DENSAV=DENOM
              DELSQ=DELTA*DELTA
              SCADEN=ZERO
              BIGLSQ=ZERO
              KNEW=0
              DO 460 K=1,NPT
              HDIAG=ZERO
              DO 440 JJ=1,NPTM
  440         HDIAG=HDIAG+ZMAT(K,JJ)**2
              DEN=BETA*HDIAG+VLAG(K)**2
              DISTSQ=ZERO
              DO 450 J=1,N
  450         DISTSQ=DISTSQ+(XPT(K,J)-XNEW(J))**2
              TEMP=DMAX1(ONE,(DISTSQ/DELSQ)**2)
              IF (TEMP*DEN .GT. SCADEN) THEN
                  SCADEN=TEMP*DEN
                  KNEW=K
                  DENOM=DEN
              END IF
  460         BIGLSQ=DMAX1(BIGLSQ,TEMP*VLAG(K)**2)
              IF (SCADEN .LE. HALF*BIGLSQ) THEN
                  KNEW=KSAV
                  DENOM=DENSAV
              END IF
          END IF
      END IF
!
!     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be
!     moved. Also update the second derivative terms of the model.
!
      CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
      IH=0
      PQOLD=PQ(KNEW)
      PQ(KNEW)=ZERO
      DO 470 I=1,N
      TEMP=PQOLD*XPT(KNEW,I)
      DO 470 J=1,I
      IH=IH+1
  470 HQ(IH)=HQ(IH)+TEMP*XPT(KNEW,J)
      DO 480 JJ=1,NPTM
      TEMP=DIFF*ZMAT(KNEW,JJ)
      DO 480 K=1,NPT
  480 PQ(K)=PQ(K)+TEMP*ZMAT(K,JJ)
!
!     Include the new interpolation point, and make the changes to GOPT at
!     the old XOPT that are caused by the updating of the quadratic model.
!
      FVAL(KNEW)=F
      DO 490 I=1,N
      XPT(KNEW,I)=XNEW(I)
  490 W(I)=BMAT(KNEW,I)
      DO 520 K=1,NPT
      SUMA=ZERO
      DO 500 JJ=1,NPTM
  500 SUMA=SUMA+ZMAT(KNEW,JJ)*ZMAT(K,JJ)
      SUMB=ZERO
      DO 510 J=1,N
  510 SUMB=SUMB+XPT(K,J)*XOPT(J)
      TEMP=SUMA*SUMB
      DO 520 I=1,N
  520 W(I)=W(I)+TEMP*XPT(K,I)
      DO 530 I=1,N
  530 GOPT(I)=GOPT(I)+DIFF*W(I)
!
!     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT.
!
      IF (F .LT. FOPT) THEN
          KOPT=KNEW
          XOPTSQ=ZERO
          IH=0
          DO 540 J=1,N
          XOPT(J)=XNEW(J)
          XOPTSQ=XOPTSQ+XOPT(J)**2
          DO 540 I=1,J
          IH=IH+1
          IF (I .LT. J) GOPT(J)=GOPT(J)+HQ(IH)*D(I)
  540     GOPT(I)=GOPT(I)+HQ(IH)*D(J)
          DO 560 K=1,NPT
          TEMP=ZERO
          DO 550 J=1,N
  550     TEMP=TEMP+XPT(K,J)*D(J)
          TEMP=PQ(K)*TEMP
          DO 560 I=1,N
  560     GOPT(I)=GOPT(I)+TEMP*XPT(K,I)
      END IF
!
!     Calculate the parameters of the least Frobenius norm interpolant to
!     the current data, the gradient of this interpolant at XOPT being put
!     into VLAG(NPT+I), I=1,2,...,N.
!
      IF (NTRITS .GT. 0) THEN
          DO 570 K=1,NPT
          VLAG(K)=FVAL(K)-FVAL(KOPT)
  570     W(K)=ZERO
          DO 590 J=1,NPTM
          SUM=ZERO
          DO 580 K=1,NPT
  580     SUM=SUM+ZMAT(K,J)*VLAG(K)
          DO 590 K=1,NPT
  590     W(K)=W(K)+SUM*ZMAT(K,J)
          DO 610 K=1,NPT
          SUM=ZERO
          DO 600 J=1,N
  600     SUM=SUM+XPT(K,J)*XOPT(J)
          W(K+NPT)=W(K)
  610     W(K)=SUM*W(K)
          GQSQ=ZERO
          GISQ=ZERO
          DO 630 I=1,N
          SUM=ZERO
          DO 620 K=1,NPT
  620     SUM=SUM+BMAT(K,I)*VLAG(K)+XPT(K,I)*W(K)
          IF (XOPT(I) .EQ. SL(I)) THEN
              GQSQ=GQSQ+DMIN1(ZERO,GOPT(I))**2
              GISQ=GISQ+DMIN1(ZERO,SUM)**2
          ELSE IF (XOPT(I) .EQ. SU(I)) THEN
              GQSQ=GQSQ+DMAX1(ZERO,GOPT(I))**2
              GISQ=GISQ+DMAX1(ZERO,SUM)**2
          ELSE
              GQSQ=GQSQ+GOPT(I)**2
              GISQ=GISQ+SUM*SUM
          END IF
  630     VLAG(NPT+I)=SUM
!
!     Test whether to replace the new quadratic model by the least Frobenius
!     norm interpolant, making the replacement if the test is satisfied.
!
          ITEST=ITEST+1
          IF (GQSQ .LT. TEN*GISQ) ITEST=0
          IF (ITEST .GE. 3) THEN
              DO 640 I=1,MAX0(NPT,NH)
              IF (I .LE. N) GOPT(I)=VLAG(NPT+I)
              IF (I .LE. NPT) PQ(I)=W(NPT+I)
              IF (I .LE. NH) HQ(I)=ZERO
              ITEST=0
  640         CONTINUE
          END IF
      END IF
!
!     If a trust region step has provided a sufficient decrease in F, then
!     branch for another trust region calculation. The case NTRITS=0 occurs
!     when the new interpolation point was reached by an alternative step.
!
      IF (NTRITS .EQ. 0) GOTO 60
      IF (F .LE. FOPT+TENTH*VQUAD) GOTO 60
!
!     Alternatively, find out if the interpolation points are close enough
!       to the best point so far.
!
      DISTSQ=DMAX1((TWO*DELTA)**2,(TEN*RHO)**2)
  650 KNEW=0
      DO 670 K=1,NPT
      SUM=ZERO
      DO 660 J=1,N
  660 SUM=SUM+(XPT(K,J)-XOPT(J))**2
      IF (SUM .GT. DISTSQ) THEN
          KNEW=K
          DISTSQ=SUM
      END IF
  670 CONTINUE
!
!     If KNEW is positive, then ALTMOV finds alternative new positions for
!     the KNEW-th interpolation point within distance ADELT of XOPT. It is
!     reached via label 90. Otherwise, there is a branch to label 60 for
!     another trust region iteration, unless the calculations with the
!     current RHO are complete.
!
      IF (KNEW .GT. 0) THEN
          DIST=DSQRT(DISTSQ)
          IF (NTRITS .EQ. -1) THEN
              DELTA=DMIN1(TENTH*DELTA,HALF*DIST)
              IF (DELTA .LE. 1.5D0*RHO) DELTA=RHO
          END IF
          NTRITS=0
          ADELT=DMAX1(DMIN1(TENTH*DIST,DELTA),RHO)
          DSQ=ADELT*ADELT
          GOTO 90
      END IF
      IF (NTRITS .EQ. -1) GOTO 680
      IF (RATIO .GT. ZERO) GOTO 60
      IF (DMAX1(DELTA,DNORM) .GT. RHO) GOTO 60
!
!     The calculations with the current value of RHO are complete. Pick the
!       next values of RHO and DELTA.
!
  680 IF (RHO .GT. RHOEND) THEN
          DELTA=HALF*RHO
          RATIO=RHO/RHOEND
          IF (RATIO .LE. 16.0D0) THEN
              RHO=RHOEND
          ELSE IF (RATIO .LE. 250.0D0) THEN
              RHO=DSQRT(RATIO)*RHOEND
          ELSE
              RHO=TENTH*RHO
          END IF
          DELTA=DMAX1(DELTA,RHO)
          IF (IPRINT .GE. 2) THEN
              IF (IPRINT .GE. 3) PRINT 690
  690         FORMAT (5X)
              PRINT 700, RHO,NF
  700         FORMAT (/4X,'New RHO =',1PD11.4,5X,'Number of', &
                ' function values =',I6)
              PRINT 710, FVAL(KOPT),(XBASE(I)+XOPT(I),I=1,N)
  710         FORMAT (4X,'Least value of F =',1PD23.15,9X, &
                'The corresponding X is:'/(2X,5D15.6))
          END IF
          NTRITS=0
          NFSAV=NF
          GOTO 60
      END IF
!
!     Return from the calculation, after another Newton-Raphson step, if
!       it is too short to have been tried before.
!
      IF (NTRITS .EQ. -1) GOTO 360
  720 IF (FVAL(KOPT) .LE. FSAVE) THEN
          DO 730 I=1,N
          X(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XOPT(I)),XU(I))
          IF (XOPT(I) .EQ. SL(I)) X(I)=XL(I)
          IF (XOPT(I) .EQ. SU(I)) X(I)=XU(I)
  730     CONTINUE
          F=FVAL(KOPT)
      END IF
      IF (IPRINT .GE. 1) THEN
          PRINT 740, NF
  740     FORMAT (/4X,'At the return from BOBYQA',5X, &
            'Number of function values =',I6)
          PRINT 710, F,(X(I),I=1,N)
!          write(99,700) rho, nf ! printa rho e numero de avalicao de fun
!          write(99,710) F,(X(I),I=1,N) ! printa residuo e valores otimos
!          write(99,"(i6,100(e15.7,2X))") nf,F,(X(I),I=1,N)
!            nst(it) = nf
!            fst(it) = F
!            do i=1,n
!            xst(it,i) = x(i)
!            end do !i
      END IF
!
      nst(ii) = nf
      fst(ii) = F
      do i=1,n
      xst(ii,i) = x(i)
      end do !i
!
      RETURN
      END

!--------------------------------------------------------------------------

      SUBROUTINE PRELIM (N,NPT,X,XL,XU,RHOBEG,IPRINT,MAXFUN,XBASE, &
        XPT,FVAL,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,KOPT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),GOPT(*), &
        HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*)
!
!     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the
!       same as the corresponding arguments in SUBROUTINE BOBYQA.
!     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU
!       are the same as the corresponding arguments in BOBYQB, the elements
!       of SL and SU being set in BOBYQA.
!     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but
!       it is set by PRELIM to the gradient of the quadratic model at XBASE.
!       If XOPT is nonzero, BOBYQB will change it to its usual value later.
!     NF is maintaned as the number of calls of CALFUN so far.
!     KOPT will be such that the least calculated value of F so far is at
!       the point XPT(KOPT,.)+XBASE in the space of the variables.
!
!     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ,
!     BMAT and ZMAT for the first iteration, and it maintains the values of
!     NF and KOPT. The vector X is also changed by PRELIM.
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      TWO=2.0D0
      ZERO=0.0D0
      RHOSQ=RHOBEG*RHOBEG
      RECIP=ONE/RHOSQ
      NP=N+1
!
!     Set XBASE to the initial vector of variables, and set the initial
!     elements of XPT, BMAT, HQ, PQ and ZMAT to zero.
!
      DO 20 J=1,N
      XBASE(J)=X(J)
      DO 10 K=1,NPT
   10 XPT(K,J)=ZERO
      DO 20 I=1,NDIM
   20 BMAT(I,J)=ZERO
      DO 30 IH=1,(N*NP)/2
   30 HQ(IH)=ZERO
      DO 40 K=1,NPT
      PQ(K)=ZERO
      DO 40 J=1,NPT-NP
   40 ZMAT(K,J)=ZERO
!
!     Begin the initialization procedure. NF becomes one more than the number
!     of function values so far. The coordinates of the displacement of the
!     next initial interpolation point from XBASE are set in XPT(NF+1,.).
!
      NF=0
   50 NFM=NF
      NFX=NF-N
      NF=NF+1
      IF (NFM .LE. 2*N) THEN
          IF (NFM .GE. 1 .AND. NFM .LE. N) THEN
              STEPA=RHOBEG
              IF (SU(NFM) .EQ. ZERO) STEPA=-STEPA
              XPT(NF,NFM)=STEPA
          ELSE IF (NFM .GT. N) THEN
              STEPA=XPT(NF-N,NFX)
              STEPB=-RHOBEG
              IF (SL(NFX) .EQ. ZERO) STEPB=DMIN1(TWO*RHOBEG,SU(NFX))
              IF (SU(NFX) .EQ. ZERO) STEPB=DMAX1(-TWO*RHOBEG,SL(NFX))
              XPT(NF,NFX)=STEPB
          END IF
      ELSE
          ITEMP=(NFM-NP)/N
          JPT=NFM-ITEMP*N-N
          IPT=JPT+ITEMP
          IF (IPT .GT. N) THEN
              ITEMP=JPT
              JPT=IPT-N
              IPT=ITEMP
          END IF
          XPT(NF,IPT)=XPT(IPT+1,IPT)
          XPT(NF,JPT)=XPT(JPT+1,JPT)
      END IF
!
!     Calculate the next value of F. The least function value so far and
!     its index are required.
!
      DO 60 J=1,N
      X(J)=DMIN1(DMAX1(XL(J),XBASE(J)+XPT(NF,J)),XU(J))
      IF (XPT(NF,J) .EQ. SL(J)) X(J)=XL(J)
      IF (XPT(NF,J) .EQ. SU(J)) X(J)=XU(J)
   60 CONTINUE
      CALL CALFUN (N,X,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 70, NF,F,(X(I),I=1,N)
   70      FORMAT (/4X,'Function number',I6,'    F =',1PD18.10, &
             '    The corresponding X is:'/(2X,5D15.6))
      END IF
      FVAL(NF)=F
      IF (NF .EQ. 1) THEN
          FBEG=F
          KOPT=1
      ELSE IF (F .LT. FVAL(KOPT)) THEN
          KOPT=NF
      END IF
!
!     Set the nonzero initial elements of BMAT and the quadratic model in the
!     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions
!     of the NF-th and (NF-N)-th interpolation points may be switched, in
!     order that the function value at the first of them contributes to the
!     off-diagonal second derivative terms of the initial quadratic model.
!
      IF (NF .LE. 2*N+1) THEN
          IF (NF .GE. 2 .AND. NF .LE. N+1) THEN
              GOPT(NFM)=(F-FBEG)/STEPA
              IF (NPT .LT. NF+N) THEN
                  BMAT(1,NFM)=-ONE/STEPA
                  BMAT(NF,NFM)=ONE/STEPA
                  BMAT(NPT+NFM,NFM)=-HALF*RHOSQ
              END IF
          ELSE IF (NF .GE. N+2) THEN
              IH=(NFX*(NFX+1))/2
              TEMP=(F-FBEG)/STEPB
              DIFF=STEPB-STEPA
              HQ(IH)=TWO*(TEMP-GOPT(NFX))/DIFF
              GOPT(NFX)=(GOPT(NFX)*STEPB-TEMP*STEPA)/DIFF
              IF (STEPA*STEPB .LT. ZERO) THEN
                  IF (F .LT. FVAL(NF-N)) THEN
                      FVAL(NF)=FVAL(NF-N)
                      FVAL(NF-N)=F
                      IF (KOPT .EQ. NF) KOPT=NF-N
                      XPT(NF-N,NFX)=STEPB
                      XPT(NF,NFX)=STEPA
                  END IF
              END IF
              BMAT(1,NFX)=-(STEPA+STEPB)/(STEPA*STEPB)
              BMAT(NF,NFX)=-HALF/XPT(NF-N,NFX)
              BMAT(NF-N,NFX)=-BMAT(1,NFX)-BMAT(NF,NFX)
              ZMAT(1,NFX)=DSQRT(TWO)/(STEPA*STEPB)
              ZMAT(NF,NFX)=DSQRT(HALF)/RHOSQ
              ZMAT(NF-N,NFX)=-ZMAT(1,NFX)-ZMAT(NF,NFX)
          END IF
!
!     Set the off-diagonal second derivatives of the Lagrange functions and
!     the initial quadratic model.
!
      ELSE
          IH=(IPT*(IPT-1))/2+JPT
          ZMAT(1,NFX)=RECIP
          ZMAT(NF,NFX)=RECIP
          ZMAT(IPT+1,NFX)=-RECIP
          ZMAT(JPT+1,NFX)=-RECIP
          TEMP=XPT(NF,IPT)*XPT(NF,JPT)
          HQ(IH)=(FBEG-FVAL(IPT+1)-FVAL(JPT+1)+F)/TEMP
      END IF
      IF (NF .LT. NPT .AND. NF .LT. MAXFUN) GOTO 50
      RETURN
      END

!--------------------------------------------------------------------------

      SUBROUTINE RESCUE (N,NPT,XL,XU,IPRINT,MAXFUN,XBASE,XPT, &
        FVAL,XOPT,GOPT,HQ,PQ,BMAT,ZMAT,NDIM,SL,SU,NF,DELTA, &
        KOPT,VLAG,PTSAUX,PTSID,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XL(*),XU(*),XBASE(*),XPT(NPT,*),FVAL(*),XOPT(*), &
        GOPT(*),HQ(*),PQ(*),BMAT(NDIM,*),ZMAT(NPT,*),SL(*),SU(*), &
        VLAG(*),PTSAUX(2,*),PTSID(*),W(*)
!
!     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT,
!       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as
!       the corresponding arguments of BOBYQB on the entry to RESCUE.
!     NF is maintained as the number of calls of CALFUN so far, except that
!       NF is set to -1 if the value of MAXFUN prevents further progress.
!     KOPT is maintained so that FVAL(KOPT) is the least calculated function
!       value. Its correct value must be given on entry. It is updated if a
!       new least function value is found, but the corresponding changes to
!       XOPT and GOPT have to be made later by the calling program.
!     DELTA is the current trust region radius.
!     VLAG is a working space vector that will be used for the values of the
!       provisional Lagrange functions at each of the interpolation points.
!       They are part of a product that requires VLAG to be of length NDIM.
!     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and
!       PTSAUX(2,J) specify the two positions of provisional interpolation
!       points when a nonzero step is taken along e_J (the J-th coordinate
!       direction) through XBASE+XOPT, as specified below. Usually these
!       steps have length DELTA, but other lengths are chosen if necessary
!       in order to satisfy the given bounds on the variables.
!     PTSID is also a working space array. It has NPT components that denote
!       provisional new positions of the original interpolation points, in
!       case changes are needed to restore the linear independence of the
!       interpolation conditions. The K-th point is a candidate for change
!       if and only if PTSID(K) is nonzero. In this case let p and q be the
!       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p
!       and q are both positive, the step from XBASE+XOPT to the new K-th
!       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise
!       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or
!       p=0, respectively.
!     The first NDIM+NPT elements of the array W are used for working space.
!     The final elements of BMAT and ZMAT are set in a well-conditioned way
!       to the values that are appropriate for the new interpolation points.
!     The elements of GOPT, HQ and PQ are also revised to the values that are
!       appropriate to the final quadratic model.
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      ZERO=0.0D0
      NP=N+1
      SFRAC=HALF/DFLOAT(NP)
      NPTM=NPT-NP
!
!     Shift the interpolation points so that XOPT becomes the origin, and set
!     the elements of ZMAT to zero. The value of SUMPQ is required in the
!     updating of HQ below. The squares of the distances from XOPT to the
!     other interpolation points are set at the end of W. Increments of WINC
!     may be added later to these squares to balance the consideration of
!     the choice of point that is going to become current.
!
      SUMPQ=ZERO
      WINC=ZERO
      DO 20 K=1,NPT
      DISTSQ=ZERO
      DO 10 J=1,N
      XPT(K,J)=XPT(K,J)-XOPT(J)
   10 DISTSQ=DISTSQ+XPT(K,J)**2
      SUMPQ=SUMPQ+PQ(K)
      W(NDIM+K)=DISTSQ
      WINC=DMAX1(WINC,DISTSQ)
      DO 20 J=1,NPTM
   20 ZMAT(K,J)=ZERO
!
!     Update HQ so that HQ and PQ define the second derivatives of the model
!     after XBASE has been shifted to the trust region centre.
!
      IH=0
      DO 40 J=1,N
      W(J)=HALF*SUMPQ*XOPT(J)
      DO 30 K=1,NPT
   30 W(J)=W(J)+PQ(K)*XPT(K,J)
      DO 40 I=1,J
      IH=IH+1
   40 HQ(IH)=HQ(IH)+W(I)*XOPT(J)+W(J)*XOPT(I)
!
!     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and
!     also set the elements of PTSAUX.
!
      DO 50 J=1,N
      XBASE(J)=XBASE(J)+XOPT(J)
      SL(J)=SL(J)-XOPT(J)
      SU(J)=SU(J)-XOPT(J)
      XOPT(J)=ZERO
      PTSAUX(1,J)=DMIN1(DELTA,SU(J))
      PTSAUX(2,J)=DMAX1(-DELTA,SL(J))
      IF (PTSAUX(1,J)+PTSAUX(2,J) .LT. ZERO) THEN
          TEMP=PTSAUX(1,J)
          PTSAUX(1,J)=PTSAUX(2,J)
          PTSAUX(2,J)=TEMP
      END IF
      IF (DABS(PTSAUX(2,J)) .LT. HALF*DABS(PTSAUX(1,J))) THEN
          PTSAUX(2,J)=HALF*PTSAUX(1,J)
      END IF
      DO 50 I=1,NDIM
   50 BMAT(I,J)=ZERO
      FBASE=FVAL(KOPT)
!
!     Set the identifiers of the artificial interpolation points that are
!     along a coordinate direction from XOPT, and set the corresponding
!     nonzero elements of BMAT and ZMAT.
!
      PTSID(1)=SFRAC
      DO 60 J=1,N
      JP=J+1
      JPN=JP+N
      PTSID(JP)=DFLOAT(J)+SFRAC
      IF (JPN .LE. NPT) THEN
          PTSID(JPN)=DFLOAT(J)/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,J)-PTSAUX(2,J))
          BMAT(JP,J)=-TEMP+ONE/PTSAUX(1,J)
          BMAT(JPN,J)=TEMP+ONE/PTSAUX(2,J)
          BMAT(1,J)=-BMAT(JP,J)-BMAT(JPN,J)
          ZMAT(1,J)=DSQRT(2.0D0)/DABS(PTSAUX(1,J)*PTSAUX(2,J))
          ZMAT(JP,J)=ZMAT(1,J)*PTSAUX(2,J)*TEMP
          ZMAT(JPN,J)=-ZMAT(1,J)*PTSAUX(1,J)*TEMP
      ELSE
          BMAT(1,J)=-ONE/PTSAUX(1,J)
          BMAT(JP,J)=ONE/PTSAUX(1,J)
          BMAT(J+NPT,J)=-HALF*PTSAUX(1,J)**2
      END IF
   60 CONTINUE
!
!     Set any remaining identifiers with their nonzero elements of ZMAT.
!
      IF (NPT .GE. N+NP) THEN
          DO 70 K=2*NP,NPT
          IW=(DFLOAT(K-NP)-HALF)/DFLOAT(N)
          IP=K-NP-IW*N
          IQ=IP+IW
          IF (IQ .GT. N) IQ=IQ-N
          PTSID(K)=DFLOAT(IP)+DFLOAT(IQ)/DFLOAT(NP)+SFRAC
          TEMP=ONE/(PTSAUX(1,IP)*PTSAUX(1,IQ))
          ZMAT(1,K-NP)=TEMP
          ZMAT(IP+1,K-NP)=-TEMP
          ZMAT(IQ+1,K-NP)=-TEMP
   70     ZMAT(K,K-NP)=TEMP
      END IF
      NREM=NPT
      KOLD=1
      KNEW=KOPT
!
!     Reorder the provisional points in the way that exchanges PTSID(KOLD)
!     with PTSID(KNEW).
!
   80 DO 90 J=1,N
      TEMP=BMAT(KOLD,J)
      BMAT(KOLD,J)=BMAT(KNEW,J)
   90 BMAT(KNEW,J)=TEMP
      DO 100 J=1,NPTM
      TEMP=ZMAT(KOLD,J)
      ZMAT(KOLD,J)=ZMAT(KNEW,J)
  100 ZMAT(KNEW,J)=TEMP
      PTSID(KOLD)=PTSID(KNEW)
      PTSID(KNEW)=ZERO
      W(NDIM+KNEW)=ZERO
      NREM=NREM-1
      IF (KNEW .NE. KOPT) THEN
          TEMP=VLAG(KOLD)
          VLAG(KOLD)=VLAG(KNEW)
          VLAG(KNEW)=TEMP
!
!     Update the BMAT and ZMAT matrices so that the status of the KNEW-th
!     interpolation point can be changed from provisional to original. The
!     branch to label 350 occurs if all the original points are reinstated.
!     The nonnegative values of W(NDIM+K) are required in the search below.
!
          CALL UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM,KNEW,W)
          IF (NREM .EQ. 0) GOTO 350
          DO 110 K=1,NPT
  110     W(NDIM+K)=DABS(W(NDIM+K))
      END IF
!
!     Pick the index KNEW of an original interpolation point that has not
!     yet replaced one of the provisional interpolation points, giving
!     attention to the closeness to XOPT and to previous tries with KNEW.
!
  120 DSQMIN=ZERO
      DO 130 K=1,NPT
      IF (W(NDIM+K) .GT. ZERO) THEN
          IF (DSQMIN .EQ. ZERO .OR. W(NDIM+K) .LT. DSQMIN) THEN
              KNEW=K
              DSQMIN=W(NDIM+K)
          END IF
      END IF
  130 CONTINUE
      IF (DSQMIN .EQ. ZERO) GOTO 260
!
!     Form the W-vector of the chosen original interpolation point.
!
      DO 140 J=1,N
  140 W(NPT+J)=XPT(KNEW,J)
      DO 160 K=1,NPT
      SUM=ZERO
      IF (K .EQ. KOPT) THEN
          CONTINUE
      ELSE IF (PTSID(K) .EQ. ZERO) THEN
          DO 150 J=1,N
  150     SUM=SUM+W(NPT+J)*XPT(K,J)
      ELSE
          IP=PTSID(K)
          IF (IP .GT. 0) SUM=W(NPT+IP)*PTSAUX(1,IP)
          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
          IF (IQ .GT. 0) THEN
              IW=1
              IF (IP .EQ. 0) IW=2
              SUM=SUM+W(NPT+IQ)*PTSAUX(IW,IQ)
          END IF
      END IF
  160 W(K)=HALF*SUM*SUM
!
!     Calculate VLAG and BETA for the required updating of the H matrix if
!     XPT(KNEW,.) is reinstated in the set of interpolation points.
!
      DO 180 K=1,NPT
      SUM=ZERO
      DO 170 J=1,N
  170 SUM=SUM+BMAT(K,J)*W(NPT+J)
  180 VLAG(K)=SUM
      BETA=ZERO
      DO 200 J=1,NPTM
      SUM=ZERO
      DO 190 K=1,NPT
  190 SUM=SUM+ZMAT(K,J)*W(K)
      BETA=BETA-SUM*SUM
      DO 200 K=1,NPT
  200 VLAG(K)=VLAG(K)+SUM*ZMAT(K,J)
      BSUM=ZERO
      DISTSQ=ZERO
      DO 230 J=1,N
      SUM=ZERO
      DO 210 K=1,NPT
  210 SUM=SUM+BMAT(K,J)*W(K)
      JP=J+NPT
      BSUM=BSUM+SUM*W(JP)
      DO 220 IP=NPT+1,NDIM
  220 SUM=SUM+BMAT(IP,J)*W(IP)
      BSUM=BSUM+SUM*W(JP)
      VLAG(JP)=SUM
  230 DISTSQ=DISTSQ+XPT(KNEW,J)**2
      BETA=HALF*DISTSQ*DISTSQ+BETA-BSUM
      VLAG(KOPT)=VLAG(KOPT)+ONE
!
!     KOLD is set to the index of the provisional interpolation point that is
!     going to be deleted to make way for the KNEW-th original interpolation
!     point. The choice of KOLD is governed by the avoidance of a small value
!     of the denominator in the updating calculation of UPDATE.
!
      DENOM=ZERO
      VLMXSQ=ZERO
      DO 250 K=1,NPT
      IF (PTSID(K) .NE. ZERO) THEN
          HDIAG=ZERO
          DO 240 J=1,NPTM
  240     HDIAG=HDIAG+ZMAT(K,J)**2
          DEN=BETA*HDIAG+VLAG(K)**2
          IF (DEN .GT. DENOM) THEN
              KOLD=K
              DENOM=DEN
          END IF
      END IF
  250 VLMXSQ=DMAX1(VLMXSQ,VLAG(K)**2)
      IF (DENOM .LE. 1.0D-2*VLMXSQ) THEN
          W(NDIM+KNEW)=-W(NDIM+KNEW)-WINC
          GOTO 120
      END IF
      GOTO 80
!
!     When label 260 is reached, all the final positions of the interpolation
!     points have been chosen although any changes have not been included yet
!     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart
!     from the shift of XBASE, the updating of the quadratic model remains to
!     be done. The following cycle through the new interpolation points begins
!     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero,
!     except that a RETURN occurs if MAXFUN prohibits another value of F.
!
  260 DO 340 KPT=1,NPT
      IF (PTSID(KPT) .EQ. ZERO) GOTO 340
      IF (NF .GE. MAXFUN) THEN
          NF=-1
          GOTO 350
      END IF
      IH=0
      DO 270 J=1,N
      W(J)=XPT(KPT,J)
      XPT(KPT,J)=ZERO
      TEMP=PQ(KPT)*W(J)
      DO 270 I=1,J
      IH=IH+1
  270 HQ(IH)=HQ(IH)+TEMP*W(I)
      PQ(KPT)=ZERO
      IP=PTSID(KPT)
      IQ=DFLOAT(NP)*PTSID(KPT)-DFLOAT(IP*NP)
      IF (IP .GT. 0) THEN
          XP=PTSAUX(1,IP)
          XPT(KPT,IP)=XP
      END IF
      IF (IQ .GT. 0) THEN
          XQ=PTSAUX(1,IQ)
          IF (IP .EQ. 0) XQ=PTSAUX(2,IQ)
          XPT(KPT,IQ)=XQ
      END IF
!
!     Set VQUAD to the value of the current model at the new point.
!
      VQUAD=FBASE
      IF (IP .GT. 0) THEN
          IHP=(IP+IP*IP)/2
          VQUAD=VQUAD+XP*(GOPT(IP)+HALF*XP*HQ(IHP))
      END IF
      IF (IQ .GT. 0) THEN
          IHQ=(IQ+IQ*IQ)/2
          VQUAD=VQUAD+XQ*(GOPT(IQ)+HALF*XQ*HQ(IHQ))
          IF (IP .GT. 0) THEN
              IW=MAX0(IHP,IHQ)-IABS(IP-IQ)
              VQUAD=VQUAD+XP*XQ*HQ(IW)
          END IF
      END IF
      DO 280 K=1,NPT
      TEMP=ZERO
      IF (IP .GT. 0) TEMP=TEMP+XP*XPT(K,IP)
      IF (IQ .GT. 0) TEMP=TEMP+XQ*XPT(K,IQ)
  280 VQUAD=VQUAD+HALF*PQ(K)*TEMP*TEMP
!
!     Calculate F at the new interpolation point, and set DIFF to the factor
!     that is going to multiply the KPT-th Lagrange function when the model
!     is updated to provide interpolation to the new function value.
!
      DO 290 I=1,N
      W(I)=DMIN1(DMAX1(XL(I),XBASE(I)+XPT(KPT,I)),XU(I))
      IF (XPT(KPT,I) .EQ. SL(I)) W(I)=XL(I)
      IF (XPT(KPT,I) .EQ. SU(I)) W(I)=XU(I)
  290 CONTINUE
      NF=NF+1
      CALL CALFUN (N,W,F)
      IF (IPRINT .EQ. 3) THEN
          PRINT 300, NF,F,(W(I),I=1,N)
  300     FORMAT (/4X,'Function number',I6,'    F =',1PD18.10, &
            '    The corresponding X is:'/(2X,5D15.6))
      END IF
      FVAL(KPT)=F
      IF (F .LT. FVAL(KOPT)) KOPT=KPT
      DIFF=F-VQUAD
!
!     Update the quadratic model. The RETURN from the subroutine occurs when
!     all the new interpolation points are included in the model.
!
      DO 310 I=1,N
  310 GOPT(I)=GOPT(I)+DIFF*BMAT(KPT,I)
      DO 330 K=1,NPT
      SUM=ZERO
      DO 320 J=1,NPTM
  320 SUM=SUM+ZMAT(K,J)*ZMAT(KPT,J)
      TEMP=DIFF*SUM
      IF (PTSID(K) .EQ. ZERO) THEN
          PQ(K)=PQ(K)+TEMP
      ELSE
          IP=PTSID(K)
          IQ=DFLOAT(NP)*PTSID(K)-DFLOAT(IP*NP)
          IHQ=(IQ*IQ+IQ)/2
          IF (IP .EQ. 0) THEN
              HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(2,IQ)**2
          ELSE
              IHP=(IP*IP+IP)/2
              HQ(IHP)=HQ(IHP)+TEMP*PTSAUX(1,IP)**2
              IF (IQ .GT. 0) THEN
                  HQ(IHQ)=HQ(IHQ)+TEMP*PTSAUX(1,IQ)**2
                  IW=MAX0(IHP,IHQ)-IABS(IQ-IP)
                  HQ(IW)=HQ(IW)+TEMP*PTSAUX(1,IP)*PTSAUX(1,IQ)
              END IF
          END IF
      END IF
  330 CONTINUE
      PTSID(KPT)=ZERO
  340 CONTINUE
  350 RETURN
      END

!--------------------------------------------------------------------------

      SUBROUTINE TRSBOX (N,NPT,XPT,XOPT,GOPT,HQ,PQ,SL,SU,DELTA, &
        XNEW,D,GNEW,XBDI,S,HS,HRED,DSQ,CRVMIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XPT(NPT,*),XOPT(*),GOPT(*),HQ(*),PQ(*),SL(*),SU(*), &
        XNEW(*),D(*),GNEW(*),XBDI(*),S(*),HS(*),HRED(*)
!
!     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same
!       meanings as the corresponding arguments of BOBYQB.
!     DELTA is the trust region radius for the present calculation, which
!       seeks a small value of the quadratic model within distance DELTA of
!       XOPT subject to the bounds on the variables.
!     XNEW will be set to a new vector of variables that is approximately
!       the one that minimizes the quadratic model within the trust region
!       subject to the SL and SU constraints on the variables. It satisfies
!       as equations the bounds that become active during the calculation.
!     D is the calculated trial step from XOPT, generated iteratively from an
!       initial value of zero. Thus XNEW is XOPT+D after the final iteration.
!     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated
!       when D is updated.
!     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is
!       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the
!       I-th variable has become fixed at a bound, the bound being SL(I) or
!       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This
!       information is accumulated during the construction of XNEW.
!     The arrays S, HS and HRED are also used for working space. They hold the
!       current search direction, and the changes in the gradient of Q along S
!       and the reduced D, respectively, where the reduced D is the same as D,
!       except that the components of the fixed variables are zero.
!     DSQ will be set to the square of the length of XNEW-XOPT.
!     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise
!       it is set to the least curvature of H that occurs in the conjugate
!       gradient searches that are not restricted by any constraints. The
!       value CRVMIN=-1.0D0 is set, however, if all of these searches are
!       constrained.
!
!     A version of the truncated conjugate gradient is applied. If a line
!     search is restricted by a constraint, then the procedure is restarted,
!     the values of the variables that are at their bounds being fixed. If
!     the trust region boundary is reached, then further changes may be made
!     to D, each one being in the two dimensional space that is spanned
!     by the current D and the gradient of Q at XOPT+D, staying on the trust
!     region boundary. Termination occurs when the reduction in Q seems to
!     be close to the greatest reduction that can be achieved.
!
!     Set some constants.
!
      HALF=0.5D0
      ONE=1.0D0
      ONEMIN=-1.0D0
      ZERO=0.0D0
!
!     The sign of GOPT(I) gives the sign of the change to the I-th variable
!     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether
!     or not to fix the I-th variable at one of its bounds initially, with
!     NACT being set to the number of fixed variables. D and GNEW are also
!     set for the first iteration. DELSQ is the upper bound on the sum of
!     squares of the free variables. QRED is the reduction in Q so far.
!
      ITERC=0
      NACT=0
      SQSTP=ZERO
      DO 10 I=1,N
      XBDI(I)=ZERO
      IF (XOPT(I) .LE. SL(I)) THEN
          IF (GOPT(I) .GE. ZERO) XBDI(I)=ONEMIN
      ELSE IF (XOPT(I) .GE. SU(I)) THEN
          IF (GOPT(I) .LE. ZERO) XBDI(I)=ONE
      END IF
      IF (XBDI(I) .NE. ZERO) NACT=NACT+1
      D(I)=ZERO
   10 GNEW(I)=GOPT(I)
      DELSQ=DELTA*DELTA
      QRED=ZERO
      CRVMIN=ONEMIN
!
!     Set the next search direction of the conjugate gradient method. It is
!     the steepest descent direction initially and when the iterations are
!     restarted because a variable has just been fixed by a bound, and of
!     course the components of the fixed variables are zero. ITERMAX is an
!     upper bound on the indices of the conjugate gradient iterations.
!
   20 BETA=ZERO
   30 STEPSQ=ZERO
      DO 40 I=1,N
      IF (XBDI(I) .NE. ZERO) THEN
          S(I)=ZERO
      ELSE IF (BETA .EQ. ZERO) THEN
          S(I)=-GNEW(I)
      ELSE
          S(I)=BETA*S(I)-GNEW(I)
      END IF
   40 STEPSQ=STEPSQ+S(I)**2
      IF (STEPSQ .EQ. ZERO) GOTO 190
      IF (BETA .EQ. ZERO) THEN
          GREDSQ=STEPSQ
          ITERMAX=ITERC+N-NACT
      END IF
      IF (GREDSQ*DELSQ .LE. 1.0D-4*QRED*QRED) GO TO 190
!
!     Multiply the search direction by the second derivative matrix of Q and
!     calculate some scalars for the choice of steplength. Then set BLEN to
!     the length of the the step to the trust region boundary and STPLEN to
!     the steplength, ignoring the simple bounds.
!
      GOTO 210
   50 RESID=DELSQ
      DS=ZERO
      SHS=ZERO
      DO 60 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          RESID=RESID-D(I)**2
          DS=DS+S(I)*D(I)
          SHS=SHS+S(I)*HS(I)
      END IF
   60 CONTINUE
      IF (RESID .LE. ZERO) GOTO 90
      TEMP=DSQRT(STEPSQ*RESID+DS*DS)
      IF (DS .LT. ZERO) THEN
          BLEN=(TEMP-DS)/STEPSQ
      ELSE
          BLEN=RESID/(TEMP+DS)
      END IF
      STPLEN=BLEN
      IF (SHS .GT. ZERO) THEN
          STPLEN=DMIN1(BLEN,GREDSQ/SHS)
      END IF

!
!     Reduce STPLEN if necessary in order to preserve the simple bounds,
!     letting IACT be the index of the new constrained variable.
!
      IACT=0
      DO 70 I=1,N
      IF (S(I) .NE. ZERO) THEN
          XSUM=XOPT(I)+D(I)
          IF (S(I) .GT. ZERO) THEN
              TEMP=(SU(I)-XSUM)/S(I)
          ELSE
              TEMP=(SL(I)-XSUM)/S(I)
          END IF
          IF (TEMP .LT. STPLEN) THEN
              STPLEN=TEMP
              IACT=I
          END IF
      END IF
   70 CONTINUE
!
!     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q.
!
      SDEC=ZERO
      IF (STPLEN .GT. ZERO) THEN
          ITERC=ITERC+1
          TEMP=SHS/STEPSQ
          IF (IACT .EQ. 0 .AND. TEMP .GT. ZERO) THEN
              CRVMIN=DMIN1(CRVMIN,TEMP)
              IF (CRVMIN .EQ. ONEMIN) CRVMIN=TEMP
          END IF
          GGSAV=GREDSQ
          GREDSQ=ZERO
          DO 80 I=1,N
          GNEW(I)=GNEW(I)+STPLEN*HS(I)
          IF (XBDI(I) .EQ. ZERO) GREDSQ=GREDSQ+GNEW(I)**2
   80     D(I)=D(I)+STPLEN*S(I)
          SDEC=DMAX1(STPLEN*(GGSAV-HALF*STPLEN*SHS),ZERO)
          QRED=QRED+SDEC
      END IF
!
!     Restart the conjugate gradient method if it has hit a new bound.
!
      IF (IACT .GT. 0) THEN
          NACT=NACT+1
          XBDI(IACT)=ONE
          IF (S(IACT) .LT. ZERO) XBDI(IACT)=ONEMIN
          DELSQ=DELSQ-D(IACT)**2
          IF (DELSQ .LE. ZERO) GOTO 90
          GOTO 20
      END IF
!
!     If STPLEN is less than BLEN, then either apply another conjugate
!     gradient iteration or RETURN.
!
      IF (STPLEN .LT. BLEN) THEN
          IF (ITERC .EQ. ITERMAX) GOTO 190
          IF (SDEC .LE. 0.01D0*QRED) GOTO 190
          BETA=GREDSQ/GGSAV
          GOTO 30
      END IF
   90 CRVMIN=ZERO
!
!     Prepare for the alternative iteration by calculating some scalars and
!     by multiplying the reduced D by the second derivative matrix of Q.
!
  100 IF (NACT .GE. N-1) GOTO 190
      DREDSQ=ZERO
      DREDG=ZERO
      GREDSQ=ZERO
      DO 110 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          DREDSQ=DREDSQ+D(I)**2
          DREDG=DREDG+D(I)*GNEW(I)
          GREDSQ=GREDSQ+GNEW(I)**2
          S(I)=D(I)
      ELSE
          S(I)=ZERO
      END IF
  110 CONTINUE
      ITCSAV=ITERC
      GOTO 210
!
!     Let the search direction S be a linear combination of the reduced D
!     and the reduced G that is orthogonal to the reduced D.
!
  120 ITERC=ITERC+1
      TEMP=GREDSQ*DREDSQ-DREDG*DREDG
      IF (TEMP .LE. 1.0D-4*QRED*QRED) GOTO 190
      TEMP=DSQRT(TEMP)
      DO 130 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          S(I)=(DREDG*D(I)-DREDSQ*GNEW(I))/TEMP
      ELSE
          S(I)=ZERO
      END IF
  130 CONTINUE
      SREDG=-TEMP
!
!     By considering the simple bounds on the variables, calculate an upper
!     bound on the tangent of half the angle of the alternative iteration,
!     namely ANGBD, except that, if already a free variable has reached a
!     bound, there is a branch back to label 100 after fixing that variable.
!
      ANGBD=ONE
      IACT=0
      DO 140 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          TEMPA=XOPT(I)+D(I)-SL(I)
          TEMPB=SU(I)-XOPT(I)-D(I)
          IF (TEMPA .LE. ZERO) THEN
              NACT=NACT+1
              XBDI(I)=ONEMIN
              GOTO 100
          ELSE IF (TEMPB .LE. ZERO) THEN
              NACT=NACT+1
              XBDI(I)=ONE
              GOTO 100
          END IF
          RATIO=ONE
          SSQ=D(I)**2+S(I)**2
          TEMP=SSQ-(XOPT(I)-SL(I))**2
          IF (TEMP .GT. ZERO) THEN
              TEMP=DSQRT(TEMP)-S(I)
              IF (ANGBD*TEMP .GT. TEMPA) THEN
                  ANGBD=TEMPA/TEMP
                  IACT=I
                  XSAV=ONEMIN
              END IF
          END IF
          TEMP=SSQ-(SU(I)-XOPT(I))**2
          IF (TEMP .GT. ZERO) THEN
              TEMP=DSQRT(TEMP)+S(I)
              IF (ANGBD*TEMP .GT. TEMPB) THEN
                  ANGBD=TEMPB/TEMP
                  IACT=I
                  XSAV=ONE
              END IF
          END IF
      END IF
  140 CONTINUE
!
!     Calculate HHD and some curvatures for the alternative iteration.
!
      GOTO 210
  150 SHS=ZERO
      DHS=ZERO
      DHD=ZERO
      DO 160 I=1,N
      IF (XBDI(I) .EQ. ZERO) THEN
          SHS=SHS+S(I)*HS(I)
          DHS=DHS+D(I)*HS(I)
          DHD=DHD+D(I)*HRED(I)
      END IF
  160 CONTINUE
!
!     Seek the greatest reduction in Q for a range of equally spaced values
!     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of
!     the alternative iteration.
!
      REDMAX=ZERO
      ISAV=0
      REDSAV=ZERO
      IU=17.0D0*ANGBD+3.1D0
      DO 170 I=1,IU
      ANGT=ANGBD*DFLOAT(I)/DFLOAT(IU)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      REDNEW=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      IF (REDNEW .GT. REDMAX) THEN
          REDMAX=REDNEW
          ISAV=I
          RDPREV=REDSAV
      ELSE IF (I .EQ. ISAV+1) THEN
          RDNEXT=REDNEW
      END IF
  170 REDSAV=REDNEW
!
!     Return if the reduction is zero. Otherwise, set the sine and cosine
!     of the angle of the alternative iteration, and calculate SDEC.
!
      IF (ISAV .EQ. 0) GOTO 190
      IF (ISAV .LT. IU) THEN
          TEMP=(RDNEXT-RDPREV)/(REDMAX+REDMAX-RDPREV-RDNEXT)
          ANGT=ANGBD*(DFLOAT(ISAV)+HALF*TEMP)/DFLOAT(IU)
      END IF
      CTH=(ONE-ANGT*ANGT)/(ONE+ANGT*ANGT)
      STH=(ANGT+ANGT)/(ONE+ANGT*ANGT)
      TEMP=SHS+ANGT*(ANGT*DHD-DHS-DHS)
      SDEC=STH*(ANGT*DREDG-SREDG-HALF*STH*TEMP)
      IF (SDEC .LE. ZERO) GOTO 190
!
!     Update GNEW, D and HRED. If the angle of the alternative iteration
!     is restricted by a bound on a free variable, that variable is fixed
!     at the bound.
!
      DREDG=ZERO
      GREDSQ=ZERO
      DO 180 I=1,N
      GNEW(I)=GNEW(I)+(CTH-ONE)*HRED(I)+STH*HS(I)
      IF (XBDI(I) .EQ. ZERO) THEN
          D(I)=CTH*D(I)+STH*S(I)
          DREDG=DREDG+D(I)*GNEW(I)
          GREDSQ=GREDSQ+GNEW(I)**2
      END IF
  180 HRED(I)=CTH*HRED(I)+STH*HS(I)
      QRED=QRED+SDEC
      IF (IACT .GT. 0 .AND. ISAV .EQ. IU) THEN
          NACT=NACT+1
          XBDI(IACT)=XSAV
          GOTO 100
      END IF
!
!     If SDEC is sufficiently small, then RETURN after setting XNEW to
!     XOPT+D, giving careful attention to the bounds.
!
      IF (SDEC .GT. 0.01D0*QRED) GOTO 120
  190 DSQ=ZERO
      DO 200 I=1,N
      XNEW(I)=DMAX1(DMIN1(XOPT(I)+D(I),SU(I)),SL(I))
      IF (XBDI(I) .EQ. ONEMIN) XNEW(I)=SL(I)
      IF (XBDI(I) .EQ. ONE) XNEW(I)=SU(I)
      D(I)=XNEW(I)-XOPT(I)
  200 DSQ=DSQ+D(I)**2
      RETURN

!     The following instructions multiply the current S-vector by the second
!     derivative matrix of the quadratic model, putting the product in HS.
!     They are reached from three different parts of the software above and
!     they can be regarded as an external subroutine.
!
  210 IH=0
      DO 220 J=1,N
      HS(J)=ZERO
      DO 220 I=1,J
      IH=IH+1
      IF (I .LT. J) HS(J)=HS(J)+HQ(IH)*S(I)
  220 HS(I)=HS(I)+HQ(IH)*S(J)
      DO 250 K=1,NPT
      IF (PQ(K) .NE. ZERO) THEN
          TEMP=ZERO
          DO 230 J=1,N
  230     TEMP=TEMP+XPT(K,J)*S(J)
          TEMP=TEMP*PQ(K)
          DO 240 I=1,N
  240     HS(I)=HS(I)+TEMP*XPT(K,I)
      END IF
  250 CONTINUE
      IF (CRVMIN .NE. ZERO) GOTO 50
      IF (ITERC .GT. ITCSAV) GOTO 150
      DO 260 I=1,N
  260 HRED(I)=HS(I)
      GOTO 120
      END

!--------------------------------------------------------------------------

      SUBROUTINE UPDATE (N,NPT,BMAT,ZMAT,NDIM,VLAG,BETA,DENOM, &
        KNEW,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION BMAT(NDIM,*),ZMAT(NPT,*),VLAG(*),W(*)
!
!     The arrays BMAT and ZMAT are updated, as required by the new position
!     of the interpolation point that has the index KNEW. The vector VLAG has
!     N+NPT components, set on entry to the first NPT and last N components
!     of the product Hw in equation (4.11) of the Powell (2006) paper on
!     NEWUOA. Further, BETA is set on entry to the value of the parameter
!     with that name, and DENOM is set to the denominator of the updating
!     formula. Elements of ZMAT may be treated as zero if their moduli are
!     at most ZTEST. The first NDIM elements of W are used for working space.
!
!     Set some constants.
!
      ONE=1.0D0
      ZERO=0.0D0
      NPTM=NPT-N-1
      ZTEST=ZERO
      DO 10 K=1,NPT
      DO 10 J=1,NPTM
   10 ZTEST=DMAX1(ZTEST,DABS(ZMAT(K,J)))
      ZTEST=1.0D-20*ZTEST
!
!     Apply the rotations that put zeros in the KNEW-th row of ZMAT.
!
      JL=1
      DO 30 J=2,NPTM
      IF (DABS(ZMAT(KNEW,J)) .GT. ZTEST) THEN
          TEMP=DSQRT(ZMAT(KNEW,1)**2+ZMAT(KNEW,J)**2)
          TEMPA=ZMAT(KNEW,1)/TEMP
          TEMPB=ZMAT(KNEW,J)/TEMP
          DO 20 I=1,NPT
          TEMP=TEMPA*ZMAT(I,1)+TEMPB*ZMAT(I,J)
          ZMAT(I,J)=TEMPA*ZMAT(I,J)-TEMPB*ZMAT(I,1)
   20     ZMAT(I,1)=TEMP
      END IF
      ZMAT(KNEW,J)=ZERO
   30 CONTINUE
!
!     Put the first NPT components of the KNEW-th column of HLAG into W,
!     and calculate the parameters of the updating formula.
!
      DO 40 I=1,NPT
      W(I)=ZMAT(KNEW,1)*ZMAT(I,1)
   40 CONTINUE
      ALPHA=W(KNEW)
      TAU=VLAG(KNEW)
      VLAG(KNEW)=VLAG(KNEW)-ONE
!
!     Complete the updating of ZMAT.
!
      TEMP=DSQRT(DENOM)
      TEMPB=ZMAT(KNEW,1)/TEMP
      TEMPA=TAU/TEMP
      DO 50 I=1,NPT
   50 ZMAT(I,1)=TEMPA*ZMAT(I,1)-TEMPB*VLAG(I)
!
!     Finally, update the matrix BMAT.
!
      DO 60 J=1,N
      JP=NPT+J
      W(JP)=BMAT(KNEW,J)
      TEMPA=(ALPHA*VLAG(JP)-TAU*W(JP))/DENOM
      TEMPB=(-BETA*W(JP)-TAU*VLAG(JP))/DENOM
      DO 60 I=1,JP
      BMAT(I,J)=BMAT(I,J)+TEMPA*VLAG(I)+TEMPB*W(I)
      IF (I .GT. NPT) BMAT(JP,I-NPT)=BMAT(I,J)
   60 CONTINUE
      RETURN
      END

!--------------------------------------------------------------------------

!      SUBROUTINE CALFUN (N,X,F)
!      IMPLICIT REAL*8 (A-H,O-Z)
!      DIMENSION X(*)
!      F=0.0D0
!      DO 10 I=4,N,2
!      DO 10 J=2,I-2,2
!      TEMP=(X(I-1)-X(J-1))**2+(X(I)-X(J))**2
!      TEMP=DMAX1(TEMP,1.0D-6)
!   10 F=F+1.0D0/DSQRT(TEMP)
!      RETURN
!      END

      subroutine calfun(n,x,f)

      use moptim

      implicit real*8 (a-h,o-z)
      dimension x(*)

      f=0.0d0

      call residual(x,f)

      end subroutine

!----------------------------------------------------------------------
!
      end module
!
!----------------------------------------------------------------------
!
      program fitting
!
!     solve mosquitoes problem with optimal control
!
      use mgeometry
      use marrays
      use mcoeficientes
      use merro
      use mgauss
      use mpdes
      use moptim
      use mbobyqa
!
      implicit none
!
      integer :: i,j,indopt
      real(8), dimension(5) :: est
      real(8) :: cflag
!
      integer :: k,l,m,n,jj,kk,ll
      real(8), dimension(100) :: x, xlo, xup
      real(8), dimension(500000) :: w
!
      integer :: iprint, maxfun, npt
      real(8) :: rhobeg, rhoend, bdl, bdu

      real(8), dimension(200,200) :: resout 
!
!     bobyqa parameters
!
      n = 5!6!3!4 ! number of parameters
      npar = n
!
      iprint=3! 1 ou 2 ou 3! output
      maxfun=500000 ! max evaluation number
      rhobeg=1.0d-4!4! ! step size
      rhoend=1.0d-6!6 ! accuracy
      bdl=0.d0+rhoend ! min restriction
      bdu=1.d0!0.1d0 ! max restriction
      npt = 2*n+1!n+6!2*n+1 ! or n+6
!
!     lower and upper bounds
!
      do i=1,n
         xlo(i) = bdl
         xup(i) = bdu
      end do !i
! 
      xlo(1) = 1.d-5 !beta 
      xlo(2) = 1.d-5 !betav
      xlo(3) = 1.d-4 !i01  
! 
      xlo(4) = 1.d-4 !i02
      xlo(5) = 1.d-4 !i03
! 
      xlo(6) = 1.d-7 !betas
      ! xlo(4) = 1.d-5 !betas
!
      xup(1) = 1.d-1!1.d0!0.75d0!5.d-1 !beta 
      xup(2) = 1.d0!0.75d0!5.d-1 !betav
      xup(3) = 3.d-3   !i01  
! 
      xup(4) = 5.d-3!1.d-1 !3.d-3 !i02
      xup(5) = 5.d-3!1.d-1 !3.d-3 !i03
      xup(6) = 1.d-2!1.d-1 !betas

      ! ! ! xlo(1) = 1.d-4 !i02
      ! ! ! xlo(2) = 1.d-4 !i03
      ! ! ! xup(1) = 3.d-2 !i02
      ! ! ! xup(2) = 3.d-2 !i03
! 
!     defining parameters
! 
      call iniparams
!
      allocate(nst(500))
      allocate(fst(500))
      allocate(xst(500,n))
      nst = 0.d0
      fst = 1.d10
      xst = 0.d0
!
      open(unit = 99,file = 'result.dat',status = 'unknown')
!
!     store data for numerical integration
!
      call gausstable
!
!     read mesh information
!
      ! call genmesh
      call readmesh
!
      ndof = nnodes !2*nelem + 1
      nmax = 3*nen*nnodes
!
!     for lapack
!
      ! call calcband
!
!     allocating
!
!     structure
!
      allocate(ilin(nmax))
      allocate(jcol(nmax))      
!
!     state solutions
!
      allocate(suh(ndof,nt+1))
      allocate(inh(ndof,nt+1))
      allocate(reh(ndof,nt+1))
      allocate(suv(ndof,nt+1))
      allocate(inv(ndof,nt+1))
      allocate(shold(ndof,nt+1))
      allocate(ihold(ndof,nt+1))
      allocate(rhold(ndof,nt+1))
      allocate(svold(ndof,nt+1))
      allocate(ivold(ndof,nt+1))
      allocate(intsuh(nt+1,1))
      allocate(intinh(nt+1,1))
      allocate(intreh(nt+1,1))
      allocate(intsuv(nt+1,1))
      allocate(intinv(nt+1,1))
      allocate(intctl(nt+1,1))
      allocate(intj(nt+1,1))
      allocate(intnnc(nt+1,1))
      allocate(intvac(nt+1,1))
      allocate(obsinh(ndof,nt+1))
!
!     adjoint solutions
!
      allocate(l1(ndof,nt+1))
      allocate(l2(ndof,nt+1))
      allocate(l3(ndof,nt+1))
      allocate(l4(ndof,nt+1))
      allocate(l5(ndof,nt+1))
      allocate(l1old(ndof,nt+1))
      allocate(l2old(ndof,nt+1))
      allocate(l3old(ndof,nt+1))
      allocate(l4old(ndof,nt+1))
      allocate(l5old(ndof,nt+1))
!
!     control
!
      allocate(u(ndof,nt+1))
      allocate(uold(ndof,nt+1))
!
!     state system matrices
!
      allocate(mesh(nmax,1))
      allocate(meih(nmax,1))
      allocate(merh(nmax,1))
      allocate(mesv(nmax,1))
      allocate(meiv(nmax,1))
      allocate(mdsh(nmax,1))
      allocate(mdih(nmax,1))
      allocate(mdrh(nmax,1))
      allocate(mdsv(nmax,1))
      allocate(mdiv(nmax,1))
!
!     loading data
!
      call loaddata
!
!...........................
!    solving pdes
! !
! x(1:6)=(/4.131237d-3,9.999437d-1,3.000000d-3,4.809473d-3,5.000000d-3,0.d0/) !2.462972d-2 noise10  5.3013e-3
! x(1:6)=(/4.297419d-3,9.613113d-1,3.000000d-3,4.423685d-3,4.950486d-3,0.d0/) !7.172323d-2 noise30  3.3787e-2
! x(1:6)=(/4.190709d-3,9.891175d-1,2.766257d-3,4.226634d-3,4.359968d-3,0.d0/) !1.564472d-1 noise50  2.6726e-2
x(1:6)=(/4.360170d-3,9.780675d-1,2.033269d-3,3.506517d-3,4.147075d-3,0.d0/) !2.002742d-1 noise70  6.8478E-2
! x(1:6)=(/0.409698D-2,0.100000D+1,0.300000D-2,0.454592D-2,0.492066D-2,0.d0/) !3.051700d-1 noise100 2.2789E-2
! x(1:6)=(/0.412368D-02,0.100000D+01,0.300000D-02,0.500000D-02,0.5d-2,0.d0/) !0.8253 *
! beta  =x(1)
! betav =x(2)
! i01   =x(3)
! i02   =x(4)!0.d0!
! i03   =x(5)!0.d0!
! betas =0.d0 !x(6)

! ! beta=0
! ! betav=0
! ! delta=0
! ! i01=0
! ! i03=0     
call iniparams
call residual(x,resout(1,1))
write(*,*) resout(1,1)
stop
! ! 
!       nt1 = 1
!       nt2 = nt!int(5*7/dt) !0!int(52/dt) !
!       call initial
!       call mosqpdes
! !
!       ! nt1 = nt2+1
!       ! nt2 = nt
!       ! call contpdes
! !
!       call integral
!       call plotsols
!       call plotints
!       stop
!...........................
!
!     days loop
!
!       do expnum=1,140
! !
!       nt1 = 1
!       nt2 = int((expnum-1)/dt) !nt! 
!       call initial
!       call mosqpdes

!       nt1 = nt2+1
!       nt2 = nt
!       call contpdes
! ! 
!       call integral
!       ! call plotsols
!       ! call plotints
! !
!       end do !expnum
!       stop
!
!...........................
!
!     loading LHS matrix
!
!       open(unit = 20, file = 'lhsmatrixctl.dat', status = 'old')

!     exp loop

!       do expnum=1,12!1,nsim
! !
! !     calling lhs parameters
! !
!       ! read(20,*) beta,delta,betav,lambv,kpv,muv,alpha(1),alpha(2),alpha(4),c1,c2,c3
!       ! alpha(3) = alpha(1)
!       ! alpha(5) = alpha(4)
! !
! !     test monotonicity
!       indpar(1:12)=(/1,2,3,4,8,9,10,11,12,13,17,18/)
!       do j=1,41
! ! 
!       paramlhs = params
! !
!       if ((indpar(expnum)==4).or.(indpar(expnum)==18)) then
! !
!       paramlhs(indpar(expnum)) = paramlhs(indpar(expnum))*(dble(j)*0.3d0 + 0.85d0*41.d0-1.15d0)/(41.d0-1.d0)
! ! 
!       else if ((indpar(expnum)==8).or.(indpar(expnum)==10).or.(indpar(expnum)==11)) then
! !
!       paramlhs(indpar(expnum)) = paramlhs(indpar(expnum))*(dble(j)*0.2d0 + 0.9d0*41.d0-1.1d0)/(41.d0-1.d0)
! !
!       else
! !
!       paramlhs(indpar(expnum)) = paramlhs(indpar(expnum))*(dble(j)*0.1d0 + 0.95d0*41.d0-1.05d0)/(41.d0-1.d0)
! !
!       end if 
! !
!       call updtparam(paramlhs)
! !
! !     solving pdes
! !
!       nt1 = 1
!       nt2 = int(5*7/dt) !nt! 
!       call initial
!       call mosqpdes

!       nt1 = nt2+1
!       nt2 = nt
!       call contpdes
! !
!       call integral
!       ! call plotsols
!       ! call plotintslhs
! !
!       end do !j
! ! 
!       end do !expnum
!       stop
!
!...........................
!
      call newseed
!     optimization exp loop
!
      do expnum=1,1!nsim
!
!       do ii=1,1!6
! !
!       x(1)=dble(ii-1)*0.15d0 !0.04d0! initial guess
!       x(2)=dble(ii-1)*0.015d0 !0.04d0! initial guess
!
      ii = 0
      do j=1,50!3!11
! 
      do k=1,1!3!11

      do l=1,1!3!1!

      do jj=1,1!3!11

      do kk=1,1!3!1!
! 
      ii = ii + 1
! 
      ! x = rhobeg
      x(1:6) = xlo(1:6)
      x(3:6) = xup(3:6)
      ! x(6) = xlo(6)
      ! x(2) = 0.75d0
      x(6) = 0.d0
      x(1:6)=(/0.412368D-02,0.100000D+01,0.300000D-02,0.500000D-02,0.5d-2,0.d0/)
!
      ! open(unit = 20, file = 'initguess.dat', status = 'old')
      ! read(20,*) x(1),x(2),x(3),x(4),x(5)
!
      call random_number(x)
      x = xlo + (xup-xlo)*x
      x(6) = 0.d0
      ! xst(ii,1:5) = x(1:5)
!
      ! x(3:5) = 0.d0
      ! x(1:2) = xup(1:2)
      ! x(j) = (xlo(j) + xup(j))/2.d0
      ! x(1:6)=(/0.586552D+00,0.275735D-02,0.300000D-02,0.500000D-02,0.500000D-02,0.d0/)  
      ! x(3:5) = xup(3:5)
      ! x(1:6) = 0.5d0*(xlo(1:6) + xup(1:6))
            ! x(1) = xlo(1)+dble(j-1)*(xup(1)-xlo(1))/(3.d0-1.d0)
            ! x(2) = xlo(2)+dble(k-1)*(xup(2)-xlo(2))/(3.d0-1.d0)
            ! x(3) = xlo(3)+dble(l-1)*(xup(3)-xlo(3))/(3.d0-1.d0)
            ! x(4) = xlo(4)+dble(jj-1)*(xup(4)-xlo(4))/(3.d0-1.d0)
            ! x(5) = xlo(5)+dble(kk-1)*(xup(5)-xlo(5))/(3.d0-1.d0)
      ! betas = 0.d0
      ! x(2) = x(2)*(dble(j)*0.4d0 + 0.8d0*11.d0-1.2d0)/(11.d0-1.d0)!

      ! x(1) = 4.3834d-1 -2.d0*(0.1d-1)+(0.1d-1)*dble(j-1)!
      ! x(2) = 2.8865d-3 -2.d0*(0.1d-3)+(0.1d-3)*dble(k-1)!
      ! x(3) = 1.9806d-3 -2.d0*(0.1d-3)+(0.1d-3)*dble(l-1)!
      ! x(1) = 0.1141947d-1 -2.d0*(0.1d-2)+(0.1d-2)*dble(j-1)!
      ! x(2) = 0.8401042d-2 -2.d0*(0.1d-3)+(0.1d-3)*dble(k-1)!
!
! call residual(x,resout(j,k))
! write(*,*) x(5), resout(j,k)
! ! write(*,*) x(1), x(2), x(3), x(4), x(5), resout(1,1)!x(2), 
! write(*,*) resout(1,1)!x(2), 
! stop
! write(*,"(100(D11.2,x))") (x(i),i=1,n)
      call bobyqa(n,npt,x,xlo,xup,rhobeg,rhoend,iprint,maxfun,w)
! 
      write(*,"(i4x,i3x,i6,100(D15.6,2X))")expnum,ii,nst(ii), &
         fst(ii),(xst(ii,i),i=1,n)
      write(99,"(i4x,i3x,i6,100(D15.7,2X))")expnum,ii,nst(ii), &
         fst(ii),(xst(ii,i),i=1,n)
!
      end do !kk
      end do !jj
      end do !l

      end do !k
! write(222,*) (resout(j,1))!(resout(j,i),i=1,30)
      end do !ii
!
      indopt = minloc(fst,1)
!      print*, indopt
!
      write(*,"(5x,i3x,i6,100(D15.6,2X))")indopt,nst(indopt), &
         fst(indopt),(xst(indopt,i),i=1,n)
      ! write(99,"(i3x,i6,100(D15.6,2X))")indopt,nst(indopt),fst(indopt),(xst(indopt,i),i=1,n)
!
      end do !expnum
      ! write(1111,*) 'fim'
!
!     deallocating
!
      deallocate(mesh,meih,merh,mesv,meiv)
      deallocate(mdsh,mdih,mdrh,mdsv,mdiv)
      deallocate(suh,inh,reh,suv,inv)
      deallocate(obsinh)
      deallocate(intsuh,intinh,intreh,intsuv,intinv,intctl,intj,intnnc,intvac)
      deallocate(shold,ihold,rhold,svold,ivold)
      deallocate(l1,l2,l3,l4,l5)
      deallocate(l1old,l2old,l3old,l4old,l5old)
      deallocate(u,uold)
      deallocate(coord,el,ed)
!
      end program
!
!----------------------------------------------------------------------
!==========================================================================
