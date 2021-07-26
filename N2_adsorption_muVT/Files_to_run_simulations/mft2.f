!----------------------------------------------------------------
! This code is modified by Dmitry N Lapshin.
! Modified on November 2020.
! The initial code is written by Peter A. Monson.
! Created on February 2007.
! Solution of the mean field kinetic equations for a 3D 
! lattice gas model in finite length slit pore with nearest 
! neighbor interactions.
! Calculations take advantage of the 2D symmetry of the slit.
! Simple cubic lattice.
!----------------------------------------------------------------


      implicit real*8(a-h,o-z)
      common/coords/ ix(1000000),iz(1000000)
      common/list/ nlist(4,1000000),ilist(1000000)
      common/params/ nsite,mx,mz,l
      common/density/ rhoold(1000000),rhonew(1000000)
      common/walls/ ieta(1000000),phi(1000000)
      character*2 t(10)
      character*12 fname1,fname2
      ! types of sites for visualizations
      t(1) = "C" 
      t(2) = "O" 
      t(3) = "H"
      t(4) = "N"
      t(5) = "S"
      t(6) = "P"
      t(7) = "Cl"
      t(8) = "Br"
      t(9) = "Na"
      t(10) = "Fe"
      ! open input file
      open(1,file='mft2.dat',status='unknown')
      ! read pore dimensions
      ! read interaction parameter and temperature
      read(1,*) y, tstar
      ! y - ratio of solid-fluid to fluid-fluid interactions
      ! tstar - kT/eff
      read(1,*) mx, mz, l, psite
      ! dimension of lattice is mx x mz x infinity (the pore width is mz-2) 
      ! l - length of the pore
      ! psite - total number of sites
      ! read simulation conditions
      read(1,*) ides
      read (1,*) ncp, cp0, dcp
      ! ncp - number of chemical potential values
      ! cp0 - initial chemical potential
      ! dcp - chemical potential increment
      read (1,*) rhog
      ! rhog - initial gas density
      ! read result convergence conditions
      read(1,*) maxit,error
      ! maxit - maximum number of iterations
      ! error - convergence criterion is that mean square change in
      ! local density between iterations is less than error
!----------------------------------------------------------------


      open(10,file='script.spt',status='unknown')
      write(10,*) "write on "
      write(10,*) "write true "
      ! ifile is a counter for determining the filenames on the density
      ! profile output files
      ifile=0
      ! coordination number 4 for 2D and 6 for 3D
      nc=4
      pi=3.141592654d0
      seed=1.986242345d0
      ! total number of sites
      nsite = mx*mz
!----------------------------------------------------------------


      ! initialize lattice coordinates and determine list of nearest neighbor 
      ! sites for each site on the lattice
      call lattice
!----------------------------------------------------------------


      ! initialize density distribution and external field
      do i=1,nsite
          is=ilist(i)
          phi(is)=0.0d0
          rhoold(is)=0.0d0
          if(ieta(is).eq.1) then
              if((iz(is).lt.27).and.(iz(is).gt.10)) then
              if((ix(is).gt.(mx-l)/2).and.(ix(is).le.(mx+l)/2)) then
                  do j=1,nc
                      jj=nlist(j,is)
                      phi(is)=phi(is)-(1.0-float(ieta(jj)))*y
                  end do
                  xi = (-phi(is)+cp0)/tstar
                  !rhoold(is) = 1.0d0/(1+exp(-xi))
                  rhoold(is) = rhog
              end if
              end if
          end if
      end do
!----------------------------------------------------------------


      ! calculate adsorption and desorption isotherms
      do iads=1,ides
          ! calculate chemical potential
          do icp = 1, ncp
              cp = cp0+(icp-1)*dcp
          ! here iteration begins
          iter = 0
1000      iter = iter + 1
          ! the number of iterations do not exceed maxit
          if(iter.gt.maxit) stop "Increase maxit"
          ! compute grand free energy and average density
          sumsq = 0.0d0
          sumrho = 0.0d0
          sumh0 = 0.0d0
          sumgp = 0.0d0
          sumrho1 = 0.0d0
          sumh01 = 0.0d0
          sumgp1 = 0.0d0
          sumf = 0.0d0
          do i=1,nsite
              is=ilist(i)
              if(ieta(is).eq.1) then
              ! Compute nearest neighbor sum for rhoj
                  if (nc.eq.4) then
                      snn1=2.0d0*rhoold(is)
                  else
                      snn=0
                  end if
                  do j=1, nc
                      jj=nlist(j,is)
                      if(ieta(jj).eq.1) then
                          snn1=snn1+rhoold(jj)
                      end if
                  end do
                  ! xi is the second term in the Hamiltonian
                  xi = (snn1-phi(is)+cp)/tstar
                  ! calculate new density
                  rhonew(is) = 1.0d0/(1.0d0+exp(-xi))
                  ! compute mean square change in local density
                  sumsq = sumsq+(rhonew(is)-rhoold(is))**2
                  ! total density and grand potential
                  ! over the whole pore length
                  if((iz(is).lt.27).and.(iz(is).gt.10)) then
                  if((ix(is).gt.(mx-l)/2).and.(ix(is).le.(mx+l)/2)) then
                      sumrho = sumrho+rhonew(is)
                      sumgp = sumgp-ieta(is)*tstar*log(1.0d0+exp(xi))
                      sumh0 = sumh0+rhonew(is)*snn1
                  end if
                  end if
                  ! total density and grand potential
                  ! at the middle (xy-plane) of the pore
                  if((iz(is).lt.27).and.(iz(is).gt.10)) then
                  if(ix(is).eq.mx/2) then
                      sumrho1 = sumrho1+rhonew(is)
                      sumgp1 = sumgp1-ieta(is)*tstar*log(1.0d0+exp(xi))
                      sumh01 = sumh01+rhonew(is)*snn1
                  end if
                  end if
              end if
          end do
          do i = 1,nsite
              is=ilist(i)
              if(ieta(is).eq.1) then
                  rhoold(is) = rhonew(is)
              endif
          end do
          ! normalize total density and grand potential
          sumsq = sumsq/nsite
          ! to the whole pore
          rhoav = sumrho/psite
          sumgp = sumgp/psite
          sumh0 = sumh0/2.0d0/psite
          gp = sumgp+sumh0
          ! to the middle (xy-plane) of the pore
          rhoav1 = sumrho1/(mz-2)
          sumgp1 = sumgp1/(mz-2)
          sumh01 = sumh01/2.0d0/(mz-2)
          gp1 = sumgp1+sumh01
          ! check for error
          if(sumsq.gt.error) go to 1000


          !print to screen
          print *, iter, sumsq
          print *, rhoav, gp, sumh0

          ! write isotherms to file
          open(2, file = 'isotherm_wpore.dat', status = 'unknown') ! for the whole pore
          open(4, file = 'isotherm_mpore.dat', status = 'unknown') ! for the middle of the pore
          write(2,'(5g16.8)') cp,exp((cp+3.0d0)/tstar),rhoav,gp
          write(4,'(5g16.8)') cp,exp((cp+3.0d0)/tstar),rhoav1,gp1

          ! write visualization files for each step
          ifile=ifile+1
          call namef(ifile,fname1,fname2)
          ! open(15, file = fname, status = 'unknown')
          open(3, file = fname1, status = 'unknown')
          ! write(15,*) nsite, rhoav ! another file
          write(3,*) psite+2*10*mz   ! the first line in rho*.dat files
          write(3,*) '#', rhoav      ! next empty line
          ! write sites to rho*.dat files
          do i=1, nsite
              is=ilist(i)
              if(ieta(is).eq.1) then
                  ! write(15,'(2i8,f16.10,i8)') ix(is),1,iz(is),rhonew(is),i
                  m = 10*rhonew(is)+1
                  write(3,*) t(m),ix(is),1,iz(is)
              end if
          end do
          ! close(15)
          close(3)
          ! script content
          write(10,*) "load xyz ", fname1
          write(10,*) "zoom 100"
          write(10,*) "select all"
          write(10,*) "spacefill 200"
          write(10,*) "select carbon"
          write(10,*) "color [250,250,250]"
          write(10,*) "select oxygen"
          write(10,*) "color [225,225,225]"
          write(10,*) "select hydrogen"
          write(10,*) "color [200,200,200]"
          write(10,*) "select nitrogen"
          write(10,*) "color [175,175,175]"
          write(10,*) "select sulfur"
          write(10,*) "color [150,150,150]"
          write(10,*) "select phosphorus"
          write(10,*) "color [125,125,125]"
          write(10,*) "select chlorine"
          write(10,*) "color [100,100,100]"
          write(10,*) "select bromine"
          write(10,*) "color [75,75,75]"
          write(10,*) "select sodium"
          write(10,*) "color [50,50,50]"
          write(10,*) "select iron"
          write(10,*) "color [25,25,25]"
          write(10,*) "set background [255,255,255]"
          write(10,*) "rotate x 90"
          write(10,*) "write ",fname2
          write(10,*) "zap"
          end do


          ! next chemical potential
          cp0 = cp-dcp
          dcp = -dcp
          ncp = ncp-1
      end do
      write(10,*) "quit"
      stop
      end
!----------------------------------------------------------------


      ! lattice model of the pore
      subroutine lattice

      implicit real*8(a-h,o-z)
      common/coords/ ix(1000000),iz(1000000)
      common/list/ nlist(4,1000000),ilist(1000000)
      common/params/ nsite,mx,mz,l
      common/density/ rhoold(1000000),rhonew(1000000)
      common/walls/ ieta(1000000),phi(1000000)
      integer, dimension(3110) :: PSD1, PSD2
      open(22,file='mft_psd.dat',status='unknown')
      read(22,*) PSD1
      read(22,*) PSD2

      nsite=0
      ! loop over mx and mz coordiates
      do i=1,mx
         do k=1,mz
           ! go through each site
           is=(i-1)*mz+k
           nsite=nsite+1
           ilist(nsite)=is
           ! write(*,*) i,k,is,nsite
           ix(is)=i
           iz(is)=k
           ! define sites of the pore
           ieta(is)=1
           do ipore=1,size(PSD1)
           upore = PSD1(ipore)
           lpore = PSD2(ipore)
           lwall = 1+(0.5*mz-lpore)
           ! if((i.ge.ipore).and.(i.lt.(1+1*ipore))) then
           ! if((k.le.(lwall-1)).or.(k.ge.(upore+lpore+lwall)))ieta(is)=0
           ! endif
           if((i.ge.ipore).and.(i.lt.(1+1*ipore))) then
           if((k.le.(lwall-1)).and.(k.ge.(11)))ieta(is)=0
           endif
           if((i.ge.ipore).and.(i.lt.(1+1*ipore))) then
           if((k.le.(26)).and.(k.ge.(upore+lpore+lwall)))ieta(is)=0
           endif
           if((i.ge.ipore).and.(i.lt.(1+1*ipore))) then
           if((k.le.(1)).or.(k.ge.(mz)))ieta(is)=0
           endif
           ! PBC
           i1=i+1
           if(i1.gt.mx) i1=1
           i2=i-1
           if(i2.lt.1) i2=mx
           k1=k+1
           if(k1.gt.mz) k1=1
           k2=k-1
           if(k2.lt.1) k2=mz
           ! list of nearest neighbors
           nlist(1,is)=(i1-1)*mz+k
           nlist(2,is)=(i2-1)*mz+k
           nlist(3,is)=(i-1)*mz+k1
           nlist(4,is)=(i-1)*mz+k2
           end do
         end do
      end do
      return
      end
!----------------------------------------------------------------


      subroutine rhobulk(cpb,rhob,tstar)
      implicit real*8(a-h,o-z)

      act=exp(cpb/tstar)
      rhob=act

      do i=1, 1000
          rhob=act/(act+exp(-6.0*rhob/tstar))
      end do
      return
      end
!----------------------------------------------------------------


      subroutine namef(i,fname1,fname2)
      character num(3)
      character*12 fname1,fname2
      character*3 string

      num(1)='0'
      num(2)='0'
      num(3)='0'

      if(i.lt.10) then
          write(string,'(i1)')i
          read(string,'(3a1)') num(3)
      else if(i.lt.100) then
          write(string,'(i2)')i
          read(string,'(3a1)') num(2),num(3)
      else
          write(string,'(i3)')i
          read(string,'(3a1)') num(1),num(2),num(3)
      endif
      ! file name
      fname1='rho'//num(1)//num(2)//num(3)//'.dat'
      fname2='rho'//num(1)//num(2)//num(3)//'.gif'
      return
      end
!----------------------------------------------------------------


      function ranf(seed)
c     random no. generator
      real*8 seed,d2p31m,d2p31
      data d2p31m/2147483647.0d0/
      data d2p31/2147483648.0d0/

      seed = dmod(16807.0d0*seed,d2p31m)
      ranf = seed/d2p31
      return
      end