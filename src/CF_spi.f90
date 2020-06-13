program CCF_spi
  implicit none
  
  real*8 :: rgrid, boxsize, vol, rhomed
  real*8 :: disx, disy, disz, para, perp
  real*8 :: xvc, yvc, zvc
  real*8 :: perpwidth, dim1_min, dim1_max
  real*8 :: parawidth, dim2_min, dim2_max
  real*8 :: pi = 4.*atan(1.)
  
  integer*8 :: ng, dim1_nbin, perpind, dim2_nbin, paraind
  integer*8 :: i, j, ii, jj, ix, iy, iz, ix2, iy2, iz2
  integer*8 :: indx, indy, indz, nrows, ncols
  integer*8 :: ipx, ipy, ipz, ndif
  integer*8 :: ngrid
  
  integer*8, dimension(:, :, :), allocatable :: lirst, nlirst
  integer*8, dimension(:), allocatable :: ll
  
  real*8, dimension(3) :: r, com
  real*8, allocatable, dimension(:,:)  :: tracers
  real*8, dimension(:, :), allocatable :: DD, delta
  real*8, dimension(:), allocatable :: perpbin, perpbin_edges, parabin, parabin_edges

  logical :: has_velocity = .false.
  
  character(20), external :: str
  character(len=500) :: data_filename, output_filename
  character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char, ngrid_char, box_char
  
  if (iargc() .ne. 7) then
      write(*,*) 'Some arguments are missing.'
      write(*,*) '1) data_filename'
      write(*,*) '2) output_filename'
      write(*,*) '3) boxsize'
      write(*,*) '4) dim1_min'
      write(*,*) '5) dim1_max'
      write(*,*) '6) dim1_nbin'
      write(*,*) '7) ngrid'
      write(*,*) ''
      stop
    end if
    
  call getarg(1, data_filename)
  call getarg(2, output_filename)
  call getarg(3, box_char)
  call getarg(4, dim1_min_char)
  call getarg(5, dim1_max_char)
  call getarg(6, dim1_nbin_char)
  call getarg(7, ngrid_char)
  
  read(box_char, *) boxsize
  read(dim1_min_char, *) dim1_min
  read(dim1_max_char, *) dim1_max
  read(dim1_nbin_char, *) dim1_nbin
  read(ngrid_char, *) ngrid
  
  write(*,*) '-----------------------'
  write(*,*) 'Running CCF_spi.exe'
  write(*,*) 'input parameters:'
  write(*,*) ''
  write(*, *) 'data_filename: ', trim(data_filename)
  write(*, *) 'boxsize: ', trim(box_char)
  write(*, *) 'output_filename: ', trim(output_filename)
  write(*, *) 'dim1_min: ', trim(dim1_min_char), ' Mpc'
  write(*, *) 'dim1_max: ', trim(dim1_max_char), ' Mpc'
  write(*, *) 'dim1_nbin: ', trim(dim1_nbin_char)
  write(*, *) 'ngrid: ', trim(ngrid_char)
  write(*,*) ''

  open(10, file=data_filename, status='old', form='unformatted')
  read(10) nrows
  read(10) ncols
  allocate(tracers(ncols, nrows))
  read(10) tracers
  close(10)
  ng = nrows
  if (ncols .eq. 6) then
    has_velocity = .true.
    write(*,*) 'Tracer file has velocity information.'
  end if
  write(*,*) 'ntracers dim: ', size(tracers, dim=1), size(tracers, dim=2)
  write(*,*) 'pos(min), pos(max) = ', minval(tracers(1,:)), maxval(tracers(1,:))

  dim2_nbin = dim1_nbin
  allocate(perpbin(dim1_nbin))
  allocate(parabin(dim2_nbin))
  allocate(perpbin_edges(dim1_nbin + 1))
  allocate(parabin_edges(dim2_nbin + 1))
  allocate(DD(dim1_nbin, dim2_nbin))
  allocate(delta(dim1_nbin, dim2_nbin))
  
  
  perpwidth = (dim1_max - dim1_min) / dim1_nbin
  do i = 1, dim1_nbin + 1
    perpbin_edges(i) = dim1_min+(i-1)*perpwidth
  end do
  do i = 1, dim1_nbin
    perpbin(i) = perpbin_edges(i+1)-perpwidth/2.
  end do

  dim2_min = dim1_min
  dim2_max = dim1_max

  parawidth = (dim2_max - dim2_min) / dim2_nbin
  do i = 1, dim2_nbin + 1
    parabin_edges(i) = dim2_min+(i-1)*parawidth
  end do
  do i = 1, dim2_nbin
    parabin(i) = parabin_edges(i+1)-parawidth/2.
  end do
  
  ! Mean density inside the box
  rhomed = ng / (boxsize ** 3)
  
  ! Construct linked list for tracers
  write(*,*) ''
  write(*,*) 'Constructing linked list...'
  allocate(lirst(ngrid, ngrid, ngrid))
  allocate(nlirst(ngrid, ngrid, ngrid))
  allocate(ll(ng))
  rgrid = (boxsize) / real(ngrid)
  
  lirst = 0
  ll = 0
  
  do i = 1, ng
    indx = int((tracers(1, i)) / rgrid + 1.)
    indy = int((tracers(2, i)) / rgrid + 1.)
    indz = int((tracers(3, i)) / rgrid + 1.)
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    indz.gt.0.and.indz.le.ngrid)lirst(indx,indy,indz)=i
  
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    indz.gt.0.and.indz.le.ngrid)nlirst(indx,indy,indz) = &
    nlirst(indx, indy, indz) + 1
  end do
  
  do i = 1, ng
    indx = int((tracers(1, i))/ rgrid + 1.)
    indy = int((tracers(2, i))/ rgrid + 1.)
    indz = int((tracers(3, i))/ rgrid + 1.)
    if(indx.gt.0.and.indx.le.ngrid.and.indy.gt.0.and.indy.le.ngrid.and.&
    &indz.gt.0.and.indz.le.ngrid) then
      ll(lirst(indx,indy,indz)) = i
      lirst(indx,indy,indz) = i
    endif
  end do
  
  write(*,*) 'Linked list successfully constructed'
  write(*,*) ''
  write(*,*) 'Starting loop over centres...'
  
  DD = 0
  delta = 0
  
  do i = 1, ng
    xvc = tracers(1, i)
    yvc = tracers(2, i)
    zvc = tracers(3, i)

    ipx = int((xvc) / rgrid + 1.)
    ipy = int((yvc) / rgrid + 1.)
    ipz = int((zvc) / rgrid + 1.)

    ndif = int(dim1_max / rgrid + 1.)
  
    do ix = ipx - ndif, ipx + ndif
      do iy = ipy - ndif, ipy + ndif
        do iz = ipz - ndif, ipz + ndif
  
          ix2 = ix
          iy2 = iy
          iz2 = iz
  
          if (ix2 .gt. ngrid) ix2 = ix2 - ngrid
          if (ix2 .lt. 1) ix2 = ix2 + ngrid
          if (iy2 .gt. ngrid) iy2 = iy2 - ngrid
          if (iy2 .lt. 1) iy2 = iy2 + ngrid
          if (iz2 .gt. ngrid) iz2 = iz2 - ngrid
          if (iz2 .lt. 1) iz2 = iz2 + ngrid
  
          ii = lirst(ix2,iy2,iz2)
          if(ii.ne.0) then
            do
              ii = ll(ii)
              disx = tracers(1, ii) - xvc
              disy = tracers(2, ii) - yvc
              disz = tracers(3, ii) - zvc

              if (disx .lt. -boxsize/2) disx = disx + boxsize
              if (disx .gt. boxsize/2) disx = disx - boxsize
              if (disy .lt. -boxsize/2) disy = disy + boxsize
              if (disy .gt. boxsize/2) disy = disy - boxsize
              if (disz .lt. -boxsize/2) disz = disz + boxsize
              if (disz .gt. boxsize/2) disz = disz - boxsize
  
              r = (/ disx, disy, disz /)
              com = (/ 0, 0, 1 /)
              para = abs(dot_product(r, com))
              perp = sqrt(norm2(r) ** 2 - para ** 2)

              if (perp .gt. dim1_min .and. perp .lt. dim1_max .and. &
                  para .gt. dim2_min .and. para .lt. dim2_max ) then
                perpind = int((perp - dim1_min) / perpwidth + 1)
                paraind = int((para - dim2_min) / parawidth + 1)
                DD(perpind, paraind) = DD(perpind, paraind) + 1
              end if

  
              if(ii.eq.lirst(ix2,iy2,iz2)) exit
  
            end do
          end if
        end do
      end do
    end do
  end do

  do ii = 1, dim1_nbin
    do jj = 1, dim2_nbin
      vol = 2 * pi * (parabin_edges(jj + 1) - parabin_edges(jj))&
                 * (perpbin_edges(ii + 1)**2 - perpbin_edges(ii)**2)
                 
      delta(ii, jj) = DD(ii, jj) / (vol * rhomed * ng) - 1
    end do
  end do
  
  write(*,*) ''
  write(*,*) 'Calculation finished. Writing output...'
  
  open(12, file=output_filename, status='replace')
  do j = 1, dim2_nbin
    do i = 1, dim1_nbin
      write(12, fmt='(3f15.5)') perpbin(i), parabin(j), delta(i, j)
    end do
  end do

  end program CCF_spi
  