program CF_monopole
  implicit none
  
  real*8 :: rgrid, boxsize, vol, int_vol, rhomed
  real*8 :: disx, disy, disz, dis, vr, vt, vlos
  real*8 :: xvc, yvc, zvc
  real*8 :: velx, vely, velz
  real*8 :: rwidth, dim1_max, dim1_min
  real*8 :: pi = 4.*atan(1.)
  
  integer*8 :: ng, dim1_nbin, rind
  integer*8 :: i, ii, ix, iy, iz, ix2, iy2, iz2
  integer*8 :: indx, indy, indz, nrows, ncols
  integer*8 :: ipx, ipy, ipz, ndif
  integer*8 :: ngrid
  
  integer*8, dimension(:, :, :), allocatable :: lirst, nlirst
  integer*8, dimension(:), allocatable :: ll
  
  real*8, dimension(3) :: r, vel, com
  real*8, allocatable, dimension(:,:)  :: tracers
  real*8, dimension(:), allocatable :: DD, int_DD, delta, int_delta
  real*8, dimension(:), allocatable :: VV_r, VV2_r, mean_vr, std_vr
  real*8, dimension(:), allocatable :: VV_t, VV2_t, mean_vt, std_vt
  real*8, dimension(:), allocatable :: VV_los, VV2_los, mean_vlos, std_vlos
  real*8, dimension(:), allocatable :: rbin, rbin_edges

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
  write(*,*) 'Running CF_monopole.exe'
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

  allocate(rbin(dim1_nbin))
  allocate(rbin_edges(dim1_nbin + 1))
  allocate(DD(dim1_nbin))
  allocate(int_DD(dim1_nbin))
  allocate(delta(dim1_nbin))
  allocate(int_delta(dim1_nbin))
  if (has_velocity) then
    allocate(VV_r(dim1_nbin))
    allocate(VV_t(dim1_nbin))
    allocate(VV_los(dim1_nbin))
    allocate(VV2_r(dim1_nbin))
    allocate(VV2_t(dim1_nbin))
    allocate(VV2_los(dim1_nbin))
    allocate(mean_vr(dim1_nbin))
    allocate(mean_vt(dim1_nbin))
    allocate(mean_vlos(dim1_nbin))
    allocate(std_vr(dim1_nbin))
    allocate(std_vt(dim1_nbin))
    allocate(std_vlos(dim1_nbin))
  end if
  
  
  rwidth = (dim1_max - dim1_min) / dim1_nbin
  do i = 1, dim1_nbin + 1
    rbin_edges(i) = dim1_min+(i-1)*rwidth
  end do
  do i = 1, dim1_nbin
    rbin(i) = rbin_edges(i+1)-rwidth/2.
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
  write(*,*) 'Starting loop over tracers...'
  
  DD = 0
  int_DD = 0
  delta = 0
  int_DD = 0
  if (has_velocity) then
    VV_r = 0
    VV_t = 0
    VV_los = 0
    VV2_r = 0
    VV2_t = 0
    VV2_los = 0
    mean_vr = 0
    mean_vt = 0
    std_vr = 0
    std_vt = 0
    mean_vlos = 0
    std_vlos = 0
  end if
  
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
              disx = tracers(1, ii) - tracers(1, i)
              disy = tracers(2, ii) - tracers(2, i)
              disz = tracers(3, ii) - tracers(3, i)

              if (disx .lt. -boxsize/2) disx = disx + boxsize
              if (disx .gt. boxsize/2) disx = disx - boxsize
              if (disy .lt. -boxsize/2) disy = disy + boxsize
              if (disy .gt. boxsize/2) disy = disy - boxsize
              if (disz .lt. -boxsize/2) disz = disz + boxsize
              if (disz .gt. boxsize/2) disz = disz - boxsize
  
              r = (/ disx, disy, disz /)
              dis = norm2(r)

              if (has_velocity) then
                velx = tracers(4, ii) - tracers(4, i)
                vely = tracers(5, ii) - tracers(5, i)
                velz = tracers(6, ii) - tracers(6, i)
                vel = (/ velx, vely, velz /)
                com = (/ 0, 0, 1 /)
                vr = dot_product(vel, r) / norm2(r)
                vt = norm2(vel - vr * r / norm2(r))
                vlos = dot_product(vel, com)
              end if

              if (dis .gt. dim1_min .and. dis .lt. dim1_max) then
                rind = int((dis - dim1_min) / rwidth + 1)
                DD(rind) = DD(rind) + 1

                if (has_velocity) then
                  VV_r(rind) = VV_r(rind) + vr
                  VV_t(rind) = VV_t(rind) + vt
                  VV_los(rind) = VV_los(rind) + vlos
                  VV2_r(rind) = VV2_r(rind) + vr**2
                  VV2_t(rind) = VV2_t(rind) + vt**2
                  VV2_los(rind) = VV2_los(rind) + vlos**2
                end if
              end if
  
              if(ii.eq.lirst(ix2,iy2,iz2)) exit
  
            end do
          end if
        end do
      end do
    end do
  end do

  int_DD(1) = DD(1)
  do i = 2, dim1_nbin
    int_DD(i) =  int_DD(i - 1) + DD(i)
  end do

  do i = 1, dim1_nbin
    vol = 4./3 * pi * (rbin_edges(i + 1) ** 3 - rbin_edges(i) ** 3)
    int_vol = 4./3 * pi * rbin_edges(i + 1) ** 3
    delta(i) = DD(i) / (vol * rhomed * ng) - 1
    int_delta(i) = int_DD(i) / (int_vol * rhomed * ng) - 1

    if (has_velocity) then
      mean_vr(i) = VV_r(i) / DD(i)
      mean_vt(i) = VV_t(i) / DD(i)
      mean_vlos(i) = VV_los(i) / DD(i)
      std_vr(i) = sqrt((VV2_r(i) - (VV_r(i) ** 2 / DD(i))) / (DD(i) - 1))
      std_vt(i) = sqrt((VV2_t(i) - (VV_t(i) ** 2 / DD(i))) / (DD(i) - 1))
      std_vlos(i) = sqrt((VV2_los(i) - (VV_los(i) ** 2 / DD(i))) / (DD(i) - 1))
    end if

  end do
  
  write(*,*) ''
  write(*,*) 'Calculation finished. Writing output...'
  
  open(12, file=output_filename, status='replace')
  do i = 1, dim1_nbin
    if (has_velocity) then
      write(12, fmt='(9f15.5)') rbin(i), delta(i), int_delta(i), mean_vr(i),&
      & std_vr(i), mean_vt(i), std_vt(i), mean_vlos(i), std_vlos(i)
    else
      write(12, fmt='(3f15.5)') rbin(i), delta(i), int_delta(i)
    end if
  end do

  end program CF_monopole
  