program tpcf
    implicit none
    
    real*8 :: rgrid_x, rgrid_y, rgrid_z
    real*8 :: boxsize_x, boxsize_y, boxsize_z
    real*8 :: disx, disy, disz, dis, vol, rhomed
    real*8 :: rwidth, dim1_max, dim1_min
    real*8 :: pi = 4.*atan(1.)
    
    integer*8 :: ntracers, ncentres, dim1_nbin, rind
    integer*8 :: i, ii, ix, iy, iz, ix2, iy2, iz2
    integer*8 :: indx, indy, indz, nrows, ncols
    integer*8 :: ipx, ipy, ipz, ndif_x, ndif_y, ndif_z
    integer*8 :: ngrid_x, ngrid_y, ngrid_z
    
    integer*8, dimension(:, :, :), allocatable :: lirst, nlirst
    integer*8, dimension(:), allocatable :: ll
    
    real*8, dimension(3) :: r
    real*8, allocatable, dimension(:,:)  :: tracers, centres
    real*8, dimension(:), allocatable :: DD, delta
    real*8, dimension(:), allocatable :: rbin, rbin_edges
  
    character(20), external :: str
    character(len=500) :: data_filename, data_filename_2, output_filename
    character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char
    character(len=10) :: boxchar_x, boxchar_y, boxchar_z
    character(len=10) :: ngridchar_x, ngridchar_y, ngridchar_z

    if (iargc() .ne. 12) then
        write(*,*) 'Some arguments are missing.'
        write(*,*) '1) data_filename'
        write(*,*) '2) data_filename_2'
        write(*,*) '3) output_filename'
        write(*,*) '4) boxsize_x'
        write(*,*) '5) boxsize_y'
        write(*,*) '6) boxsize_z'
        write(*,*) '7) dim1_min'
        write(*,*) '8) dim1_max'
        write(*,*) '9) dim1_nbin'
        write(*,*) '10) ngrid_x'
        write(*,*) '11) ngrid_y'
        write(*,*) '12) ngrid_z'
        write(*,*) ''
        stop
      end if
      
    call getarg(1, data_filename)
    call getarg(2, data_filename_2)
    call getarg(3, output_filename)
    call getarg(4, boxchar_x)
    call getarg(5, boxchar_y)
    call getarg(6, boxchar_z)
    call getarg(7, dim1_min_char)
    call getarg(8, dim1_max_char)
    call getarg(9, dim1_nbin_char)
    call getarg(10, ngridchar_x)
    call getarg(11, ngridchar_y)
    call getarg(12, ngridchar_z)
    
    read(boxchar_x, *) boxsize_x
    read(boxchar_y, *) boxsize_y
    read(boxchar_z, *) boxsize_z
    read(dim1_min_char, *) dim1_min
    read(dim1_max_char, *) dim1_max
    read(dim1_nbin_char, *) dim1_nbin
    read(ngridchar_x, *) ngrid_x
    read(ngridchar_y, *) ngrid_y
    read(ngridchar_z, *) ngrid_z
    
    write(*,*) '-----------------------'
    write(*,*) 'Running tpcf.exe'
    write(*,*) 'input parameters:'
    write(*,*) ''
    write(*, *) 'data_filename: ', trim(data_filename)
    write(*, *) 'data_filename_2: ', trim(data_filename_2)
    write(*, *) 'boxsize_x: ', trim(boxchar_x)
    write(*, *) 'boxsize_y: ', trim(boxchar_y)
    write(*, *) 'boxsize_z: ', trim(boxchar_z)
    write(*, *) 'output_filename: ', trim(output_filename)
    write(*, *) 'dim1_min: ', trim(dim1_min_char), ' Mpc'
    write(*, *) 'dim1_max: ', trim(dim1_max_char), ' Mpc'
    write(*, *) 'dim1_nbin: ', trim(dim1_nbin_char)
    write(*, *) 'ngrid_x: ', trim(ngridchar_x)
    write(*, *) 'ngrid_y: ', trim(ngridchar_y)
    write(*, *) 'ngrid_z: ', trim(ngridchar_z)
    write(*,*) ''
  
    open(10, file=data_filename, status='old', form='unformatted')
    read(10) nrows
    read(10) ncols
    allocate(tracers(ncols, nrows))
    read(10) tracers
    close(10)
    ntracers = nrows
    write(*,*) 'ntracers dim: ', size(tracers, dim=1), size(tracers, dim=2)
    write(*,*) 'pos(min), pos(max) = ', minval(tracers(1,:)), maxval(tracers(1,:))

    open(11, file=data_filename_2, status='old', form='unformatted')
    read(11) nrows
    read(11) ncols
    allocate(centres(ncols, nrows))
    read(11) centres
    close(11)
    ncentres = nrows
    write(*,*) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)
  
    allocate(rbin(dim1_nbin))
    allocate(rbin_edges(dim1_nbin + 1))
    allocate(DD(dim1_nbin))
    allocate(delta(dim1_nbin))
    
    rwidth = (dim1_max - dim1_min) / dim1_nbin
    do i = 1, dim1_nbin + 1
      rbin_edges(i) = dim1_min+(i-1)*rwidth
    end do
    do i = 1, dim1_nbin
      rbin(i) = rbin_edges(i+1)-rwidth/2.
    end do
    
    ! Mean density inside the box
    rhomed = ntracers / (boxsize_x * boxsize_y * boxsize_z)
    
    ! Construct linked list for tracers
    write(*,*) ''
    write(*,*) 'Constructing linked list...'
    allocate(lirst(ngrid_x, ngrid_y, ngrid_z))
    allocate(nlirst(ngrid_x, ngrid_y, ngrid_z))
    allocate(ll(ntracers))
    rgrid_x = boxsize_x / real(ngrid_x)
    rgrid_y = boxsize_y / real(ngrid_y)
    rgrid_z = boxsize_z / real(ngrid_z)
    
    lirst = 0
    ll = 0
    
    do i = 1, ntracers
      indx = int((tracers(1, i)) / rgrid_x + 1.)
      indy = int((tracers(2, i)) / rgrid_y + 1.)
      indz = int((tracers(3, i)) / rgrid_z + 1.)
    
      if(indx.gt.0.and.indx.le.ngrid_x.and.indy.gt.0.and.indy.le.ngrid_y.and.&
      indz.gt.0.and.indz.le.ngrid_z)lirst(indx,indy,indz)=i
    
      if(indx.gt.0.and.indx.le.ngrid_x.and.indy.gt.0.and.indy.le.ngrid_y.and.&
      indz.gt.0.and.indz.le.ngrid_z)nlirst(indx,indy,indz) = &
      nlirst(indx, indy, indz) + 1
    end do
    
    do i = 1, ntracers
      indx = int((tracers(1, i))/ rgrid_x + 1.)
      indy = int((tracers(2, i))/ rgrid_y + 1.)
      indz = int((tracers(3, i))/ rgrid_z + 1.)
      if(indx.gt.0.and.indx.le.ngrid_x.and.indy.gt.0.and.indy.le.ngrid_y.and.&
      &indz.gt.0.and.indz.le.ngrid_z) then
        ll(lirst(indx,indy,indz)) = i
        lirst(indx,indy,indz) = i
      endif
    end do
    
    write(*,*) 'Linked list successfully constructed'
    write(*,*) ''
    write(*,*) 'Starting loop over tracers...'
    
    DD = 0
    delta = 0
    
    do i = 1, ncentres
  
      ipx = int(centres(1, i) / rgrid_x + 1.)
      ipy = int(centres(2, i) / rgrid_y + 1.)
      ipz = int(centres(3, i) / rgrid_z + 1.)
  
      ndif_x = int(dim1_max / rgrid_x + 1.)
      ndif_y = int(dim1_max / rgrid_y + 1.)
      ndif_z = int(dim1_max / rgrid_z + 1.)
    
      do ix = ipx - ndif_x, ipx + ndif_x
        do iy = ipy - ndif_y, ipy + ndif_y
          do iz = ipz - ndif_z, ipz + ndif_z
    
            ix2 = ix
            iy2 = iy
            iz2 = iz
    
            if (ix2 .gt. ngrid_x) ix2 = ix2 - ngrid_x
            if (ix2 .lt. 1) ix2 = ix2 + ngrid_x
            if (iy2 .gt. ngrid_y) iy2 = iy2 - ngrid_y
            if (iy2 .lt. 1) iy2 = iy2 + ngrid_y
            if (iz2 .gt. ngrid_z) iz2 = iz2 - ngrid_z
            if (iz2 .lt. 1) iz2 = iz2 + ngrid_z
    
            ii = lirst(ix2,iy2,iz2)
            if(ii.ne.0) then
              do
                ii = ll(ii)
                disx = tracers(1, ii) - centres(1, i)
                disy = tracers(2, ii) - centres(2, i)
                disz = tracers(3, ii) - centres(3, i)
  
                if (disx .lt. -boxsize_x/2) disx = disx + boxsize_x
                if (disx .gt. boxsize_x/2) disx = disx - boxsize_x
                if (disy .lt. -boxsize_y/2) disy = disy + boxsize_y
                if (disy .gt. boxsize_y/2) disy = disy - boxsize_y
                if (disz .lt. -boxsize_z/2) disz = disz + boxsize_z
                if (disz .gt. boxsize_z/2) disz = disz - boxsize_z
    
                r = (/ disx, disy, disz /)
                dis = norm2(r)
  
                if (dis .gt. dim1_min .and. dis .lt. dim1_max) then
                  rind = int((dis - dim1_min) / rwidth + 1)
                  DD(rind) = DD(rind) + 1
                end if
    
                if(ii.eq.lirst(ix2,iy2,iz2)) exit
    
              end do
            end if
          end do
        end do
      end do
    end do
  
    do i = 1, dim1_nbin
      vol = 4./3 * pi * (rbin_edges(i + 1) ** 3 - rbin_edges(i) ** 3)
      delta(i) = DD(i) / (vol * rhomed * ncentres) - 1
    end do
    
    write(*,*) ''
    write(*,*) 'Calculation finished. Writing output...'
    
    open(12, file=output_filename, status='replace')
    do i = 1, dim1_nbin
      write(12, fmt='(2f15.5)') rbin(i), delta(i)
    end do
  
    end program tpcf
    