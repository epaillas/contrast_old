program s_mu_tpcf
    implicit none

    real*8 :: rgrid, boxsize, vol, rhomed
    real*8 :: disx, disy, disz, dis, mu
    real*8 :: rwidth, dim1_max, dim1_min
    real*8 :: muwidth, mumin, mumax
    real*8 :: pi = 4.*atan(1.)

    integer*8 :: ntracers, ncentres, dim1_nbin, rind, dim2_nbin, muind
    integer*8 :: i, j, ii, ix, iy, iz, ix2, iy2, iz2
    integer*8 :: indx, indy, indz, nrows, ncols
    integer*8 :: ipx, ipy, ipz, ndif
    integer*8 :: ngrid

    integer*8, dimension(:, :, :), allocatable :: lirst, nlirst
    integer*8, dimension(:), allocatable :: ll

    real*8, dimension(3) :: r, com
    real*8, allocatable, dimension(:, :)  :: tracers, centres
    real*8, dimension(:, :), allocatable :: DD, delta
    real*8, dimension(:), allocatable :: rbin, rbin_edges, mubin, mubin_edges

    character(20), external :: str
    character(len=500) :: input_tracers, centres_filename, output_filename
    character(len=10) :: dim1_max_char, dim1_min_char, dim1_nbin_char, dim2_nbin_char, ngrid_char, box_char

    if (iargc() .ne. 9) then
        write (*, *) 'Some arguments are missing.'
        write (*, *) '1) data_filename'
        write (*, *) '2) centres_filename'
        write (*, *) '3) output_filename'
        write (*, *) '4) boxsize'
        write (*, *) '5) dim1_min'
        write (*, *) '6) dim1_max'
        write (*, *) '7) dim1_nbin'
        write (*, *) '8) dim2_nbin'
        write (*, *) '9) ngrid'
        write (*, *) ''
        stop
    end if

    call getarg(1, input_tracers)
    call getarg(2, centres_filename)
    call getarg(3, output_filename)
    call getarg(4, box_char)
    call getarg(5, dim1_min_char)
    call getarg(6, dim1_max_char)
    call getarg(7, dim1_nbin_char)
    call getarg(8, dim2_nbin_char)
    call getarg(9, ngrid_char)

    read (box_char, *) boxsize
    read (dim1_min_char, *) dim1_min
    read (dim1_max_char, *) dim1_max
    read (dim1_nbin_char, *) dim1_nbin
    read (dim2_nbin_char, *) dim2_nbin
    read (ngrid_char, *) ngrid

    write (*, *) '-----------------------'
    write (*, *) 'Running s_mu_tpcf.exe'
    write (*, *) 'input parameters:'
    write (*, *) ''
    write (*, *) 'input_tracers: ', trim(input_tracers)
    write (*, *) 'centres_filename: ', trim(centres_filename)
    write (*, *) 'boxsize: ', trim(box_char)
    write (*, *) 'output_filename: ', trim(output_filename)
    write (*, *) 'dim1_min: ', trim(dim1_min_char), ' Mpc'
    write (*, *) 'dim1_max: ', trim(dim1_max_char), ' Mpc'
    write (*, *) 'dim1_nbin: ', trim(dim1_nbin_char)
    write (*, *) 'ngrid: ', trim(ngrid_char)
    write (*, *) ''

    open (10, file=input_tracers, status='old', form='unformatted')
    read (10) nrows
    read (10) ncols
    allocate (tracers(ncols, nrows))
    read (10) tracers
    close (10)
    ntracers = nrows
    write (*, *) 'ntracers dim: ', size(tracers, dim=1), size(tracers, dim=2)
    write (*, *) 'pos(min), pos(max) = ', minval(tracers(1, :)), maxval(tracers(1, :))

    open (11, file=centres_filename, status='old', form='unformatted')
    read (11) nrows
    read (11) ncols
    allocate (centres(ncols, nrows))
    read (11) centres
    close (11)
    ncentres = nrows
    write (*, *) 'ncentres dim: ', size(centres, dim=1), size(centres, dim=2)

    allocate (rbin(dim1_nbin))
    allocate (mubin(dim2_nbin))
    allocate (rbin_edges(dim1_nbin + 1))
    allocate (mubin_edges(dim2_nbin + 1))
    allocate (DD(dim1_nbin, dim2_nbin))
    allocate (delta(dim1_nbin, dim2_nbin))

    rwidth = (dim1_max - dim1_min)/dim1_nbin
    do i = 1, dim1_nbin + 1
        rbin_edges(i) = dim1_min + (i - 1)*rwidth
    end do
    do i = 1, dim1_nbin
        rbin(i) = rbin_edges(i + 1) - rwidth/2.
    end do

    mumin = -1
    mumax = 1

    muwidth = (mumax - mumin)/dim2_nbin
    do i = 1, dim2_nbin + 1
        mubin_edges(i) = mumin + (i - 1)*muwidth
    end do
    do i = 1, dim2_nbin
        mubin(i) = mubin_edges(i + 1) - muwidth/2.
    end do

    ! Mean density inside the box
    rhomed = ntracers/(boxsize**3)

    ! Construct linked list for tracers
    write (*, *) ''
    write (*, *) 'Constructing linked list...'
    allocate (lirst(ngrid, ngrid, ngrid))
    allocate (nlirst(ngrid, ngrid, ngrid))
    allocate (ll(ntracers))
    rgrid = (boxsize)/real(ngrid)

    lirst = 0
    ll = 0

    do i = 1, ntracers
        indx = int((tracers(1, i))/rgrid + 1.)
        indy = int((tracers(2, i))/rgrid + 1.)
        indz = int((tracers(3, i))/rgrid + 1.)

        if (indx .gt. 0 .and. indx .le. ngrid .and. indy .gt. 0 .and. indy .le. ngrid .and. &
            indz .gt. 0 .and. indz .le. ngrid) lirst(indx, indy, indz) = i

        if (indx .gt. 0 .and. indx .le. ngrid .and. indy .gt. 0 .and. indy .le. ngrid .and. &
            indz .gt. 0 .and. indz .le. ngrid) nlirst(indx, indy, indz) = &
            nlirst(indx, indy, indz) + 1
    end do

    do i = 1, ntracers
        indx = int((tracers(1, i))/rgrid + 1.)
        indy = int((tracers(2, i))/rgrid + 1.)
        indz = int((tracers(3, i))/rgrid + 1.)
        if (indx .gt. 0 .and. indx .le. ngrid .and. indy .gt. 0 .and. indy .le. ngrid .and.&
        &indz .gt. 0 .and. indz .le. ngrid) then
            ll(lirst(indx, indy, indz)) = i
            lirst(indx, indy, indz) = i
        endif
    end do

    write (*, *) 'Linked list successfully constructed'
    write (*, *) ''
    write (*, *) 'Starting loop over centres...'

    DD = 0
    delta = 0

    do i = 1, ncentres

        ipx = int(centres(1, i)/rgrid + 1.)
        ipy = int(centres(2, i)/rgrid + 1.)
        ipz = int(centres(3, i)/rgrid + 1.)

        ndif = int(dim1_max/rgrid + 1.)

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

                    ii = lirst(ix2, iy2, iz2)
                    if (ii .ne. 0) then
                        do
                            ii = ll(ii)
                            disx = tracers(1, ii) - centres(1, i)
                            disy = tracers(2, ii) - centres(2, i)
                            disz = tracers(3, ii) - centres(3, i)

                            if (disx .lt. -boxsize/2) disx = disx + boxsize
                            if (disx .gt. boxsize/2) disx = disx - boxsize
                            if (disy .lt. -boxsize/2) disy = disy + boxsize
                            if (disy .gt. boxsize/2) disy = disy - boxsize
                            if (disz .lt. -boxsize/2) disz = disz + boxsize
                            if (disz .gt. boxsize/2) disz = disz - boxsize

                            r = (/disx, disy, disz/)
                            com = (/0, 0, 1/)
                            dis = norm2(r)
                            mu = dot_product(r, com)/(norm2(r)*norm2(com))

                            if (dis .gt. dim1_min .and. dis .lt. dim1_max) then
                                rind = int((dis - dim1_min)/rwidth + 1)
                                muind = int((mu - mumin)/muwidth + 1)
                                DD(rind, muind) = DD(rind, muind) + 1
                            end if

                            if (ii .eq. lirst(ix2, iy2, iz2)) exit

                        end do
                    end if
                end do
            end do
        end do
    end do

    do i = 1, dim1_nbin
        do j = 1, dim2_nbin
            vol = 4./3*pi*(rbin_edges(i + 1)**3 - rbin_edges(i)**3)/(dim2_nbin)
            delta(i, j) = DD(i, j)/(vol*rhomed*ncentres) - 1
        end do
    end do

    write (*, *) ''
    write (*, *) 'Calculation finished. Writing output...'

    open (12, file=output_filename, status='replace')
    do j = 1, dim2_nbin
        do i = 1, dim1_nbin
            write (12, fmt='(3f15.5)') rbin(i), mubin(j), delta(i, j)
        end do
    end do

end program s_mu_tpcf
