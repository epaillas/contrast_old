module procedures
    implicit none
contains

    subroutine detect_format(filename, formatted)
        character(*), intent(in) :: filename
        logical, intent(out) :: formatted
        integer :: fid, stat
        character :: c

        stat = 0
        formatted = .true. !assume formatted
        open (newunit=fid, file=filename, status='old', form='unformatted', recl=1)

        do while ((stat == 0) .and. formatted)
            read (fid, iostat=stat) c
            formatted = formatted .and. (iachar(c) <= 127)
        end do
        if (formatted) then
            print *, trim(filename), ' is a formatted file'
        else
            print *, trim(filename), ' is an unformatted file'
        end if
        close (fid)

    end subroutine detect_format

end module procedures
