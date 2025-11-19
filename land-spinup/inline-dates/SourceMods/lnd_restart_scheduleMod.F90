module lnd_restart_scheduleMod
   implicit none
   save

   integer, parameter :: n_restart_times = 1000
   integer :: restart_times(n_restart_times)
   integer :: n_loaded = 0

   public :: load_restart_schedule
   public :: restart_times
   public :: n_loaded
   public :: n_restart_times

contains

   subroutine load_restart_schedule(filename, iulog)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: iulog
      integer :: ierr, unit

      restart_times(:) = -1
      n_loaded = 0

      inquire(file=filename, exist=ierr)
      if (.not. ierr) then
         write(iulog,*) "Restart schedule: file not found: ", filename
         return
      end if

      open(newunit=unit, file=filename, status="old", action="read")

      do
         read(unit, *, iostat=ierr) restart_times(n_loaded+1)
         if (ierr /= 0) exit
         n_loaded = n_loaded + 1
         if (n_loaded == n_restart_times) exit
      end do

      close(unit)
   end subroutine load_restart_schedule

end module lnd_restart_scheduleMod
