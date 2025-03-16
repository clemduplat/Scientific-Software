module timings
    implicit none
    save
    private
    public tic,toc,startClock,stopClock
    integer, parameter :: dp = selected_real_kind(15,307)
    integer, parameter :: sp = selected_real_kind(6,37)
    real(sp) :: start_CPU_time=0.0_sp
    real(sp) :: start_real_time=0.0_sp
    
contains

    subroutine tic(startTime)
        ! tic() starts the cpu timer. The subroutines should record the processor time at execution of this command in a module variable. Use toc()/ toc(elapsedTime) to print/get the elapsed time since the last call of tic().
         ! tic(startTime) returns the value of the processor time at execution of this command in the real(kind(0.e0)) variable startTime. Calling this command shall not change the module variable associated with calling tic() (ie. without argument).
        !-----Optional start time----!
        real(sp), optional, intent(out) :: startTime
        !----Start CPU timer----!
        if(present(startTime)) then
            call cpu_time(startTime)
        else
            !--should record the CPU time at the moment of executing this comment in a module variable---!
            call cpu_time(start_CPU_time)
        endif
        
    end subroutine tic

    subroutine toc(elapsedTime, startTime)
	! toc() prints the elapsed cpu time in seconds since the most recent call of tic() (ie. tic called without output argument).
        ! toc(elapsedTime) returns the elapsed cpu time in seconds since the most recent call of tic() (ie. tic called without output argument) in the real(kind(0.e0)) variable elapsedTime.
        ! if toc() or toc(elapsedTime) is called without first calling tic(), a meaningless result may be returned.
        ! toc(startTime=startTime) prints the elapsed cpu time in seconds since the call of the tic command corresponding to startTime.
        ! toc(elapsedTime, startTime) returns the elapsed cpu time in seconds since the call of the tic command corresponding to startTime in the real(kind(0.e0)) variable elapsedTime.
        !-----Two optional arguments----!
        real(sp), optional, intent(out) :: elapsedTime
        real(sp), optional, intent(in) ::startTime 
        real(sp) :: endTime
        call cpu_time(endTime)
        if(present(startTime)) then
            
            if(present(elapsedTime)) then
                !---Case elapsedTime present: value should not be printed but returned via variable elapsedTime-------!
                elapsedTime=endTime-startTime
            else
                !-----elapsedTime since last tic() call corresponding to startTime should be returned-----!
                print*,"ElapsedTime since last tic() call corresponding to startTime", endTime-startTime
            endif
        else
            !----Case where startTime is not present-> using of start_CPU_time-----!
            if(present(elapsedTime)) then
                elapsedTime=endTime-start_CPU_time
            else
                !----Case where neither arguments are present---!
                !----Print elapsedCPU time since most recent call of tic()----!
                print*,"ElapsedCPU time since most recent call of tic()",endTime-start_CPU_time
            endif
        endif
        if(elapsedTime<0) then
            print*,"ElapsedTime negatif"
        endif
    end subroutine toc

    subroutine startClock(startTime)
        ! startClock() starts the wall clock timer. The subroutine should record the value of a real-time clock at the execution of this command in a module variable. Use stopClock()/stopClock(elapsedTime) to print/get the elapsed time.
        ! startClock(startTime) stores the value of a real-time clock at execution of this command in the integer(kind(0)) variable startTime. Calling this command shall not change the module variable associated with calling startClock() (ie.startClock called without argument).
        integer(sp), optional, intent(out) :: startTime
        integer :: COUNT
        real(sp) :: COUNT_RATE ! determines the number of clock ticks per seconds
        call system_clock(COUNT, COUNT_RATE)
        if(present(startTime)) then
            startTime=real(COUNT,kind=sp)/COUNT_RATE
        else
            !----Case where StartTime is not present-> use of module variable----!
            start_real_time=real(COUNT,kind=sp)/COUNT_RATE
        endif
    end subroutine startClock
    
    subroutine stopClock(elapsedTime, startTime)
	 ! stopClock() prints the elapsed wall clock time in milliseconds since the most recent call of startClock() (ie. startClock called without output argument).
        ! stopClock(elapsedTime) returns the elapsed wall clock time in milliseconds since the most recent call of startClock() (ie. startClock called without output argument) in the integer(kind(0)) variable  elapsedTime.
        ! if stopClock() or stopClock(elapsedTime) is called without first calling startClock() a meaningless result may be returned.
        ! stopClock(startTime=startTime) prints the elapsed wall clock time in milliseconds since the call of the startClock command corresponding to startTime.
        ! stopClock(elapsedTime, startTime) returns the elapsed wall clock time in milliseconds since the call of the startClock command corresponding to startTime in the integer(kind(0)) variable elapsedTime.
        integer(sp), optional, intent(out) :: elapsedTime
        real(sp), optional, intent(in) ::startTime 
        integer :: COUNT 
        real(sp) :: COUNT_RATE ! determines the number of clock ticks per seconds
        real(sp) :: endTime
         
        !-----Get current time----!
        call system_clock(COUNT,COUNT_RATE)
        endTime=real(COUNT,kind=sp)/COUNT_RATE
        if(present(startTime)) then
            
            if(present(elapsedTime)) then
                !---Case elapsedTime present: value should not be printed but returned via variable elapsedTime-------!
                elapsedTime=(endTime-startTime)*1000
            else
                !-----elapsedTime since last tic() call corresponding to startTime should be returned-----!
                print*,"ElapsedTime since last tic() call corresponding to startTime", (endTime-startTime)*1000
            endif
        else
            !----Case where startTime is not present-> using of start_real_time-----!
            if(present(elapsedTime)) then
                elapsedTime=(endTime-start_real_time)*1000
            else
                !----Case where neither arguments are present---!
                !----Print elapsedCPU time since most recent call of tic()----!
                print*,"Elapsed real time since most recent call of tic()",(endTime-start_real_time)*1000
            endif
        endif

    end subroutine stopClock
end