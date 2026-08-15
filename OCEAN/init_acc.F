        subroutine initialisation_openacc
#include "cppdefs.h"
#ifdef OPENACC
# ifdef JEANZAY
        USE openacc
        implicit none
#  include "param.h"
 
        character(len=6) :: local_rank_env
        integer          :: local_rank_env_status
        integer          :: ngpus, local_rank, devicenum
 
        !  OpenACC initialization
        !$acc init
 
        ! Retrieve local process rank via environment variable
        ! environment variable set by Slurm, as MPI_Comm_rank is not yet available.
        ! possible as this routine is used BEFORE MPI initialization.
!        call get_environment_variable(name="SLURM_LOCALID", &
        call get_environment_variable(name="SLURM_PROCID", 
     &      value=local_rank_env, status=local_rank_env_status)
 
        if (local_rank_env_status == 0) then
            read(local_rank_env, *) local_rank
                ngpus=acc_get_num_devices( acc_device_nvidia )
                !choose device to use by this rank
                devicenum = MOD( local_rank, ngpus )
            ! Defining the GPU to be used via OpenACC
            call acc_set_device_num( devicenum, acc_get_device_type())
            my_acc_device=devicenum
            compute_on_device=.TRUE.
!            write(100+local_rank,*) devicenum,ngpus
        else
            print *,    
     &   "Error: unable to determine local process rank"
            stop 1
        end if
# endif
#endif        
        end subroutine initialisation_openacc
