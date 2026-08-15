MODULE in_out_manager


! Empty module
   TYPE :: sn_ctl                !: structure for control over output selection
      LOGICAL :: l_runstat = .FALSE.  !: Produce/do not produce run.stat file (T/F)
      LOGICAL :: l_trcstat = .FALSE.  !: Produce/do not produce tracer.stat file (T/F)
      LOGICAL :: l_oceout  = .FALSE.  !: Produce all ocean.outputs    (T) or just one (F)
      LOGICAL :: l_layout  = .FALSE.  !: Produce all layout.dat files (T) or just one (F)
      LOGICAL :: l_prtctl  = .FALSE.  !: Produce/do not produce mpp.output_XXXX files (T/F)
      LOGICAL :: l_prttrc  = .FALSE.  !: Produce/do not produce mpp.top.output_XXXX files (T/F)
      LOGICAL :: l_oasout  = .FALSE.  !: Produce/do not write oasis setup info to ocean.output (T/F)
                                      !  Optional subsetting of processor report files
                                      !  Default settings of 0/1000000/1 should ensure all areas report.
                                      !  Set to a more restrictive range to select specific areas
      INTEGER :: procmin   = 0        !: Minimum narea to output
      INTEGER :: procmax   = 1000000  !: Maximum narea to output
      INTEGER :: procincr  = 1        !: narea increment to output
      INTEGER :: ptimincr  = 1        !: timestep increment to output (time.step and run.stat)
   END TYPE
   TYPE(sn_ctl), SAVE :: sn_cfctl     !: run control structure for selective output, must have SAVE for default init. of sn_ctl


END MODULE in_out_manager
