! Created by albpl on 3/31/2025.

program parallel_example
    USE OMP_LIB


    INTEGER :: thread_id

    !$OMP PARALLEL

    PRINT *, "Hello from process: ", OMP_GET_THREAD_NUM()

    !$OMP END PARALLEL
end program parallel_example