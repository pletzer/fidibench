from psyclone.psyir.frontend.fortran import FortranReader
from psyclone.psyir.backend.fortran import FortranWriter
from psyclone.transformations import OMPLoopTrans, OMPParallelTrans
from psyclone.psyir.nodes import Loop

def trans(psyir):
    
    for loop in psyir.walk(Loop):

        # Apply OpenMP transformations
        omp_parallel = OMPParallelTrans()
        omp_loop = OMPLoopTrans()

        omp_parallel.apply(loop)
        omp_loop.apply(loop)

    return psyir