#include <mpi.h>

extern "C" {

int MPI_init() { return MPI_Init(0, 0); }

MPI_Op get_mpi_max() { return MPI_MAX; }
MPI_Op get_mpi_sum() { return MPI_SUM; }
MPI_Op get_mpi_min() { return MPI_MIN; }
MPI_Op get_mpi_prod(){ return MPI_PROD; }

MPI_Datatype get_mpi_int() { return MPI_INT; }
MPI_Datatype get_mpi_double() { return MPI_DOUBLE; }
MPI_Datatype get_mpi_char() { return MPI_CHAR; }
MPI_Datatype get_mpi_byte() { return MPI_BYTE; }
MPI_Datatype get_mpi_float() { return MPI_FLOAT; }

MPI_Comm get_mpi_comm_world() { return MPI_COMM_WORLD; }

MPI_Status* get_mpi_status_ignore() { return MPI_STATUS_IGNORE; }
MPI_Info get_mpi_info_null() {return MPI_INFO_NULL; }

}
