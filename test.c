#include <stdio.h>

#include "include/CGNSReader.h"
#include "include/Mesh.h"

#include <mpi.h>
#include <parmetis.h>

int main(int argc, char *argv[]) {
    ReaderInterface *reader;
    Mesh *mesh;

    MPI_Init(&argc, &argv);

    reader = CGNSReader_Create();

    mesh = Mesh_Create(reader, MPI_COMM_WORLD);

    Mesh_Destroy(mesh);
    CGNSReader_Destroy(reader);

    MPI_Finalize();

    return 0;
}
