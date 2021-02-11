#include "../include/Mesh.h"

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

const int MeshElemTypeNVerts[7] = {
    [ELEMTYPE_SEG] = 2,
    [ELEMTYPE_TRI] = 3,
    [ELEMTYPE_QUAD] = 4,
    [ELEMTYPE_TETRA] = 4,
    [ELEMTYPE_PYRA] = 5,
    [ELEMTYPE_PRISM] = 6,
    [ELEMTYPE_HEXA] = 8,
};

const int MeshElemTypeNFaces[7] = {
    [ELEMTYPE_SEG] = 1,
    [ELEMTYPE_TRI] = 3,
    [ELEMTYPE_QUAD] = 4,
    [ELEMTYPE_TETRA] = 4,
    [ELEMTYPE_PYRA] = 5,
    [ELEMTYPE_PRISM] = 5,
    [ELEMTYPE_HEXA] = 6,
};

Mesh *Mesh_Create(ReaderInterface *reader, MPI_Comm comm) {
    Mesh *mesh;
    int rank;

    mesh = malloc(sizeof(*mesh));

    mesh->reader = reader;
    mesh->comm = comm;
    mesh->verts = NULL;
    mesh->elems = NULL;

    /* Get rank. */
    MPI_Comm_rank (mesh->comm, &rank);

    /* Process 0 reads CGNS file. */
    if (rank == 0) {
        if (CGNSReader_Read(reader))
            MPI_Abort(mesh->comm, -1);
        ReadCGNSReader(mesh);
    }

    return mesh;
}

void Mesh_Destroy(Mesh *mesh) {
    free(mesh->verts); free(mesh->elems);
    free(mesh->sect_name);

    free(mesh);
}
