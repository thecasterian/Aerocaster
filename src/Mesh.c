#include "../include/Mesh.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "../include/MetisPartitioner.h"

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
    MetisPartitioner *metis;
    int rank;

    mesh = malloc(sizeof(*mesh));

    mesh->comm = comm;
    mesh->verts = NULL;
    mesh->elems = NULL;
    mesh->sect_name = NULL;

    /* Get rank. */
    MPI_Comm_rank(mesh->comm, &rank);

    /* Process 0 reads CGNS file. */
    if (rank == 0)
        if (ReaderInterface_ReadFile(reader)
            || ReaderInterface_WriteToMesh(reader, mesh))
            MPI_Abort(mesh->comm, -1);

    metis = MetisPartitioner_Create(mesh, true);
    if (MetisPartitioner_PartitionMesh(metis))
        MPI_Abort(mesh->comm, -1);
    MetisPartitioner_Destroy(metis);

    return mesh;
}

size_t Mesh_AdjacencySize(Mesh *mesh) {
    size_t adjacency_size;

    adjacency_size = 0;
    for (size_t i = 0; i < mesh->nelems; i++)
        for (int j = 0; j < MeshElemTypeNFaces[mesh->elems[i].type]; j++)
            if (mesh->elems[i].idx_adj[j] != IDX_ADJ_NO_ADJ)
                adjacency_size++;

    return adjacency_size;
}

void Mesh_Destroy(Mesh *mesh) {
    free(mesh->verts);
    free(mesh->elems);
    free(mesh->sect_name);

    free(mesh);
}
