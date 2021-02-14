#include "../include/MetisPartitioner.h"

#include <stdio.h>
#include <stdlib.h>

#include "../include/Mesh.h"

static void CalcAdjncy(MetisPartitioner *metis);
static int ReorderMesh(MetisPartitioner *metis);
static int CalcPartition(MetisPartitioner *metis);
static void Distribute(MetisPartitioner *metis);

MetisPartitioner *MetisPartitioner_Create(Mesh *mesh, bool reorder) {
    MetisPartitioner *metis;

    metis = malloc(sizeof(*metis));
    metis->mesh = mesh;
    metis->reorder = reorder;
    metis->xadj = NULL;
    metis->adjncy = NULL;
    metis->part = NULL;

    return metis;
}

int MetisPartitioner_PartitionMesh(MetisPartitioner *metis) {
    int comm_size, rank;
    idx_t adjncy_size;

    /* Get the number of processes in the communicator, i.e. the number of
       partitions. */
    MPI_Comm_size(metis->mesh->comm, &comm_size);
    if (comm_size == 1 && !metis->reorder) return 0;

    /* Get the rank. */
    MPI_Comm_rank(metis->mesh->comm, &rank);

    if (rank == 0) {
        /* Calculate the size of adjncy array. */
        adjncy_size = Mesh_AdjacencySize(metis->mesh);

        /* Initialize the partitioner of process 0. */
        metis->nvtxs = metis->mesh->nelems;
        metis->xadj = malloc((metis->mesh->nelems+1) * sizeof(idx_t));
        metis->adjncy = malloc(adjncy_size * sizeof(idx_t));
        metis->nparts = comm_size;
        metis->part = malloc(metis->mesh->nelems * sizeof(idx_t));

        /* Calculate the adjacency of the mesh. */
        CalcAdjncy(metis);

        /* If `reorder` is true, reorder the mesh and calculate the adjacency
           again. */
        if (metis->reorder) {
            if (ReorderMesh(metis))
                return PARTITIONER_ERROR;
            CalcAdjncy(metis);
        }

        /* Calculate the partition. */
        if (comm_size != 1)
            if (CalcPartition(metis))
                return PARTITIONER_ERROR;
    }

    Distribute(metis);

    return 0;
}

void MetisPartitioner_Destroy(MetisPartitioner *metis) {
    free(metis->xadj);
    free(metis->adjncy);
    free(metis->part);

    free(metis);
}

static void CalcAdjncy(MetisPartitioner *metis) {
    Mesh *const mesh = metis->mesh;
    idx_t cnt;

    metis->xadj[0] = 0;
    cnt = 0;
    for (size_t i = 0; i < mesh->nelems; i++) {
        for (int j = 0; j < MeshElemTypeNFaces[mesh->elems[i].type]; j++)
            if (mesh->elems[i].idx_adj[j] != IDX_ADJ_NO_ADJ)
                metis->adjncy[cnt++] = mesh->elems[i].idx_adj[j];
        metis->xadj[i+1] = cnt;
    }
}

static int ReorderMesh(MetisPartitioner *metis) {
    idx_t *perm, *iperm;
    idx_t options[METIS_NOPTIONS];
    int err;

    perm = malloc(metis->mesh->nelems * sizeof(idx_t));
    iperm = malloc(metis->mesh->nelems * sizeof(idx_t));

    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;

    err = METIS_NodeND(&metis->nvtxs,   /* nvtxs */
                       metis->xadj,     /* xadj */
                       metis->adjncy,   /* adjncy */
                       NULL,            /* vwgt */
                       options,         /* options */
                       perm,            /* perm */
                       iperm            /* iperm */
                       );
    if (err != METIS_OK) {
        printf("error: node reordering failed\n");
        return PARTITIONER_ERROR;
    }

    // TODO: permute the elements in-place.

    return 0;
}

static int CalcPartition(MetisPartitioner *metis) {
    idx_t ncon = 1;
    idx_t options[METIS_NOPTIONS];
    int err;

    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0;

    err = METIS_PartGraphKway(&metis->nvtxs,    /* nvtxs */
                              &ncon,            /* ncon */
                              metis->xadj,      /* xadj */
                              metis->adjncy,    /* adjncy */
                              NULL,             /* vwgt */
                              NULL,             /* vsize */
                              NULL,             /* adjwgt */
                              &metis->nparts,   /* nparts */
                              NULL,             /* tpwgts */
                              NULL,             /* ubvec */
                              options,          /* options */
                              &metis->objval,   /* objval */
                              metis->part       /* part */
                              );
    if (err != METIS_OK) {
        printf("error: mesh partitioning failed\n");
        return PARTITIONER_ERROR;
    }

    return 0;
}

static void Distribute(MetisPartitioner *metis) {
    // TODO: distribute mesh based on the calculated partition.
}
