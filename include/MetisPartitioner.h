#ifndef METIS_PARTITIONER_H
#define METIS_PARTITIONER_H

#include <stdbool.h>
#include <metis.h>
#include <mpi.h>

#define PARTITIONER_ERROR -1

typedef struct _mesh Mesh;

typedef struct _metis_partitioner {
    Mesh *mesh;                         /* Mesh. */
    bool reorder;                       /* Reorder before partitioning? */

    idx_t nvtxs;                        /* Number of vertices in the graph. */
    idx_t *xadj, *adjncy;               /* Adjacency data. */
    idx_t nparts;                       /* Number of partitions. */
    idx_t objval;                       /* Number of edge cuts. */
    idx_t *part;                        /* Partition vector. */
} MetisPartitioner;

MetisPartitioner *MetisPartitioner_Create(Mesh *mesh, bool reorder);
int MetisPartitioner_PartitionMesh(MetisPartitioner *metis);
void MetisPartitioner_Destroy(MetisPartitioner *metis);

#endif