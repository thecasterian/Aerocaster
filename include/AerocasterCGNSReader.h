#ifndef AEROCASTER_MESH_READER_H
#define AEROCASTER_MESH_READER_H

#include <cgnslib.h>
#include <stdbool.h>

#define CGNS_STR_LEN 50

typedef struct {
    int fn;                             /* CGNS file number. */

    char base_name[CGNS_STR_LEN];       /* Name of base. */
    int cell_dim;                       /* Cell dimension (2-D or 3-D). */
    int phys_dim;                       /* Physical dimension (2-D or 3-D). */

    char zone_name[CGNS_STR_LEN];       /* Name of zone. */
    int nverts;                         /* Number of vertices. */
    int nelems_internal;                /* Total number of internal (non-boundary) elements. */

    double *x, *y, *z;                  /* Coordinates of vertices. */

    int nsects;                         /* Number of sections. */
    char (*sect_name)[CGNS_STR_LEN];    /* Name of each section. */
    int *elem_idx_start;                /* First element index of each section. */
    int *elem_idx_end;                  /* Last element index of each section. */
    int *nelems;                        /* Number of elements in each section. */
    ElementType_t *elem_type;           /* Element type of each section. */
    int *elem_conn_size;                /* Size of element connectivity array of each section. */
    int **elem_conn;                    /* Element connectivity array of each section. */
    int **elem_offset;                  /* Element offset of each section; ignored for non-mixed sections. */

    bool *is_internal;                  /* Is this setion internal? (T/F) */
} AerocasterCGNSMeshReader;

AerocasterCGNSMeshReader *AerocasterCGNSMeshReader_Create(void);
void AerocasterCGNSMeshReader_Read(AerocasterCGNSMeshReader *reader, const char *file_name);
void AerocasterCGNSMeshReader_Destroy(AerocasterCGNSMeshReader *reader);

#endif
