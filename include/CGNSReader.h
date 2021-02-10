#ifndef AEROCASTER_MESH_READER_H
#define AEROCASTER_MESH_READER_H

#include <cgnslib.h>
#include <stdbool.h>

#define NAME_MAX_LEN 50

#define CGNSREADER_ERROR -1

typedef struct _cgns_reader {
    char file_name[NAME_MAX_LEN];       /* CGNS file name. */
    int fn;                             /* CGNS file number. */

    char base_name[NAME_MAX_LEN];       /* Name of base. */
    int cell_dim;                       /* Cell dimension (2-D or 3-D). */
    int phys_dim;                       /* Physical dimension (2-D or 3-D). */

    char zone_name[NAME_MAX_LEN];       /* Name of zone. */
    int nverts;                         /* Number of vertices. */
    int nelems_internal;                /* Total number of internal (non-boundary) elements. */

    double *x, *y, *z;                  /* Coordinates of vertices. */

    int nsects;                         /* Number of sections. */
    char (*sect_name)[NAME_MAX_LEN];    /* Name of each section. */
    int *elem_idx_start;                /* First element index of each section. */
    int *elem_idx_end;                  /* Last element index of each section. */
    int *nelems;                        /* Number of elements in each section. */
    ElementType_t *elem_type;           /* Element type of each section. */
    int *elem_conn_size;                /* Size of element connectivity array of each section. */
    int **elem_conn;                    /* Element connectivity array of each section. */
    int **elem_offset;                  /* Element offset of each section; ignored for non-mixed sections. */

    bool *is_internal;                  /* Is this setion internal? (T/F) */
} CGNSReader;

CGNSReader *CGNSReader_Create(const char *filename);
int CGNSReader_Read(CGNSReader *reader);
void CGNSReader_Destroy(CGNSReader *reader);

#endif
