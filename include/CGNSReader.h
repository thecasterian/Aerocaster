#ifndef AEROCASTER_MESH_READER_H
#define AEROCASTER_MESH_READER_H

#include <cgnslib.h>
#include <stdbool.h>

#include "ReaderInterface.h"

typedef struct _cgns_reader {
    ReaderInterface interface;

    int fn;                             /* CGNS file number. */

    char base_name[SECTNAME_MAX_LEN];       /* Name of base. */
    int cell_dim;                       /* Cell dimension (2-D or 3-D). */
    int phys_dim;                       /* Physical dimension (2-D or 3-D). */

    char zone_name[SECTNAME_MAX_LEN];       /* Name of zone. */
    cgsize_t nverts;                    /* Number of vertices. */
    cgsize_t nelems_internal;           /* Total number of internal (non-boundary) elements. */

    double *x, *y, *z;                  /* Coordinates of vertices. */

    int nsects;                         /* Number of sections. */
    char (*sect_name)[SECTNAME_MAX_LEN];    /* Name of each section. */
    cgsize_t *elem_idx_start;           /* First element index of each section. */
    cgsize_t *elem_idx_end;             /* Last element index of each section. */
    cgsize_t *nelems;                   /* Number of elements in each section. */
    ElementType_t *elem_type;           /* Element type of each section. */
    cgsize_t *elem_conn_size;           /* Size of element connectivity array of each section. */
    cgsize_t **elem_conn;               /* Element connectivity array of each section. */
    cgsize_t **elem_offset;             /* Element offset of each section; ignored for non-mixed sections. */

    bool *is_internal;                  /* Is this setion internal? (T/F) */
} CGNSReader;

ReaderInterface *ReaderInterface_CreateCGNSReader(const char *file_name);

#endif
