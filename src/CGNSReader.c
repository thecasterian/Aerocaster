#include "../include/CGNSReader.h"

#include <stdio.h>
#include <stdlib.h>

static int OpenFile(CGNSReader *reader);
static int ReadBase(CGNSReader *reader);
static int ReadZone(CGNSReader *reader);
static int ReadCoord(CGNSReader *reader);
static int ReadSect(CGNSReader *reader);
static int DimElem(ElementType_t elem_type);

CGNSReader *CGNSReader_Create(const char *filename) {
    CGNSReader *reader = malloc(sizeof(*reader));

    /* Copy file name. */
    snprintf(reader->file_name, NAME_MAX_LEN, "%s", filename);

    /* Nullify all pointers for memory safety. */
    reader->x = reader->y = reader->z = NULL;
    reader->sect_name = NULL;
    reader->elem_idx_start = reader->elem_idx_end = reader->nelems = NULL;
    reader->elem_type = NULL;
    reader->elem_conn_size = NULL;
    reader->elem_conn = reader->elem_offset = NULL;

    return reader;
}

int CGNSReader_Read(CGNSReader *reader) {
    if (OpenFile(reader)) return CGNSREADER_ERROR;
    if (ReadBase(reader)) return CGNSREADER_ERROR;
    if (ReadZone(reader)) return CGNSREADER_ERROR;
    if (ReadCoord(reader)) return CGNSREADER_ERROR;
    if (ReadSect(reader)) return CGNSREADER_ERROR;
    printf("read done\n");
    return 0;
}

void CGNSReader_Destroy(CGNSReader *reader) {
    free(reader->x); free(reader->y); free(reader->z);
    free(reader->sect_name);
    free(reader->elem_idx_start); free(reader->elem_idx_end);
    free(reader->nelems);
    free(reader->elem_type);
    free(reader->elem_conn_size);
    free(reader->elem_conn);
    free(reader->elem_offset);

    cg_close(reader->fn);

    free(reader);
}

static int OpenFile(CGNSReader *reader) {
    int file_type;

    /* Check validity. */
    if (cg_is_cgns(reader->file_name, &file_type)) {
        printf("error: invalid file\n");
        return CGNSREADER_ERROR;
    }

    /* Open CGNS file for read-only. */
    if (cg_open(reader->file_name, CG_MODE_READ, &reader->fn)) {
        printf("error: cannot open file\n");
        return CGNSREADER_ERROR;
    }
    printf("read file: %s\n", reader->file_name);

    return 0;
}

static int ReadBase(CGNSReader *reader) {
    int nbases;

    /* Read the number of bases. */
    if (cg_nbases(reader->fn, &nbases))
        return CGNSREADER_ERROR;

    /* Check if only one base exists. */
    if (nbases > 1) {
        printf("error: expected only one base\n");
        return CGNSREADER_ERROR;
    }

    /* Read the base. */
    if (cg_base_read(reader->fn, 1,
                     reader->base_name, &reader->cell_dim, &reader->phys_dim))
        return CGNSREADER_ERROR;
    printf("base \"%s\": cell dimension %d, physical dimension %d\n",
           reader->base_name, reader->cell_dim, reader->phys_dim);

    return 0;
}

static int ReadZone(CGNSReader *reader) {
    int nzones;
    ZoneType_t zone_type;
    int zone_size[3];

    /* Read the number of zones. */
    if (cg_nzones(reader->fn, 1, &nzones))
        return CGNSREADER_ERROR;

    /* Check if only one zone exists. */
    if (nzones > 1) {
        printf("error: expected only one zone\n");
        return CGNSREADER_ERROR;
    }

    /* Read the zone type. */
    if (cg_zone_type(reader->fn, 1, 1, &zone_type))
        return CGNSREADER_ERROR;

    /* Check if the zone type is `Unstructured`. */
    if (zone_type != CGNS_ENUMV(Unstructured)) {
        printf("error: expected unstructured mesh\n");
        return CGNSREADER_ERROR;
    }

    /* Read the zone name and sizes. */
    if (cg_zone_read(reader->fn, 1, 1, reader->zone_name, zone_size))
        return CGNSREADER_ERROR;
    reader->nverts = zone_size[0];
    reader->nelems_internal = zone_size[1];
    printf("zone \"%s\": #vertices %d, #internal elements %d\n",
           reader->zone_name, reader->nverts, reader->nelems_internal);

    return 0;
}

static int ReadCoord(CGNSReader *reader) {
    int rmin, rmax;

    /* Allocate the coordinate arrays. */
    reader->x = malloc(sizeof(*reader->x) * reader->nverts);
    reader->y = malloc(sizeof(*reader->y) * reader->nverts);
    if (reader->cell_dim > 2)
        reader->z = malloc(sizeof(*reader->z) * reader->nverts);

    /* Read the coordinates. Note that the coordinates are able to be read in
       `RealDouble` even if they are actually written in `RealSingle`. */
    rmin = 1;
    rmax = reader->nverts;
    if (cg_coord_read(reader->fn, 1, 1, "CoordinateX",
                      CGNS_ENUMV(RealDouble), &rmin, &rmax, reader->x)) {
        printf("error: reading x coordinates failed\n");
        return CGNSREADER_ERROR;
    }
    if (cg_coord_read(reader->fn, 1, 1, "CoordinateY",
                      CGNS_ENUMV(RealDouble), &rmin, &rmax, reader->y)) {
        printf("error: reading y coordinates failed\n");
        return CGNSREADER_ERROR;
    }
    if (reader->cell_dim > 2)
        if (cg_coord_read(reader->fn, 1, 1, "CoordinateZ",
                        CGNS_ENUMV(RealDouble), &rmin, &rmax, reader->z)) {
            printf("error: reading z coordinates failed\n");
            return CGNSREADER_ERROR;
        }

    return 0;
}

static int ReadSect(CGNSReader *reader) {
    int nbndry, parent_flag, parent_data;
    ElementType_t first_elem_type = CGNS_ENUMV(ElementTypeNull);

    /* Read the nubmer of sections. */
    if (cg_nsections(reader->fn, 1, 1, &reader->nsects))
        return CGNSREADER_ERROR;
    printf("total %d sections\n", reader->nsects);

    /* Allocate the arrays. */
    reader->sect_name = malloc(sizeof(*reader->sect_name) * reader->nsects);
    reader->elem_idx_start = malloc(sizeof(*reader->elem_idx_start) * reader->nsects);
    reader->elem_idx_end = malloc(sizeof(*reader->elem_idx_end) * reader->nsects);
    reader->nelems = malloc(sizeof(*reader->nelems) * reader->nsects);
    reader->elem_type = malloc(sizeof(*reader->elem_type) * reader->nsects);
    reader->elem_conn_size = malloc(sizeof(*reader->elem_conn_size) * reader->nsects);
    reader->elem_conn = malloc(sizeof(*reader->elem_conn) * reader->nsects);
    reader->elem_offset = malloc(sizeof(*reader->elem_offset) * reader->nsects);
    reader->is_internal = malloc(sizeof(*reader->is_internal) * reader->nsects);

    for (int s = 0; s < reader->nsects; s++) {
        /* Read the section info. */
        if (cg_section_read(reader->fn, 1, 1, s+1,
                            reader->sect_name[s], &reader->elem_type[s],
                            &reader->elem_idx_start[s], &reader->elem_idx_end[s],
                            &nbndry, &parent_flag))
            return CGNSREADER_ERROR;
        reader->nelems[s] = reader->elem_idx_end[s] - reader->elem_idx_start[s] + 1;

        /* Read the connectivity data size. */
        if (cg_ElementDataSize(reader->fn, 1, 1, s+1, &reader->elem_conn_size[s]))
            return CGNSREADER_ERROR;

        /* Read the connectivity. */
        reader->elem_conn[s] = malloc(sizeof(*reader->elem_conn[s])
                                      * reader->elem_conn_size[s]);
        switch (reader->elem_type[s]) {
        case CGNS_ENUMV(MIXED):
            /* `MIXED` type required the offset array for output. */
            reader->elem_offset[s] = malloc(sizeof(*reader->elem_offset[s])
                                            * (reader->nelems[s] + 1));
            if (cg_poly_elements_read(reader->fn, 1, 1, s+1,
                                      reader->elem_conn[s], reader->elem_offset[s],
                                      &parent_data))
                return CGNSREADER_ERROR;
            first_elem_type = reader->elem_conn[s][0];
            break;
        case CGNS_ENUMV(NGON_n):
        case CGNS_ENUMV(NFACE_n):
            printf("error: arbitrary polyhedral elements are not supported\n");
            return CGNSREADER_ERROR;
            break;
        default:
            if (cg_elements_read(reader->fn, 1, 1, s+1,
                                 reader->elem_conn[s], NULL))
                return CGNSREADER_ERROR;
            first_elem_type = reader->elem_type[s];
        }

        /* Check if the section contains internal elements or not from the
           dimension of the first element. */
        if ((reader->cell_dim == 2 && DimElem(first_elem_type) == 2)
            || (reader->cell_dim == 3 && DimElem(first_elem_type) == 3))
            reader->is_internal[s] = true;
        else
            reader->is_internal[s] = false;

        /* Elements with different dimensions in one section is not allowed. */
        if (reader->elem_type[s] == CGNS_ENUMV(MIXED))
            for (int i = 1; i < reader->nelems[s]; i++)
                if (DimElem(first_elem_type)
                    != DimElem(reader->elem_conn[s][reader->elem_offset[s][i]])) {
                    printf("error: elements dimensions are different\n");
                    return CGNSREADER_ERROR;
                }

        printf("section \"%s\": type %s, #elements %d, %s\n",
               reader->sect_name[s], ElementTypeName[reader->elem_type[s]],
               reader->nelems[s],
               reader->is_internal[s] ? "internal" : "boundary");
    }

    return 0;
}

static int DimElem(ElementType_t elem_type) {
    switch (elem_type) {
    case CGNS_ENUMV(ElementTypeNull):
    case CGNS_ENUMV(ElementTypeUserDefined):
    case CGNS_ENUMV(MIXED):
    case CGNS_ENUMV(NGON_n):
    case CGNS_ENUMV(NFACE_n):
        return -1;
    case CGNS_ENUMV(NODE):
        return 0;
    case CGNS_ENUMV(BAR_2):
    case CGNS_ENUMV(BAR_3):
    case CGNS_ENUMV(BAR_4):
    case CGNS_ENUMV(BAR_5):
        return 1;
    case CGNS_ENUMV(TRI_3):
    case CGNS_ENUMV(TRI_6):
    case CGNS_ENUMV(TRI_9):
    case CGNS_ENUMV(TRI_10):
    case CGNS_ENUMV(TRI_12):
    case CGNS_ENUMV(TRI_15):
    case CGNS_ENUMV(QUAD_4):
    case CGNS_ENUMV(QUAD_8):
    case CGNS_ENUMV(QUAD_9):
    case CGNS_ENUMV(QUAD_12):
    case CGNS_ENUMV(QUAD_16):
    case CGNS_ENUMV(QUAD_P4_16):
    case CGNS_ENUMV(QUAD_25):
        return 2;
    default:
        return 3;
    }
}
