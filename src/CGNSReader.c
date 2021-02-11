#include "../include/CGNSReader.h"

#include <stdio.h>
#include <stdlib.h>

#include "../include/Mesh.h"

#include <glib.h>

typedef struct {
    int n;              /* Number of vertices. */
    size_t idx[4];      /* Indices of vertices. */
} FaceKey;

typedef struct {
    size_t elem_idx[2];
    size_t face_idx[2];
} FaceValue;

static const int ElemTypeDim[7] = {
    [ELEMTYPE_SEG] = 1,
    [ELEMTYPE_TRI] = 2,
    [ELEMTYPE_QUAD] = 2,
    [ELEMTYPE_TETRA] = 3,
    [ELEMTYPE_PYRA] = 3,
    [ELEMTYPE_PRISM] = 3,
    [ELEMTYPE_HEXA] = 3,
};

static int CGNSReader_Read(void *_reader, const char *file_name);
static void CGNSReader_WriteToMesh(void *_reader, Mesh *mesh);
static void CGNSReader_Destroy(void *_reader);

static int OpenFile(CGNSReader *reader, const char *file_name);
static int ReadBase(CGNSReader *reader);
static int ReadZone(CGNSReader *reader);
static int ReadCoord(CGNSReader *reader);
static int ReadSect(CGNSReader *reader);
static int DimElem(ElementType_t elem_type);

static void ReadInternal(CGNSReader *reader, Mesh *mesh, int s);
static GTree *GetAdjacency(Mesh *mesh);
static void ReadBoundary(CGNSReader *reader, Mesh *mesh,
                         GTree *face_tree, int s);

static MeshElemType CGNSToElemType(ElementType_t type);
static void StoreToFaceTree(GTree *tree, int n, ...);

static gint CompareFace(gconstpointer a, gconstpointer b,
                        gpointer user_data G_GNUC_UNUSED);
static gboolean FindAdjElem(gpointer key G_GNUC_UNUSED, gpointer value,
                            gpointer data);

ReaderInterface *CGNSReader_Create(void) {
    CGNSReader *reader = malloc(sizeof(*reader));

    /* Initialize the interface. */
    reader->interface.reader = reader;
    reader->interface.read_file = CGNSReader_Read;
    reader->interface.write_to_mesh = CGNSReader_WriteToMesh;
    reader->interface.destroy = CGNSReader_Destroy;

    /* Nullify all pointers for memory safety. */
    reader->x = reader->y = reader->z = NULL;
    reader->sect_name = NULL;
    reader->elem_idx_start = reader->elem_idx_end = reader->nelems = NULL;
    reader->elem_type = NULL;
    reader->elem_conn_size = NULL;
    reader->elem_conn = reader->elem_offset = NULL;

    return &reader->interface;
}

static int CGNSReader_Read(void *_reader, const char *file_name) {
    CGNSReader *const reader = _reader;

    if (OpenFile(reader, file_name))
        return READER_ERROR;
    if (ReadBase(reader))
        return READER_ERROR;
    if (ReadZone(reader))
        return READER_ERROR;
    if (ReadCoord(reader))
        return READER_ERROR;
    if (ReadSect(reader))
        return READER_ERROR;
    printf("read done\n");

    return 0;
}

static int OpenFile(CGNSReader *reader, const char *file_name) {
    int file_type;

    /* Check validity. */
    if (cg_is_cgns(file_name, &file_type)) {
        printf("error: invalid file\n");
        return READER_ERROR;
    }

    /* Open CGNS file for read-only. */
    if (cg_open(file_name, CG_MODE_READ, &reader->fn)) {
        printf("error: cannot open file\n");
        return READER_ERROR;
    }
    printf("read file: %s\n", file_name);

    return 0;
}

static int ReadBase(CGNSReader *reader) {
    int nbases;

    /* Read the number of bases. */
    if (cg_nbases(reader->fn, &nbases))
        return READER_ERROR;

    /* Check if only one base exists. */
    if (nbases > 1) {
        printf("error: expected only one base\n");
        return READER_ERROR;
    }

    /* Read the base. */
    if (cg_base_read(reader->fn, 1,
                     reader->base_name, &reader->cell_dim, &reader->phys_dim))
        return READER_ERROR;
    printf("base \"%s\": cell dimension %d, physical dimension %d\n",
           reader->base_name, reader->cell_dim, reader->phys_dim);

    return 0;
}

static int ReadZone(CGNSReader *reader) {
    int nzones;
    ZoneType_t zone_type;
    cgsize_t zone_size[3];

    /* Read the number of zones. */
    if (cg_nzones(reader->fn, 1, &nzones))
        return READER_ERROR;

    /* Check if only one zone exists. */
    if (nzones > 1) {
        printf("error: expected only one zone\n");
        return READER_ERROR;
    }

    /* Read the zone type. */
    if (cg_zone_type(reader->fn, 1, 1, &zone_type))
        return READER_ERROR;

    /* Check if the zone type is `Unstructured`. */
    if (zone_type != CGNS_ENUMV(Unstructured)) {
        printf("error: expected unstructured mesh\n");
        return READER_ERROR;
    }

    /* Read the zone name and sizes. */
    if (cg_zone_read(reader->fn, 1, 1, reader->zone_name, zone_size))
        return READER_ERROR;
    reader->nverts = zone_size[0];
    reader->nelems_internal = zone_size[1];
    printf("zone \"%s\": #vertices %ld, #internal elements %ld\n",
           reader->zone_name,
           (long)reader->nverts, (long)reader->nelems_internal);

    return 0;
}

static int ReadCoord(CGNSReader *reader) {
    cgsize_t rmin, rmax;

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
        return READER_ERROR;
    }
    if (cg_coord_read(reader->fn, 1, 1, "CoordinateY",
                      CGNS_ENUMV(RealDouble), &rmin, &rmax, reader->y)) {
        printf("error: reading y coordinates failed\n");
        return READER_ERROR;
    }
    if (reader->cell_dim > 2)
        if (cg_coord_read(reader->fn, 1, 1, "CoordinateZ",
                        CGNS_ENUMV(RealDouble), &rmin, &rmax, reader->z)) {
            printf("error: reading z coordinates failed\n");
            return READER_ERROR;
        }

    return 0;
}

static int ReadSect(CGNSReader *reader) {
    int nbndry, parent_flag;
    cgsize_t parent_data;
    ElementType_t first_elem_type = CGNS_ENUMV(ElementTypeNull);

    /* Read the nubmer of sections. */
    if (cg_nsections(reader->fn, 1, 1, &reader->nsects))
        return READER_ERROR;
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
            return READER_ERROR;
        reader->nelems[s] = reader->elem_idx_end[s] - reader->elem_idx_start[s] + 1;

        /* Read the connectivity data size. */
        if (cg_ElementDataSize(reader->fn, 1, 1, s+1, &reader->elem_conn_size[s]))
            return READER_ERROR;

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
                return READER_ERROR;
            first_elem_type = reader->elem_conn[s][0];
            break;
        case CGNS_ENUMV(NGON_n):
        case CGNS_ENUMV(NFACE_n):
            printf("error: arbitrary polyhedral elements are not supported\n");
            return READER_ERROR;
            break;
        default:
            if (cg_elements_read(reader->fn, 1, 1, s+1,
                                 reader->elem_conn[s], NULL))
                return READER_ERROR;
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
                    return READER_ERROR;
                }

        printf("section \"%s\": type %s, #elements %ld, %s\n",
               reader->sect_name[s], ElementTypeName[reader->elem_type[s]],
               (long)reader->nelems[s],
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

static void CGNSReader_WriteToMesh(void *_reader, Mesh *mesh) {
    CGNSReader *const reader = _reader;
    GTree *face_tree;

    /* Read metadata. */
    mesh->dim = reader->cell_dim;
    mesh->nverts = reader->nverts;
    mesh->nelems = 0;
    mesh->nsects = reader->nsects;

    /* Allocate arrays. */
    mesh->verts = malloc(sizeof(*mesh->verts) * reader->nverts);
    mesh->elems = malloc(sizeof(*mesh->elems) * reader->nelems_internal);
    mesh->sect_name = malloc(sizeof(*mesh->sect_name) * reader->nsects);

    /* Initialize some member variables. */
    for (int i = 0; i < reader->nelems_internal; i++)
        for (int j = 0; j < 6; j++)
            mesh->elems[i].face_section[j] = FACE_SECT_UNSPEC_INTER;

    /* Read coordinates and section names. */
    for (int i = 0; i < reader->nverts; i++) {
        mesh->verts[i].x = reader->x[i];
        mesh->verts[i].y = reader->y[i];
        if (mesh->dim == 3)
            mesh->verts[i].z = reader->z[i];
    }
    memcpy(mesh->sect_name, reader->sect_name,
           sizeof(*mesh->sect_name) * reader->nsects);

    /* Read internal elements. */
    for (int s = 0; s < reader->nsects; s++)
        if (reader->is_internal[s])
            ReadInternal(reader, mesh, s);

    /* Calculate the adjacency info between elements. */
    face_tree = GetAdjacency(mesh);

    /* Read boundary elements. */
    for (int s = 0; s < reader->nsects; s++)
        if (!reader->is_internal[s])
            ReadBoundary(reader, mesh, face_tree, s);

    g_tree_destroy(face_tree);

    for (int i = 0; i < reader->nelems_internal; i++) {
        for (int j = 0; j < MeshElemTypeNFaces[mesh->elems[i].type]; j++) {
            if (mesh->elems[i].face_section[j] == FACE_SECT_UNSPEC_INTER
                && mesh->elems[i].idx_adj[j] == IDX_ADJ_NO_ADJ)
                mesh->elems[i].face_section[j] = FACE_SECT_UNSPEC_BNDRY;
        }
    }
}

static void ReadInternal(CGNSReader *reader, Mesh *mesh, int s) {
    cgsize_t *start;
    int npe;

    if (reader->elem_type[s] == CGNS_ENUMV(MIXED)) {
        /* In the element type of a section is `MIXED`, the connectivity of each
           element is represented as:
               ElementType VertexIdx1 VertexIdx2 ... VertexIdxN.
           The index of ElementType in the connectivity array is stored in the
           offset array. */

        for (int i = 0; i < reader->nelems[s]; i++) {
            /* Pointer to ElementType. */
            start = &(reader->elem_conn[s][reader->elem_offset[s][i]]);

            /* Read the type as MeshElemType. */
            mesh->elems[mesh->nelems].type = CGNSToElemType(start[0]);
            npe = MeshElemTypeNVerts[mesh->elems[mesh->nelems].type];

            /* Read the indices of vertices. Additional nodes in a quadratic,
               cubic, or quartic element are ignored. */
            for (int j = 0; j < npe; j++)
                mesh->elems[mesh->nelems].idx_verts[j] = start[j+1];

            /* Set its section. */
            mesh->elems[mesh->nelems].section = s;

            mesh->nelems++;
        }
    } else {
        /* For the other element types, the connectivity array contains the
           indices of vertices only. */

        /* Get the number of vertices in each element. */
        cg_npe(reader->elem_type[s], &npe);

        for (int i = 0; i < reader->nelems[s]; i++) {
            /* Pointer to the index of the first vertex. */
            start = &(reader->elem_conn[s][npe * i]);

            /* All elements have same type. */
            mesh->elems[mesh->nelems].type = CGNSToElemType(reader->elem_type[s]);

            /* Read the indices of vertices. */
            for (int j = 0; j < npe; j++)
                mesh->elems[mesh->nelems].idx_verts[j] = start[j];

            /* Set its section. */
            mesh->elems[mesh->nelems].section = s;

            mesh->nelems++;
        }
    }
}

static GTree *GetAdjacency(Mesh *mesh) {
    GTree *face_tree;
    size_t *v;

    face_tree = g_tree_new_full(CompareFace, NULL, g_free, g_free);

    /* Build the face tree. */
    for (size_t i = 0; i < mesh->nelems; i++) {
        v = mesh->elems[i].idx_verts;

        switch (mesh->elems[i].type) {
        case ELEMTYPE_TRI:
            StoreToFaceTree(face_tree, 2, v[0], v[1], i, 0);
            StoreToFaceTree(face_tree, 2, v[1], v[2], i, 1);
            StoreToFaceTree(face_tree, 2, v[2], v[0], i, 2);
            break;
        case ELEMTYPE_QUAD:
            StoreToFaceTree(face_tree, 2, v[0], v[1], i, 0);
            StoreToFaceTree(face_tree, 2, v[1], v[2], i, 1);
            StoreToFaceTree(face_tree, 2, v[2], v[3], i, 2);
            StoreToFaceTree(face_tree, 2, v[3], v[0], i, 3);
            break;
        case ELEMTYPE_TETRA:
            StoreToFaceTree(face_tree, 3, v[0], v[2], v[1], i, 0);
            StoreToFaceTree(face_tree, 3, v[0], v[1], v[3], i, 1);
            StoreToFaceTree(face_tree, 3, v[1], v[2], v[3], i, 2);
            StoreToFaceTree(face_tree, 3, v[2], v[0], v[3], i, 3);
            break;
        case ELEMTYPE_PYRA:
            StoreToFaceTree(face_tree, 4, v[0], v[3], v[2], v[1], i, 0);
            StoreToFaceTree(face_tree, 3, v[0], v[1], v[4],       i, 1);
            StoreToFaceTree(face_tree, 3, v[1], v[2], v[4],       i, 2);
            StoreToFaceTree(face_tree, 3, v[2], v[3], v[4],       i, 3);
            StoreToFaceTree(face_tree, 3, v[3], v[0], v[4],       i, 4);
            break;
        case ELEMTYPE_PRISM:
            StoreToFaceTree(face_tree, 4, v[0], v[1], v[4], v[3], i, 0);
            StoreToFaceTree(face_tree, 4, v[1], v[2], v[5], v[4], i, 1);
            StoreToFaceTree(face_tree, 4, v[2], v[0], v[3], v[5], i, 2);
            StoreToFaceTree(face_tree, 3, v[0], v[2], v[1],       i, 3);
            StoreToFaceTree(face_tree, 3, v[3], v[4], v[5],       i, 4);
            break;
        case ELEMTYPE_HEXA:
            StoreToFaceTree(face_tree, 4, v[0], v[3], v[2], v[1], i, 0);
            StoreToFaceTree(face_tree, 4, v[0], v[1], v[5], v[4], i, 1);
            StoreToFaceTree(face_tree, 4, v[1], v[2], v[6], v[5], i, 2);
            StoreToFaceTree(face_tree, 4, v[2], v[3], v[7], v[6], i, 3);
            StoreToFaceTree(face_tree, 4, v[0], v[4], v[7], v[3], i, 4);
            StoreToFaceTree(face_tree, 4, v[4], v[5], v[6], v[7], i, 5);
            break;
        default:
            cg_error_exit();
        }
    }

    g_tree_foreach(face_tree, FindAdjElem, mesh);

    return face_tree;
}

static void ReadBoundary(CGNSReader *reader, Mesh *mesh,
                         GTree *face_tree, int s) {
    cgsize_t *start, tmp;
    MeshElemType type;
    FaceKey key;
    FaceValue *value;

    if (reader->elem_type[s] == CGNS_ENUMV(MIXED)) {
        for (int i = 0; i < reader->nelems[s]; i++) {
            /* Pointer to ElementType. */
            start = &(reader->elem_conn[s][reader->elem_offset[s][i]]);

            /* Read the type as MeshElemType. */
            type = CGNSToElemType(start[0]);
            key.n = MeshElemTypeNVerts[type];

            /* The dimension of a boundary element must be 1-D if the mesh is
               2-D or 2-D if the mesh is 3-D. Otherwise, ignore it. */
            if ((mesh->dim == 2 && ElemTypeDim[type] != 1)
                || (mesh->dim == 3 && ElemTypeDim[type] != 2))
                continue;

            /* Read the indices of vertices. */
            for (int j = 0; j < key.n; j++)
                key.idx[j] = start[j+1];

            /* Sort the indices using the bubble sort. */
            for (int j = 0; j < key.n; j++)
                for (int k = 0; k < key.n-1; k++)
                    if (key.idx[k] > key.idx[k+1]) {
                        tmp = key.idx[k];
                        key.idx[k] = key.idx[k+1];
                        key.idx[k+1] = tmp;
                    }

            /* Find the internal elements containing this boundary element. */
            value = g_tree_lookup(face_tree, &key);
            if (!value) cg_error_exit();

            /* Set the face section of the element. */
            mesh->elems[value->elem_idx[0]].face_section[value->face_idx[0]] = s;
            if (value->elem_idx[1] != (size_t)(-1))
                mesh->elems[value->elem_idx[1]].face_section[value->face_idx[1]] = s;
        }
    }

    else {
        /* Read the type as MeshElemType. */
        type = CGNSToElemType(reader->elem_type[s]);
        key.n = MeshElemTypeNVerts[type];

        /* The dimension of a boundary element must be 1-D if the mesh is
           2-D or 2-D if the mesh is 3-D. Otherwise, ignore it. */
        if ((mesh->dim == 2 && ElemTypeDim[type] != 1)
            || (mesh->dim == 3 && ElemTypeDim[type] != 2))
            return;

        /* Get the number of vertices in each element. */
        cg_npe(reader->elem_type[s], &key.n);

        for (int i = 0; i < reader->nelems[s]; i++) {
            /* Pointer to the index of the first vertex. */
            start = &(reader->elem_conn[s][key.n * i]);

            /* Read the indices of vertices. */
            for (int j = 0; j < key.n; j++)
                key.idx[j] = start[j];

            /* Sort the indices using the bubble sort. */
            for (int j = 0; j < key.n; j++)
                for (int k = 0; k < key.n-1; k++)
                    if (key.idx[k] > key.idx[k+1]) {
                        tmp = key.idx[k];
                        key.idx[k] = key.idx[k+1];
                        key.idx[k+1] = tmp;
                    }

            /* Find the internal elements containing this boundary element. */
            value = g_tree_lookup(face_tree, &key);
            if (!value) cg_error_exit();

            /* Set the face section of the element. */
            mesh->elems[value->elem_idx[0]].face_section[value->face_idx[0]] = s;
            if (value->elem_idx[1] != (size_t)(-1))
                mesh->elems[value->elem_idx[1]].face_section[value->face_idx[1]] = s;
        }
    }
}

static MeshElemType CGNSToElemType(ElementType_t type) {
    switch (type) {
    case CGNS_ENUMV(BAR_2):
    case CGNS_ENUMV(BAR_3):
    case CGNS_ENUMV(BAR_4):
    case CGNS_ENUMV(BAR_5):
        return ELEMTYPE_SEG;
    case CGNS_ENUMV(TRI_3):
    case CGNS_ENUMV(TRI_6):
    case CGNS_ENUMV(TRI_9):
    case CGNS_ENUMV(TRI_10):
    case CGNS_ENUMV(TRI_12):
    case CGNS_ENUMV(TRI_15):
        return ELEMTYPE_TRI;
    case CGNS_ENUMV(QUAD_4):
    case CGNS_ENUMV(QUAD_8):
    case CGNS_ENUMV(QUAD_9):
    case CGNS_ENUMV(QUAD_12):
    case CGNS_ENUMV(QUAD_16):
    case CGNS_ENUMV(QUAD_P4_16):
    case CGNS_ENUMV(QUAD_25):
        return ELEMTYPE_QUAD;
    case CGNS_ENUMV(TETRA_4):
    case CGNS_ENUMV(TETRA_10):
    case CGNS_ENUMV(TETRA_16):
    case CGNS_ENUMV(TETRA_20):
    case CGNS_ENUMV(TETRA_22):
    case CGNS_ENUMV(TETRA_34):
    case CGNS_ENUMV(TETRA_35):
        return ELEMTYPE_TETRA;
    case CGNS_ENUMV(PYRA_5):
    case CGNS_ENUMV(PYRA_13):
    case CGNS_ENUMV(PYRA_14):
    case CGNS_ENUMV(PYRA_21):
    case CGNS_ENUMV(PYRA_29):
    case CGNS_ENUMV(PYRA_30):
    case CGNS_ENUMV(PYRA_P4_29):
    case CGNS_ENUMV(PYRA_50):
    case CGNS_ENUMV(PYRA_55):
        return ELEMTYPE_PYRA;
    case CGNS_ENUMV(PENTA_6):
    case CGNS_ENUMV(PENTA_15):
    case CGNS_ENUMV(PENTA_18):
    case CGNS_ENUMV(PENTA_24):
    case CGNS_ENUMV(PENTA_38):
    case CGNS_ENUMV(PENTA_40):
    case CGNS_ENUMV(PENTA_33):
    case CGNS_ENUMV(PENTA_66):
    case CGNS_ENUMV(PENTA_75):
        return ELEMTYPE_PRISM;
    case CGNS_ENUMV(HEXA_8):
    case CGNS_ENUMV(HEXA_20):
    case CGNS_ENUMV(HEXA_27):
    case CGNS_ENUMV(HEXA_32):
    case CGNS_ENUMV(HEXA_56):
    case CGNS_ENUMV(HEXA_64):
    case CGNS_ENUMV(HEXA_44):
    case CGNS_ENUMV(HEXA_98):
    case CGNS_ENUMV(HEXA_125):
        return ELEMTYPE_HEXA;
    default:
        cg_error_exit();
        return -1;
    }
}

static void StoreToFaceTree(GTree *tree, int n, ...) {
    va_list ap;
    FaceKey key_static, *key;
    FaceValue *value;
    size_t tmp;

    va_start(ap, n);

    key_static.n = n;
    for (int i = 0; i < n; i++) key_static.idx[i] = va_arg(ap, size_t);

    /* Sort the indices using the bubble sort. */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n-1; j++)
            if (key_static.idx[j] > key_static.idx[j+1]) {
                tmp = key_static.idx[j];
                key_static.idx[j] = key_static.idx[j+1];
                key_static.idx[j+1] = tmp;
            }

    value = g_tree_lookup(tree, &key_static);

    if (value) {
        if (value->elem_idx[1] != (size_t)(-1)) cg_error_exit();

        value->elem_idx[1] = va_arg(ap, int);
        value->face_idx[1] = va_arg(ap, int);
    }
    else {
        key = malloc(sizeof(*key));
        value = malloc(sizeof(*value));

        memcpy(key, &key_static, sizeof(key_static));
        value->elem_idx[0] = va_arg(ap, int);
        value->face_idx[0] = va_arg(ap, int);
        value->elem_idx[1] = value->face_idx[1] = (size_t)(-1);

        g_tree_insert(tree, key, value);
    }

    va_end(ap);
}

static gint CompareFace(gconstpointer a, gconstpointer b,
                        gpointer user_data G_GNUC_UNUSED) {
    FaceKey *fa, *fb;

    fa = (FaceKey *)a;
    fb = (FaceKey *)b;

    if (fa->n != fb->n) return fa->n - fb->n;
    for (int i = 0; i < fa->n; i++) {
        if (fa->idx[i] < fb->idx[i]) return -1;
        else if (fa->idx[i] > fb->idx[i]) return 1;
    }
    return 0;
}

static gboolean FindAdjElem(gpointer key G_GNUC_UNUSED, gpointer value,
                            gpointer data) {
    FaceValue *fv;
    Mesh *mesh;

    fv = (FaceValue *)value;
    mesh = (Mesh *)data;

    mesh->elems[fv->elem_idx[0]].idx_adj[fv->face_idx[0]]
        = fv->elem_idx[1] != (size_t)(-1) ? fv->elem_idx[1] : IDX_ADJ_NO_ADJ;
    if (fv->elem_idx[1] != (size_t)(-1))
        mesh->elems[fv->elem_idx[1]].idx_adj[fv->face_idx[1]] = fv->elem_idx[0];

    return false;
}

static void CGNSReader_Destroy(void *_reader) {
    CGNSReader *const reader = _reader;

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
