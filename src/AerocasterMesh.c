#include "../include/AerocasterMesh.h"

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include <cgnslib.h>
#include <glib.h>

typedef struct {
    int n;              /* Number of vertices. */
    int idx[4];         /* Indices of vertices. */
} FaceKey;

typedef struct {
    int elem_idx[2];
    int face_idx[2];
} FaceValue;

static void ReadInternal(AerocasterMesh *mesh, AerocasterCGNSMeshReader *reader, int s);
static GTree *GetAdjacency(AerocasterMesh *mesh);
static void ReadBoundary(AerocasterMesh *mesh, AerocasterCGNSMeshReader *reader, GTree *face_tree, int s);

static AerocasterMeshElementType CGNSToAerocasterType(ElementType_t type);
static void StoreToFaceTree(GTree *tree, int n, ...);

static gint CompareFace(gconstpointer a, gconstpointer b,
                        gpointer user_data G_GNUC_UNUSED);
static gboolean FindAdjElem(gpointer key G_GNUC_UNUSED, gpointer value,
                            gpointer data);

const int AerocasterMeshElementTypeNVerts[7] = {
    [AEROCASTER_SEG] = 2,
    [AEROCASTER_TRI] = 3,
    [AEROCASTER_QUAD] = 4,
    [AEROCASTER_TETRA] = 4,
    [AEROCASTER_PYRA] = 5,
    [AEROCASTER_PRISM] = 6,
    [AEROCASTER_HEXA] = 8,
};

const int AerocasterMeshElementTypeNFaces[7] = {
    [AEROCASTER_SEG] = 1,
    [AEROCASTER_TRI] = 3,
    [AEROCASTER_QUAD] = 4,
    [AEROCASTER_TETRA] = 4,
    [AEROCASTER_PYRA] = 5,
    [AEROCASTER_PRISM] = 5,
    [AEROCASTER_HEXA] = 6,
};

static const int ElemTypeDim[7] = {
    [AEROCASTER_SEG] = 1,
    [AEROCASTER_TRI] = 2,
    [AEROCASTER_QUAD] = 2,
    [AEROCASTER_TETRA] = 3,
    [AEROCASTER_PYRA] = 3,
    [AEROCASTER_PRISM] = 3,
    [AEROCASTER_HEXA] = 3,
};

AerocasterMesh *AerocasterMesh_Create(void) {
    AerocasterMesh *mesh = malloc(sizeof(*mesh));

    mesh->verts = NULL;
    mesh->elems = NULL;

    return mesh;
}

void AerocasterMesh_ReadCGNSMeshReader(AerocasterMesh *mesh,
                                       AerocasterCGNSMeshReader *reader) {
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
            mesh->elems[i].face_section[j] = AEROCASTER_UNSPEC_INTER;

    /* Read coordinates and section names. */
    for (int i = 0; i < reader->nverts; i++) {
        mesh->verts[i].x = reader->x[i];
        mesh->verts[i].y = reader->y[i];
        if (mesh->dim == 3)
            mesh->verts[i].z = reader->z[i];
    }
    memcpy(mesh->sect_name, reader->sect_name, sizeof(*mesh->sect_name) * reader->nsects);

    /* Read internal elements. */
    for (int s = 0; s < reader->nsects; s++)
        if (reader->is_internal[s])
            ReadInternal(mesh, reader, s);

    /* Calculate the adjacency info between elements. */
    face_tree = GetAdjacency(mesh);

    /* Read boundary elements. */
    for (int s = 0; s < reader->nsects; s++)
        if (!reader->is_internal[s])
            ReadBoundary(mesh, reader, face_tree, s);

    g_tree_destroy(face_tree);

    for (int i = 0; i < reader->nelems_internal; i++) {
        for (int j = 0; j < AerocasterMeshElementTypeNFaces[mesh->elems[i].type]; j++) {
            if (mesh->elems[i].face_section[j] == AEROCASTER_UNSPEC_INTER
                && mesh->elems[i].idx_adj[j] == AEROCASTER_NO_ADJ)
                mesh->elems[i].face_section[j] = AEROCASTER_UNSPEC_BNDRY;
        }
    }
}

void AerocasterMesh_Destroy(AerocasterMesh *mesh) {
    free(mesh->verts); free(mesh->elems);
    free(mesh->sect_name);

    free(mesh);
}

static void ReadInternal(AerocasterMesh *mesh,
                         AerocasterCGNSMeshReader *reader,
                         int s) {
    int *start, npe;

    if (reader->elem_type[s] == CGNS_ENUMV(MIXED)) {
        /* In the element type of a section is `MIXED`, the connectivity of each
           element is represented as:
               ElementType VertexIdx1 VertexIdx2 ... VertexIdxN.
           The index of ElementType in the connectivity array is stored in the
           offset array. */

        for (int i = 0; i < reader->nelems[s]; i++) {
            /* Pointer to ElementType. */
            start = &(reader->elem_conn[s][reader->elem_offset[s][i]]);

            /* Read the type as AerocasterMeshElementType. */
            mesh->elems[mesh->nelems].type = CGNSToAerocasterType(start[0]);
            npe = AerocasterMeshElementTypeNVerts[mesh->elems[mesh->nelems].type];

            /* Read the indices of vertices. Additional nodes in a quadratic,
               cubic, or quartic element are ignored. */
            for (int j = 0; j < npe; j++)
                mesh->elems[mesh->nelems].idx_verts[j] = start[j+1];

            /* Set its section. */
            mesh->elems[mesh->nelems].section = s;

            mesh->nelems++;
        }
    }
    else {
        /* For the other element types, the connectivity array contains only
           the indices of vertices. */

        /* Get the number of vertices in each element. */
        cg_npe(reader->elem_type[s], &npe);

        for (int i = 0; i < reader->nelems[s]; i++) {
            /* Pointer to the index of the first vertex. */
            start = &(reader->elem_conn[s][npe * i]);

            /* All elements have same type. */
            mesh->elems[mesh->nelems].type = CGNSToAerocasterType(reader->elem_type[s]);

            /* Read the indices of vertices. */
            for (int j = 0; j < npe; j++)
                mesh->elems[mesh->nelems].idx_verts[j] = start[j];

            /* Set its section. */
            mesh->elems[mesh->nelems].section = s;

            mesh->nelems++;
        }
    }
}

static GTree *GetAdjacency(AerocasterMesh *mesh) {
    GTree *face_tree;
    int *v;

    face_tree = g_tree_new_full(CompareFace, NULL, g_free, g_free);

    /* Build the face tree. */
    for (int i = 0; i < mesh->nelems; i++) {
        v = mesh->elems[i].idx_verts;

        switch (mesh->elems[i].type) {
        case AEROCASTER_TRI:
            StoreToFaceTree(face_tree, 2, v[0], v[1], i, 0);
            StoreToFaceTree(face_tree, 2, v[1], v[2], i, 1);
            StoreToFaceTree(face_tree, 2, v[2], v[0], i, 2);
            break;
        case AEROCASTER_QUAD:
            StoreToFaceTree(face_tree, 2, v[0], v[1], i, 0);
            StoreToFaceTree(face_tree, 2, v[1], v[2], i, 1);
            StoreToFaceTree(face_tree, 2, v[2], v[3], i, 2);
            StoreToFaceTree(face_tree, 2, v[3], v[0], i, 3);
            break;
        case AEROCASTER_TETRA:
            StoreToFaceTree(face_tree, 3, v[0], v[2], v[1], i, 0);
            StoreToFaceTree(face_tree, 3, v[0], v[1], v[3], i, 1);
            StoreToFaceTree(face_tree, 3, v[1], v[2], v[3], i, 2);
            StoreToFaceTree(face_tree, 3, v[2], v[0], v[3], i, 3);
            break;
        case AEROCASTER_PYRA:
            StoreToFaceTree(face_tree, 4, v[0], v[3], v[2], v[1], i, 0);
            StoreToFaceTree(face_tree, 3, v[0], v[1], v[4],       i, 1);
            StoreToFaceTree(face_tree, 3, v[1], v[2], v[4],       i, 2);
            StoreToFaceTree(face_tree, 3, v[2], v[3], v[4],       i, 3);
            StoreToFaceTree(face_tree, 3, v[3], v[0], v[4],       i, 4);
            break;
        case AEROCASTER_PRISM:
            StoreToFaceTree(face_tree, 4, v[0], v[1], v[4], v[3], i, 0);
            StoreToFaceTree(face_tree, 4, v[1], v[2], v[5], v[4], i, 1);
            StoreToFaceTree(face_tree, 4, v[2], v[0], v[3], v[5], i, 2);
            StoreToFaceTree(face_tree, 3, v[0], v[2], v[1],       i, 3);
            StoreToFaceTree(face_tree, 3, v[3], v[4], v[5],       i, 4);
            break;
        case AEROCASTER_HEXA:
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

static void ReadBoundary(AerocasterMesh *mesh,
                         AerocasterCGNSMeshReader *reader,
                         GTree *face_tree, int s) {
    int *start, tmp;
    AerocasterMeshElementType type;
    FaceKey key;
    FaceValue *value;

    if (reader->elem_type[s] == CGNS_ENUMV(MIXED)) {
        for (int i = 0; i < reader->nelems[s]; i++) {
            /* Pointer to ElementType. */
            start = &(reader->elem_conn[s][reader->elem_offset[s][i]]);

            /* Read the type as AerocasterMeshElementType. */
            type = CGNSToAerocasterType(start[0]);
            key.n = AerocasterMeshElementTypeNVerts[type];

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
            mesh->elems[value->elem_idx[0]].face_section[value->face_idx[0]]
                = s != -1 ? s : AEROCASTER_NO_ADJ;
            if (value->elem_idx[1] != -1)
                mesh->elems[value->elem_idx[1]].face_section[value->face_idx[1]] = s;
        }
    }

    else {
        /* Read the type as AerocasterMeshElementType. */
        type = CGNSToAerocasterType(reader->elem_type[s]);
        key.n = AerocasterMeshElementTypeNVerts[type];

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
            if (value->elem_idx[1] != -1)
                mesh->elems[value->elem_idx[1]].face_section[value->face_idx[1]] = s;
        }
    }
}

static AerocasterMeshElementType CGNSToAerocasterType(ElementType_t type) {
    switch (type) {
    case CGNS_ENUMV(BAR_2):
    case CGNS_ENUMV(BAR_3):
    case CGNS_ENUMV(BAR_4):
    case CGNS_ENUMV(BAR_5):
        return AEROCASTER_SEG;
    case CGNS_ENUMV(TRI_3):
    case CGNS_ENUMV(TRI_6):
    case CGNS_ENUMV(TRI_9):
    case CGNS_ENUMV(TRI_10):
    case CGNS_ENUMV(TRI_12):
    case CGNS_ENUMV(TRI_15):
        return AEROCASTER_TRI;
    case CGNS_ENUMV(QUAD_4):
    case CGNS_ENUMV(QUAD_8):
    case CGNS_ENUMV(QUAD_9):
    case CGNS_ENUMV(QUAD_12):
    case CGNS_ENUMV(QUAD_16):
    case CGNS_ENUMV(QUAD_P4_16):
    case CGNS_ENUMV(QUAD_25):
        return AEROCASTER_QUAD;
    case CGNS_ENUMV(TETRA_4):
    case CGNS_ENUMV(TETRA_10):
    case CGNS_ENUMV(TETRA_16):
    case CGNS_ENUMV(TETRA_20):
    case CGNS_ENUMV(TETRA_22):
    case CGNS_ENUMV(TETRA_34):
    case CGNS_ENUMV(TETRA_35):
        return AEROCASTER_TETRA;
    case CGNS_ENUMV(PYRA_5):
    case CGNS_ENUMV(PYRA_13):
    case CGNS_ENUMV(PYRA_14):
    case CGNS_ENUMV(PYRA_21):
    case CGNS_ENUMV(PYRA_29):
    case CGNS_ENUMV(PYRA_30):
    case CGNS_ENUMV(PYRA_P4_29):
    case CGNS_ENUMV(PYRA_50):
    case CGNS_ENUMV(PYRA_55):
        return AEROCASTER_PYRA;
    case CGNS_ENUMV(PENTA_6):
    case CGNS_ENUMV(PENTA_15):
    case CGNS_ENUMV(PENTA_18):
    case CGNS_ENUMV(PENTA_24):
    case CGNS_ENUMV(PENTA_38):
    case CGNS_ENUMV(PENTA_40):
    case CGNS_ENUMV(PENTA_33):
    case CGNS_ENUMV(PENTA_66):
    case CGNS_ENUMV(PENTA_75):
        return AEROCASTER_PRISM;
    case CGNS_ENUMV(HEXA_8):
    case CGNS_ENUMV(HEXA_20):
    case CGNS_ENUMV(HEXA_27):
    case CGNS_ENUMV(HEXA_32):
    case CGNS_ENUMV(HEXA_56):
    case CGNS_ENUMV(HEXA_64):
    case CGNS_ENUMV(HEXA_44):
    case CGNS_ENUMV(HEXA_98):
    case CGNS_ENUMV(HEXA_125):
        return AEROCASTER_HEXA;
    default:
        cg_error_exit();
        return -1;
    }
}

static void StoreToFaceTree(GTree *tree, int n, ...) {
    va_list ap;
    FaceKey key_static, *key;
    FaceValue *value;
    int tmp;

    va_start(ap, n);

    key_static.n = n;
    for (int i = 0; i < n; i++) key_static.idx[i] = va_arg(ap, int);

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
        if (value->elem_idx[1] != -1) cg_error_exit();

        value->elem_idx[1] = va_arg(ap, int);
        value->face_idx[1] = va_arg(ap, int);
    }
    else {
        key = malloc(sizeof(*key));
        value = malloc(sizeof(*value));

        memcpy(key, &key_static, sizeof(key_static));
        value->elem_idx[0] = va_arg(ap, int);
        value->face_idx[0] = va_arg(ap, int);
        value->elem_idx[1] = value->face_idx[1] = -1;

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
    for (int i = 0; i < fa->n; i++)
        if (fa->idx[i] != fb->idx[i]) return fa->idx[i] - fb->idx[i];
    return 0;
}

static gboolean FindAdjElem(gpointer key G_GNUC_UNUSED, gpointer value,
                            gpointer data) {
    FaceValue *fv;
    AerocasterMesh *mesh;

    fv = (FaceValue *)value;
    mesh = (AerocasterMesh *)data;

    mesh->elems[fv->elem_idx[0]].idx_adj[fv->face_idx[0]]
        = fv->elem_idx[1] != -1 ? fv->elem_idx[1] : AEROCASTER_NO_ADJ;
    if (fv->elem_idx[1] != -1)
        mesh->elems[fv->elem_idx[1]].idx_adj[fv->face_idx[1]] = fv->elem_idx[0];

    return false;
}
