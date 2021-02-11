#ifndef AEROCASTER_MESH_H
#define AEROCASTER_MESH_H

#include "ReaderInterface.h"
#include <mpi.h>

typedef struct _mesh_vert {
    double x, y, z;                     /* Coordinates. */
} MeshVert;

typedef enum _mesh_elem_type {
    ELEMTYPE_SEG,                       /* Line segment. */
    ELEMTYPE_TRI,                       /* Triangle. */
    ELEMTYPE_QUAD,                      /* Quadrilateral. */
    ELEMTYPE_TETRA,                     /* Tetrahedron. */
    ELEMTYPE_PYRA,                      /* Pyramid. */
    ELEMTYPE_PRISM,                     /* Triangular prism. */
    ELEMTYPE_HEXA,                      /* Hexahedron. */
} MeshElemType;

extern const int MeshElemTypeNVerts[7];
extern const int MeshElemTypeNFaces[7];

/* The value of idx_adj of an element if the adjacent element does not exist. */
#define IDX_ADJ_NO_ADJ ((size_t)(-1))
/* The value of face_section of an element if the section is unspecified and the
   adjacent element does not exist. */
#define FACE_SECT_UNSPEC_BNDRY -1
/* The value of face_section of an element if the section is unspecified and the
   adjacent element exists. */
#define FACE_SECT_UNSPEC_INTER -2

typedef struct _mesh_elem {
    MeshElemType type;                  /* Element type. */
    int section;                        /* Section number this element belongs to. */

    size_t idx_verts[8];                /* Indices of vertices. */
    size_t idx_adj[6];                  /* Indices of adjacent elements. */
    int face_section[6];                /* Section number each face belongs to. */

    double cx, cy, cz;                  /* Coordinates of centroid. */
    double v;                           /* Volume (area for 2-D). */
    double face_area[6];                /* Area (length for 2-D) of faces. */
} MeshElem;

typedef struct _mesh {
    ReaderInterface *reader;            /* Mesh file reader. */
    MPI_Comm comm;                      /* MPI communicator. */

    int dim;                            /* Mesh dimension. */
    size_t nverts;                      /* Number of vertices. */
    size_t nelems;                      /* Number of elements. */
    MeshVert *verts;                    /* Vertices. */
    MeshElem *elems;                    /* Elements. */

    int nsects;                         /* Number of sections. */
    char (*sect_name)[NAME_MAX_LEN];    /* Name of sections. */
} Mesh;

Mesh *Mesh_Create(ReaderInterface *reader, MPI_Comm comm);
void Mesh_Destroy(Mesh *mesh);

#endif
