#ifndef AEROCASTER_MESH_H
#define AEROCASTER_MESH_H

#include "AerocasterCGNSReader.h"

typedef struct {
    double x, y, z;                     /* Coordinates. */
} AerocasterMeshVertex;

typedef enum {
    AEROCASTER_SEG,                     /* Line segment. */
    AEROCASTER_TRI,                     /* Triangle. */
    AEROCASTER_QUAD,                    /* Quadrilateral. */
    AEROCASTER_TETRA,                   /* Tetrahedron. */
    AEROCASTER_PYRA,                    /* Pyramid. */
    AEROCASTER_PRISM,                   /* Triangular prism. */
    AEROCASTER_HEXA,                    /* Hexahedron. */
} AerocasterMeshElementType;

extern const int AerocasterMeshElementTypeNVerts[7];
extern const int AerocasterMeshElementTypeNFaces[7];

/* The value of idx_adj of an element if the adjacent element does not exist. */
#define AEROCASTER_NO_ADJ -1
/* The value of face_section of an element if the section is unspecified and the
adjacent element does not exist. */
#define AEROCASTER_UNSPEC_BNDRY -1
/* The value of face_section of an element if the section is unspecified and the
adjacent element exists. */
#define AEROCASTER_UNSPEC_INTER -2

typedef struct {
    AerocasterMeshElementType type;     /* Element type. */
    int section;                        /* Section number this element belongs to. */

    int idx_verts[8];                   /* Indices of vertices. */
    int idx_adj[6];                     /* Indices of adjacent elements. */
    int face_section[6];                /* Section number each face belongs to. */

    double cx, cy, cz;                  /* Coordinates of centroid. */
    double v;                           /* Volume (area for 2-D). */
    double face_area[6];                /* Area (length for 2-D) of faces. */
} AerocasterMeshElement;

typedef struct {
    int dim;                            /* Mesh dimension. */
    int nverts;                         /* Number of vertices. */
    int nelems;                         /* Number of elements. */
    AerocasterMeshVertex *verts;        /* Vertices. */
    AerocasterMeshElement *elems;       /* Elements. */

    int nsects;                         /* Number of sections. */
    char (*sect_name)[NAME_MAX_LEN];    /* Name of sections. */
} AerocasterMesh;

AerocasterMesh *AerocasterMesh_Create(void);
void AerocasterMesh_ReadCGNSMeshReader(AerocasterMesh *mesh,
                                       AerocasterCGNSMeshReader *reader);
void AerocasterMesh_Destroy(AerocasterMesh *mesh);

#endif