#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cgnslib.h>

int main(void) {
    int fn;
    int nbases, nzones, nsections, nelements;

    char base_name[50];
    int cell_dim, phys_dim;

    ZoneType_t zone_type;
    char zone_name[50];
    int zone_size[3];

    DataType_t coord_type;
    char coord_name[50];
    int rmin, rmax;
    double *x, *y, *z;

    ElementType_t elem_type;
    char sect_name[50];
    int start, end, nbndry, parent_flag;
    int elem_data_size;
    int **elems, **conn_offset;
    int parent_data;

    /* Read CGNS file. */
    if (cg_open("../example/StaticMixer.cgns", CG_MODE_READ, &fn)) {
        printf("error: file does not exist\n");
        cg_error_exit();
    }

    /* Read number of bases. It must be 1. */
    cg_nbases(fn, &nbases);
    if (nbases > 1) {
        printf("error: the number of bases must be 1\n");
        cg_error_exit();
    }

    /* Read base name and dimension. */
    cg_base_read(fn, 1, base_name, &cell_dim, &phys_dim);
    printf("base %s: cell_dim %d, phys_dim %d\n", base_name, cell_dim, phys_dim);

    /* Read number of zones. It must be 1. */
    cg_nzones(fn, 1, &nzones);
    if (nzones > 1) {
        printf("error: the number of zones must be 1\n");
        cg_error_exit();
    }

    /* Read zone type. It must be 'Unstructured'. */
    cg_zone_type(fn, 1, 1, &zone_type);
    if (zone_type != CGNS_ENUMV(Unstructured)) {
        printf("error: zone type must be unstructured\n");
        cg_error_exit();
    }

    /* Read zone name and sizes. Each size is the number of vertices, elements,
       and boundary vertices, respectively. */
    cg_zone_read(fn, 1, 1, zone_name, zone_size);
    printf("zone %s: #verts %d, #elems %d, #boundverts %d\n\n",
            zone_name, zone_size[0], zone_size[1], zone_size[2]);

    /* Read vertex coordinates. */
    rmin = 1;
    rmax = zone_size[0];
    x = malloc(sizeof(double) * rmax);
    y = malloc(sizeof(double) * rmax);
    z = malloc(sizeof(double) * rmax);
    cg_coord_read(fn, 1, 1, "CoordinateX", CGNS_ENUMV(RealDouble), &rmin, &rmax, x);
    cg_coord_read(fn, 1, 1, "CoordinateY", CGNS_ENUMV(RealDouble), &rmin, &rmax, y);
    cg_coord_read(fn, 1, 1, "CoordinateZ", CGNS_ENUMV(RealDouble), &rmin, &rmax, z);

    printf("first 10 points:\n");
    for (int i = 0; i < 10 && i < zone_size[0]; i++) {
        printf("  %8.4lf %8.4lf %8.4lf\n", x[i], y[i], z[i]);
    }
    printf("\n");

    /* Read number of sections. */
    cg_nsections(fn, 1, 1, &nsections);
    printf("#sections: %d\n", nsections);

    elems = malloc(sizeof(int *) * (nsections+1));
    conn_offset = malloc(sizeof(int *) * (nsections+1));
    for (int S = 1; S <= nsections; S++) {
        cg_section_read(fn, 1, 1, S, sect_name,
                        &elem_type, &start, &end, &nbndry, &parent_flag);
        printf("section %s: type %s, start %d, end %d, parent_flag %d\n",
                sect_name, ElementTypeName[elem_type], start, end, parent_flag);

        nelements = end - start + 1;

        cg_ElementDataSize(fn, 1, 1, S, &elem_data_size);
        printf("  element data size: %d\n", elem_data_size);

        elems[S] = malloc(sizeof(int) * elem_data_size);
        if (elem_type == CGNS_ENUMV(MIXED)) {
            conn_offset[S] = malloc(sizeof(int) * (nelements+1));
            cg_poly_elements_read(fn, 1, 1, S, elems[S], conn_offset[S], &parent_data);
        }
        else {
            cg_elements_read(fn, 1, 1, S, elems[S], NULL);
        }

        printf("  element data: ");
        for (int i = 0; i < 10; i++) {
            printf("%d ", elems[S][i]);
        }
        printf("\n");
        if (elem_type == CGNS_ENUMV(MIXED)) {
            printf("  offset: ");
            for (int i = 0; i < 10; i++) {
                printf("%d ", conn_offset[S][i]);
            }
            printf("\n");
        }
    }

    free(x); free(y); free(z);
    cg_close(fn);
    return 0;
}
