#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cgnslib.h>

int main(void) {
    int fn;
    int nbases, nzones, nsections;

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

    /* Read CGNS file. */
    if (cg_open("../example/StaticMixer.cgns", CG_MODE_READ, &fn)) {
        printf("error: file does not exist\n");
        return -1;
    }

    /* Read number of bases. */
    cg_nbases(fn, &nbases);
    printf("#bases: %d\n", nbases);

    /* Read number of zones. */
    cg_nzones(fn, 1, &nzones);
    printf("#zones: %d\n", nzones);

    for (int Z = 1; Z <= nzones; Z++) {
        /* Read zone type, which must be 'Unstructured'. */
        cg_zone_type(fn, 1, Z, &zone_type);
        if (zone_type != CGNS_ENUMV(Unstructured)) {
            printf("error: zone type must be unstructured\n");
        }

        /* Read zone name and size. */
        cg_zone_read(fn, 1, Z, zone_name, zone_size);
        printf("zone %s: #verts %d, #cells %d, #boundverts %d\n",
               zone_name, zone_size[0], zone_size[1], zone_size[2]);

        /* Read vertex coordinates. */
        rmin = 1;
        rmax = zone_size[0];
        x = malloc(sizeof(double) * rmax);
        y = malloc(sizeof(double) * rmax);
        z = malloc(sizeof(double) * rmax);
        cg_coord_read(fn, 1, Z, "CoordinateX", CGNS_ENUMV(RealDouble), &rmin, &rmax, x);
        cg_coord_read(fn, 1, Z, "CoordinateY", CGNS_ENUMV(RealDouble), &rmin, &rmax, y);
        cg_coord_read(fn, 1, Z, "CoordinateZ", CGNS_ENUMV(RealDouble), &rmin, &rmax, z);

        /* Read number of sections. */
        cg_nsections(fn, 1, Z, &nsections);
        printf("#sections: %d\n", nsections);

        for (int S = 1; S <= nsections; S++) {
            cg_section_read(fn, 1, Z, S, sect_name,
                            &elem_type, &start, &end, &nbndry, &parent_flag);
            printf("section %s: type %s, start %d, end %d\n",
                   sect_name, ElementTypeName[elem_type], start, end);


        }
    }

    cg_close(fn);
    return 0;
}
