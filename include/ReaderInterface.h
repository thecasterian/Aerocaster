#ifndef READER_INTERFACE_H
#define READER_INTERFACE_H

#define NAME_MAX_LEN 50
#define READER_ERROR -1

typedef struct _mesh Mesh;

typedef struct _reader_interface {
    void *reader;                               /* Pointer to a reader implementing this interface. */
    int (*read_file)(void *, const char *);     /* Function reading a file. */
    void (*write_to_mesh)(void *, Mesh *);      /* Function writing reader's data to the mesh. */
    void (*destroy)(void *);                    /* Function destroying the reader. */
} ReaderInterface;

void ReaderInterface_ReadFile(ReaderInterface *interface,
                              const char *file_name);
void ReaderInterface_WriteToMesh(ReaderInterface *interface, Mesh *mesh);
void ReaderInterface_Destroy(ReaderInterface *interface);

#endif
