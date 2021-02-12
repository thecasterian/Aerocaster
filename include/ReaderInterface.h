#ifndef READER_INTERFACE_H
#define READER_INTERFACE_H

#define FILENAME_MAX_LEN 100
#define SECTNAME_MAX_LEN 50
#define READER_INTERFACE_ERROR -1

typedef struct _mesh Mesh;

typedef int (*ReaderInterface_ReadFileFunc)(void *reader);
typedef int (*ReaderInterface_WriteToMeshFunc)(void *reader, Mesh *mesh);
typedef void (*ReaderInterface_DestroyFunc)(void *reader);

typedef struct _reader_interface {
    void *reader;                                   /* Pointer to a reader implementing this interface. */
    char file_name[FILENAME_MAX_LEN];               /* File name to read. */

    ReaderInterface_ReadFileFunc read_file;         /* Function reading a file. */
    ReaderInterface_WriteToMeshFunc write_to_mesh;  /* Function writing reader's data to the mesh. */
    ReaderInterface_DestroyFunc destroy;            /* Function destroying the reader. */
} ReaderInterface;

int ReaderInterface_ReadFile(ReaderInterface *interface);
int ReaderInterface_WriteToMesh(ReaderInterface *interface, Mesh *mesh);
void ReaderInterface_Destroy(ReaderInterface *interface);

#endif
