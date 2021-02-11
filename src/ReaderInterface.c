#include "../include/ReaderInterface.h"

void ReaderInterface_ReadFile(ReaderInterface *interface,
                              const char *file_name) {
    interface->read_file(interface->reader, file_name);
}

void ReaderInterface_WriteToMesh(ReaderInterface *interface, Mesh *mesh) {
    interface->write_to_mesh(interface->reader, mesh);
}

void ReaderInterface_Destroy(ReaderInterface *interface) {
    interface->destroy(interface->reader);
}
