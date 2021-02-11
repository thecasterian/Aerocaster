#include "../include/ReaderInterface.h"

int ReaderInterface_ReadFile(ReaderInterface *interface) {
    return interface->read_file(interface->reader);
}

void ReaderInterface_WriteToMesh(ReaderInterface *interface, Mesh *mesh) {
    interface->write_to_mesh(interface->reader, mesh);
}

void ReaderInterface_Destroy(ReaderInterface *interface) {
    interface->destroy(interface->reader);
}
