#ifndef DEFORMATION_TRANSFER_H
#define DEFORMATION_TRANSFER_H
#include <vector>
#include <utility>

class TriMesh;

//Before call this function, make sure the initial aligment is done
void deform_transfer(TriMesh& src, const TriMesh& dst, 
    const std::vector<std::pair<int, int> >& corres);

#endif/*DEFORMATION_TRANSFER_H*/
