/*
 * data_class.cpp
 *
 *  Created on: Oct 7, 2013
 *      Author: john
 */
#include "Dataclass.hpp"

using namespace std;




NodeID_Dim::NodeID_Dim ()
{

}

NodeID_Dim::~NodeID_Dim ()
{

}

void NodeID_Dim::resize (unsigned nx, unsigned ny, unsigned nz, unsigned nxx, unsigned nyy, unsigned nzz, unsigned DumID)
{
    size_x = nx + nxx;
    size_y = ny + nyy;
    size_z = nz + nzz;

    arrNodeID.resize(size_x*size_y*size_z,DumID);
}

void NodeID_Dim::assign(unsigned x, unsigned y, unsigned z, unsigned NodeID)
{
    arrNodeID[x*(size_y*size_z)+ y*size_z + z] = NodeID;
}
unsigned NodeID_Dim::get(unsigned x, unsigned y, unsigned z)
{
    return arrNodeID[x*(size_y*size_z)+ y*size_z + z];
}

bool NodeLists::removeRestrain (const unsigned ID)
{
    list<Restrain>::iterator it;

    for (it = restrain.begin(); it != restrain.end(); ++it)
    {
        if ((*it).NodeID == ID)
        {
            restrain.erase(it);
            return true;
        }
    }
    return false;
}


//----------Lattice--------------//
