/*==============================================================================*
 * MESHES
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef MESHES_H
#define MESHES_H

#include <vector>
#include <string>

#include "particles.h"

/*node type*/
enum NodeType {NORMAL,OPEN,INLET,SPHERE};

/*definition of a node*/
struct Node {
    Node(double x, double y, double z) {pos[0]=x;pos[1]=y;pos[2]=z;type=NORMAL;}
    double pos[3];    /*node position*/
    NodeType type;
    double volume;    /*node volume*/
};

/*definition of a tetrahedron*/
struct Tetra {
    int con[4];
    double volume;
    Tetra (int n1, int n2, int n3, int n4) {con[0]=n1;con[1]=n2;con[2]=n3;con[3]=n4;}

    /*data structures to hold precomputed 3x3 determinants*/
    double alpha[4], beta[4], gamma[4], delta[4];

    /*cell connectivity*/
    int cell_con[4];    /*index corresponds to the face opposite the i-th node*/
};

/*definition of a volume*/
struct Volume {
    std::vector <Node> nodes;
    std::vector <Tetra> elements;
};

bool LoadVolumeMesh(const std::string file_name, Volume &volume);
bool LoadSurfaceMesh(const std::string file_name, Volume &volume, NodeType node_type);
void OutputMesh(int ts, Volume &volume, double *phi, double **ef, double *ion_den);

bool XtoLtet(Particle &part, Volume &volume, bool search=true);

#endif /* !MESHES_H */
