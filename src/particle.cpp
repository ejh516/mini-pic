/*==============================================================================*
 * PARTICLE
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-11-03
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#include "particle.h"

#include "element.h"

Particle::Particle(double pos[3]){

}

bool Particle::getElement(Mesh mesh, bool search) {

    /*first try the current element*/
    Element &elem = mesh.elements[cell_index];

    bool inside = true;

    try { weights = elem.getWeights(pos); }
    catch(bool in) { inside = in; } 

    if (inside) return true;

    if (!search) return false;

    /*we are outside the last known tet, find most negative weight*/
    int min_i=0;
    double min_weight = weights[0];

    for (int i=1;i<4;i++) {
        if (weights[i] < min_weight) {
            min_weight = weights[i];
            min_i = i;
        }
    }

    /*is there a neighbor in this direction?*/
    if (elem.cell_con[min_i]>=0) {
        cell_index = elem.cell_con[min_i];
        return getElement(mesh, true);
    }

    return false;
}
