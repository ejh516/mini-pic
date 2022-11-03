/*==============================================================================*
 * PARTICLE_H
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-26
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include "mesh.h"

class Particle {
    public:
        // Attributes
        double pos[3];
        double vel[3];
        std::vector<double> weights;    /*particle's weights*/
        int cell_index;    /*last cell known to contain this particle*/

        // Methods
        Particle(double pos[3]);
        bool getElement(Mesh mesh, bool search);
};

#endif /* !PARTICLE_H */
