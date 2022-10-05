/*==============================================================================*
 * PARTICLES
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>

/*particle*/
struct Particle {
    double pos[3];
    double vel[3];
    double lc[4];    /*particle's weights*/
    int cell_index;    /*last cell known to contain this particle*/
};

/*species class*/
class Species {
public:
    std::vector<Particle> particles;
    double *den;
    double spwt;
    double mass;
    double charge;
    double rem;        /*fractional particle reminder*/

    Species (int n_nodes) {den=new double[n_nodes];rem=0;}
    ~Species () {delete[] den;}
};

void OutputParticles(std::vector<Particle> &particles);

#endif /* !PARTICLES_H */
