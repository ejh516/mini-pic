/*==============================================================================*
 * SPECIES
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-11-03
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef SPECIES_H
#define SPECIES_H


#include "particle.h"

#include  <vector>

#include "mesh.h"

/*species class*/
class Species {
    public:
        std::vector<Particle> particles;
        double* density;
        double spwt;
        double mass;
        double charge;
        double remainder;        /*fractional particle reminder*/

        Species (int n_nodes);
        ~Species ();

        void inject(Mesh mesh, double dt);

        void updatePositions(Mesh mesh, double** e_field);

        void updateDensity(Mesh mesh);

};

#endif /* !SPECIES_H */
