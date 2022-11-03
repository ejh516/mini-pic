/*==============================================================================*
 * SPECIES
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-11-03
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#include "species.h"

#include "particle.h"

const double PLASMA_DEN = 1e10;
const double ION_VELOCITY = 7000;

Species::Species (int n_nodes) {
    density = new double[n_nodes];
    remainder = 0;
}

Species::~Species () {
    delete[] density;
}

void Species::inject(Mesh mesh, double dt) {
    /*set area of the k=0 face, this should be the sum of triangle areas on the inlet*/
    double area = 0.2*0.2;

    /*number of real ions per sec, given prescribed density and velocity*/
    double num_per_sec = PLASMA_DEN*ION_VELOCITY*area;

    /*number of ions to generate in this time step*/
    double num_real = num_per_sec*dt;

    /*fraction number of macroparticles*/
    double fnum_mp = num_real/spwt + remainder;

    /*integer number of macroparticles*/
    int num_mp = (int)fnum_mp;

    /*update reminder*/
    remainder = fnum_mp-num_mp;

    /*sample particles*/
    for (int p=0;p<num_mp;p++) {
        /*new particle*/
        Particle part;

        /*sample random position on the inlet face*/
        part.pos[0] = -0.1 + 0.2*rnd();
        part.pos[1] = -0.1 + 0.2*rnd();
        part.pos[2] = 0;

        /*injecting cold beam*/
        part.vel[0] = 0;
        part.vel[1] = 0;
        part.vel[2] = ION_VELOCITY;

        /*set initial tetrahedron*/
        for (part.cell_index = 0; part.cell_index<(int)mesh.elements.size(); part.cell_index++) {
            bool inside = part.getElement(mesh,false);
            if (!inside) continue;    /*go to next tetrahedron*/
            else break; /*break out once we found the tet*/
        }

        /*sanity check, should not happen*/
        if (part.cell_index>=(int)mesh.elements.size()) {
            throw "Failed to find initial element";
            continue;
        }

        /*rewind velocity*/
        double ef_part[3];
        solver.evalEf(ef_part, part.cell_index, part.weights);

        for (int i=0;i<3;i++)
            part.vel[i] -= charge/mass*ef_part[i]*(0.5*dt);

        /*add to list*/
        particles.push_back(part);
    }
}


void Species::updatePositions(Mesh mesh, double** e_field) {

}

void Species::updateDensity(Mesh mesh){

}
