/*==============================================================================*
 * SIMULATION
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-26
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#include <iostream>
#include "simulation.h"


int Simulation::initialise(std::string filename) {
r   return 0;
}

void Simulation::step() {
    /*sample new particles*/
    InjectIons(ions, volume, solver, dt);

    /*update velocity and move particles*/
    MoveParticles(ions, volume, solver, dt);

    /*check values*/
    double max_den=0;
    for (int n=0;n<n_nodes;n++) if (ions.den[n]>max_den) max_den=ions.den[n];

    /*call potential solver*/
    solver.computePhi(ions.den);

    solver.updateEf();

    simulation.time += simulation.dt;
}

void Simulation::run() {
    int timestep = 0;
    while (time < max_time) {
        simulation.step();
        if ((timestep+1)%10==0) OutputMesh(timestep, volume, solver.uh, solver.ef, ions.den);

        std::cout<<"Timestep: " << timestep
                 <<"\t np:"     << ions.particles.size()
                 <<"\t max den:"<< max_den<<std::endl;

        timestep++;
    }


}

