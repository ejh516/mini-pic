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


Simulation::Simulation(std::string filename) {
    throw "Simulation::initialise not yet implemented";
}

Simulation::~Simulation() {

}

void Simulation::step() {
    ions.inject(mesh);

    ions.updatePositions(mesh, e_field);

    ions.updateDensity(mesh)
    
    potential.calculate(mesh, ions.density);

    e_field.calculate(mesh, potential);

}

void Simulation::writeOutputs() {
    throw "Simulation::writeOutputs not yet implemented";
}

void Simulation::run() {
    int timestep = 0;
    while (time < max_time) {
        step();
        if ((timestep+1)%10==0) writeOutputs();

        std::cout<<"Timestep: " << timestep
                 <<"\t np:"     << ions.particles.size()
                 <<"\t max den:"<< ions.max_den << std::endl;

        timestep++;
    }
}

