/*==============================================================================*
 * SIMULATION
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-26
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>

#include "mesh.h"
#include "species.h"

class Simulation {
    public:
        // Attributes
        double time = 0; // seconds
        double max_time = 1e-9; // seconds
                                //
        Mesh mesh;
        Species ions;
        // Field<vector<double>> e_field;
        // Field<double> potential;

        // Methods
        Simulation(std::string filename);
        void run();
        void step();
        void writeOutputs();
        void ~Simulation(std::string filename);
};


#endif /* !SIMULATION_H */
