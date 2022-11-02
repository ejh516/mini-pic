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

#include "volume.h"

class Simulation {
    public:
        // Attributes
        double time = 0; // seconds
        double max_time = 1e-9; // seconds
        Volume volume;

        // Methods
        int initialise(std::string filename);
        void run();
        void step();
        int write_outputs();
};


#endif /* !SIMULATION_H */
