/*==============================================================================*
 * MINI-PIC
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-26
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/
#include <string>
#include <iostream>

#include "simulation.h"

int main(int argc, char **argv) {
    // Get the input filename
    std::string filename;
    if (argc == 2) {
        filename = argv[1];
    } else {
        std::cerr << "Usage: fem-pic <FILENAME>" << std::endl;
        std::abort();
    }
    // Read in the simulation parameters from file
    //
    try {
        Simulation simulation;
        simulation.initialise(filename);

        simulation.run();

        simulation.write_outputs();

        return 0;
    } catch (char* msg) {
        std::cerr << "ERROR: " << msg << std::endl;
        std::abort();
    }

}
