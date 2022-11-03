/*==============================================================================*
 * ELEMENT
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-11-03
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>

/*definition of a tetrahedral element*/
class Element {
    public:
        // Attributes
        int con[4];
        double volume;

        /*data structures to hold precomputed 3x3 determinants*/
        double alpha[4], beta[4], gamma[4], delta[4];

        /*cell connectivity*/
        int cell_con[4];    /*index corresponds to the face opposite the i-th node*/

        // Methods
        Element (int n1, int n2, int n3, int n4);
        std::vector<double> getWeights(double pos[3]);
};


#endif /* !ELEMENT_H */
