/*==============================================================================*
 * ELEMENT
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-11-03
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#include "element.h"

Element::Element(int n1, int n2, int n3, int n4) {
    con[0]=n1;
    con[1]=n2;
    con[2]=n3;
    con[3]=n4;
}

std::vector<double> Element::getWeights(double pos[3]) {
    std::vector<double> weights(4);

    /*loop over vertices*/
    for (int i=0;i<4;i++) {
        weights[i] = (1.0/6.0)*(alpha[i]
                              - pos[0]*beta[i]
                              + pos[1]*gamma[i]
                              - pos[2]*delta[i])/volume;

        if (weights[i]<0 || weights[i]>1.0) throw false;
    }

    return weights;
}
