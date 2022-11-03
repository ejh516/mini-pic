/*==============================================================================*
 * MESH
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-11-03
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef MESH_H
#define MESH_H

#include <vector>

#include "node.h"
#include "element.h"

class Mesh {
    public:
        std::vector <Node> nodes;
        std::vector <Element> elements;
};


#endif /* !MESH_H */
