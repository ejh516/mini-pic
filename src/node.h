/*==============================================================================*
 * NODE
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-11-03
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef NODE_H
#define NODE_H


/*definition of a node*/
class Node {
    /*node type*/
    enum Type {NORMAL,OPEN,INLET,SPHERE};

    double pos[3];    /*node position*/
    double volume;    /*node volume*/
    Type type;

    Node(double x, double y, double z);
};


#endif /* !NODE_H */
