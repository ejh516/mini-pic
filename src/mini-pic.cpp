/* Example of a Finite Element ES-PIC Code

  Written by Lubos Brieda for Advanced PIC 2015 Lesson 8
  See https://www.particleincell.com/2015/fem-pic/ for more information

  To compile and run:
    g++ -std=c++10 -O2 fem-pic.cpp -o fem-pic
    ./fem-pic

*/

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include "trace.h"
#include "maths.h"
#include "particles.h"
#include "meshes.h"
#include "FESolver.h"

/*constants*/
const double PLASMA_DEN = 1e10;
const double ION_VELOCITY = 7000;

int trace::enabled = 1;
Trace trace::current = Trace("__TRACE_BASE__");


/*PROTOTYPES*/
void InjectIons(Species &ions, Volume &volume, FESolver &solver, double dt);
void MoveParticles(Species &ions, Volume &volume, FESolver &solver, double dt);


/**************** MAIN **************************/
int main() {

    /*instantiate volume*/
    Volume volume;
    if (!LoadVolumeMesh("mesh.dat",volume) ||
        !LoadSurfaceMesh("inlet.dat",volume,INLET) ||
        !LoadSurfaceMesh("fixed.dat",volume,FIXED)) return -1;

    /*instantiate solver*/
    FESolver solver(volume);

    /*set reference paramaters*/
    solver.phi0 = 0;
    solver.n0 = PLASMA_DEN;
    solver.kTe = 2;

    int n_nodes = volume.nodes.size();


    /*initialize solver "g" array*/
    for (int n=0;n<n_nodes;n++) {
        if (volume.nodes[n].type==INLET) solver.g[n]=0;    /*phi_inlet*/
        else if (volume.nodes[n].type==FIXED) solver.g[n]=-100; /*fixed phi points*/
        else solver.g[n]=0;    /*default*/
    }

    /*sample assembly code*/
    solver.startAssembly();
    solver.preAssembly();    /*this will form K and F0*/

    double dt = 1e-7;

    /*ions species*/
    Species ions(n_nodes);
    ions.charge=1*QE;
    ions.mass = 16*AMU;
    ions.spwt = 2e2;

    /*main loop*/
    int ts;
    for (ts=0;ts<100;ts++) {
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

        if ((ts+1)%10==0) OutputMesh(ts,volume, solver.uh, solver.ef, ions.den);

        std::cout<<"ts: "<<ts<<"\t np:"<<ions.particles.size()<<"\t max den:"<<max_den<<std::endl;
    }

    /*output mesh*/
    OutputMesh(ts,volume, solver.uh, solver.ef, ions.den);

    /*output particles*/
    OutputParticles(ions.particles);
    if (trace::enabled) trace::current.write_profile("trace.csv");

    return 0;
}


/*** FUNCTIONS ***/



/*** Samples particle on the inlet surface
here we just sample on a known plane. A generic code should instead sample
from the surface triangles making up the inlet face*/
 void InjectIons(Species &ions, Volume &volume, FESolver &solver, double dt) { TRACE_ME;
    /*set area of the k=0 face, this should be the sum of triangle areas on the inlet*/
    for (auto face: volume.inlet_faces) {

        /*number of real ions per sec, given prescribed density and velocity*/
        double num_per_sec = PLASMA_DEN*ION_VELOCITY*face.area;

        /*number of ions to generate in this time step*/
        double num_real = num_per_sec*dt;

        /*fraction number of macroparticles*/
        double fnum_mp = num_real/ions.spwt + ions.rem;

        /*integer number of macroparticles*/
        int num_mp = (int)fnum_mp;

        /*update reminder*/
        ions.rem = fnum_mp-num_mp;

        /*sample particles*/
        for (int p=0;p<num_mp;p++) {
            /*new particle*/
            Particle part;

            /*sample random position on the inlet face*/
            double a = rnd();
            double b = rnd();
            if ((a+b) > 1)  {
                a = 1-a;
                b = 1-b;
            }

            for (int i=0; i<3; i++) {
                part.pos[i] = a*face.u[i] + b*face.v[i] + volume.nodes[face.con[0]].pos[i];

                /*injecting cold beam*/
                part.vel[i] = face.normal[i] * ION_VELOCITY;
            }

            /*set initial tetrahedron*/
            part.cell_index = face.cell_con;

            /*rewind velocity*/
            double ef_part[3];
            solver.evalEf(ef_part, part.cell_index);

            for (int i=0;i<3;i++)
                part.vel[i] -= ions.charge/ions.mass*ef_part[i]*(0.5*dt);

            /*add to list*/
            ions.particles.push_back(part);
        }
    }
}

/*updates ion velocities and positions*/
void MoveParticles(Species &ions, Volume &volume, FESolver &solver, double dt) {    TRACE_ME;
    int n_nodes = (int) volume.nodes.size();

    /*reset ion density*/
    for (int i=0;i<n_nodes;i++) ions.den[i] = 0;

    /*move particles*/
    auto part_it = ions.particles.begin();
    while(part_it != ions.particles.end()) {
        Particle &part = *part_it;

        /*update particle velocity*/
        double ef_part[3];
        solver.evalEf(ef_part, part.cell_index);

        for (int i=0;i<3;i++)
            part.vel[i] += ions.charge/ions.mass*ef_part[i]*dt;

        /*update particle positions*/
        for (int i=0;i<3;i++) part.pos[i]+=part.vel[i]*dt;

        bool inside = XtoLtet(part,volume);

        if (inside) {
            Tetra &tet = volume.elements[part.cell_index];
            /*now we know that we are inside this tetrahedron, scatter*/
            double sum=0;
            for (int v=0;v<4;v++) {
                ions.den[tet.con[v]]+=part.lc[v];
                sum+=part.lc[v];    /*for testing*/
            }

            /*testing*/
            if (std::abs(sum-1.0)>0.001) std::cout<<sum<<std::endl;

            part_it++;
        }
        else part_it = ions.particles.erase(part_it);    /*outside the mesh*/
    }

    /*convert to ion density*/
    for (int n=0;n<n_nodes;n++) ions.den[n] *= ions.spwt/volume.nodes[n].volume;
}

