/*==============================================================================*
 * FESOLVER
 *------------------------------------------------------------------------------*
 * Author:  Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/


#include <iostream>
#include <iomanip>
#include <math.h>

#include "trace.h"
#include "maths.h"
#include "FESolver.h"

/*FESolver*/
FESolver::FESolver(Volume &volume):volume(volume) { TRACE_ME;
    /*count number of unknowns*/
    neq = 0;

    /*OPEN nodes are "h" nodes*/
    for (std::size_t i=0;i<volume.nodes.size();i++)
        if (volume.nodes[i].type==NORMAL ||
            volume.nodes[i].type==OPEN) neq++;

    std::cout<<"There are "<<neq<<" unknowns"<<std::endl;
    /*allocate neq*neq K matrix*/
    K = new double*[neq];
    for (int i=0;i<neq;i++) K[i] = new double[neq];
    std::cout<<"Allocated "<<neq<<"x"<<neq<<" stiffness matrix"<<std::endl;

    /*allocate neq*neq J matrix*/
    J = new double*[neq];
    for (int i=0;i<neq;i++) J[i] = new double[neq];
    std::cout<<"Allocated "<<neq<<"x"<<neq<<" Jacobian matrix"<<std::endl;

    /*allocate F0 and F1 vectors*/
    F0 = new double[neq];
    F1 = new double[neq];
    std::cout<<"Allocated two "<<neq<<"x1 force vectors"<<std::endl;

    n_nodes = volume.nodes.size();
    n_elements = volume.elements.size();

    /*allocate ID vector*/
    ID = new int[n_nodes];
    std::cout<<"Allocated "<<n_nodes<<"x1 ID vector"<<std::endl;

    /*allocate location matrix, n_elements*4 */
    LM = new int*[n_elements];
    for (int e=0;e<n_elements;e++) LM[e] = new int[4];
    std::cout<<"Allocated "<<n_elements<<"x4 location matrix"<<std::endl;

    /*allocate NX matrix*/
    NX = new double**[n_elements];
    for (int e=0;e<n_elements;e++) {
        NX[e] = new double*[4];
        for (int a=0;a<4;a++) NX[e][a] = new double[3];
    }
    std::cout<<"Allocated "<<n_elements<<"x4x3 NX matrix"<<std::endl;

    /*solution array*/
    d = new double[neq];
    for (int n=0;n<neq;n++) d[n]=0;    /*initial guess*/

    /*allocate memory for g and uh arrays*/
    g = new double[n_nodes];
    uh = new double[n_nodes];
    std::cout<<"Allocated "<<n_nodes<<"x1 g and uh vector"<<std::endl;

    detJ = new double[n_elements];
    std::cout<<"Allocated "<<n_elements<<"x1 detJ vector"<<std::endl;

    /*electric field*/
    ef = new double*[n_elements];
    for (int n=0;n<n_elements;n++)
        ef[n] = new double[3];

    /*set up the ID array
    note valid values are 0 to neq-1 and -1 indicates "g" node*/
    int P=0;
    for (int n=0;n<n_nodes;n++)
        if (volume.nodes[n].type==NORMAL ||
            volume.nodes[n].type==OPEN) {ID[n]=P;P++;}
        else ID[n]=-1;    /*dirichlet node*/

    /*now set up the LM matrix*/
    for (int e=0;e<n_elements;e++)
        for (int a=0;a<4;a++)    /*tetrahedra*/ {
            LM[e][a] = ID[volume.elements[e].con[a]];
        }
    std::cout<<"Built ID and LM matrix"<<std::endl;

    /*set quadrature points*/
    l[0]=-sqrt(1.0/3.0); l[1]=sqrt(1.0/3.0);
    W[0]=1; W[1]=1;
    n_int = 2;

    /*compute NX matrix*/
    computeNX();
}

/*~FESolver, frees memory*/
FESolver::~FESolver() { TRACE_ME;
    for (int i=0;i<neq;i++) {delete[] K[i]; delete[] J[i];}
    for (int e=0;e<n_elements;e++) delete[] LM[e];
    for (int e=0;e<n_elements;e++) {
        for (int a=0;a<4;a++) delete[] NX[e][a];
        delete NX[e];
    }

    for (int e=0;e<n_elements;e++)    delete[] ef[e];
    delete[] ef;

    delete[] K;
    delete[] J;
    delete[] LM;
    delete[] F0;
    delete[] F1;
    delete[] NX;
    delete[] ID;
    delete[] uh;
    delete[] d;
    delete[] g;
    delete[] detJ;
}

/*clears K and F*/
void FESolver::startAssembly() { TRACE_ME;
    for (int i=0;i<neq;i++)
        for (int j=0;j<neq;j++) K[i][j] = 0;

    for (int i=0;i<neq;i++) {
        F0[i]=0;
        F1[i]=0;
    }
}

/*adds contributions from element stiffness matrix*/
void FESolver::addKe(int e, double ke[4][4]) { TRACE_ME;
    for (int a=0;a<4;a++)    /*tetrahedra*/
        for (int b=0;b<4;b++) {
            int P = LM[e][a];
            int Q = LM[e][b];
            if (P<0 || Q<0) continue;    /*skip g nodes*/

            K[P][Q] += ke[a][b];
        }
}

/*adds contributions from element force vector to a global F vector*/
void FESolver::addFe(double *F, int e, double fe[4]) { TRACE_ME;
    for (int a=0;a<4;a++)    /*tetrahedra*/ {
        int P = LM[e][a];
        if (P<0) continue;    /*skip g nodes*/

        F[P] += fe[a];
    }
}

/*evaluates shape function a at position (xi,eta,zeta)*/
double FESolver::evalNa(int a, double xi, double eta, double zeta) {
    switch(a) {
    case 0: return xi; break;
    case 1: return eta; break;
    case 2: return zeta; break;
    case 3: return 1-xi-eta-zeta; break;
    default: return 0;    /*shouldn't happen*/
    }
}

/*returns derivative of N[a] at some logical point
since we are using linear elements, these are constant in each element*/
void FESolver::getNax(double nx[3], int e, int a) {
    for (int d=0;d<3;d++)
        nx[d] = NX[e][a][d];
}

/*computes derivatives of the shape functions for all elements
constants since using linear elements*/
void FESolver::computeNX() { TRACE_ME;
    /*derivatives of the shape functions vs. xi*/
    double na_xi[4][3] = {{1,0,0}, {0,1,0}, {0,0,1}, {-1,-1,-1}};

    for (int e=0;e<n_elements;e++) {
        /*node positions*/
        Tetra &tet = volume.elements[e];

        double x[4][3];
        for (int a=0;a<4;a++) {
            double *pos = volume.nodes[tet.con[a]].pos;
            for (int d=0;d<3;d++) x[a][d] = pos[d];    /*copy*/
        }

        /*compute x_xi matrix*/
        double x_xi[3][3];

        for (int i=0;i<3;i++)    /*x/y/z*/
            for (int j=0;j<3;j++) /*xi/eta/zeta*/ {
                x_xi[i][j] = 0;
                for (int a=0; a<4; a++)    /*tet node*/
                    x_xi[i][j] += na_xi[a][j]*x[a][i];
            }

        /*save det(x_xi)*/
        detJ[e] = det3(x_xi);

        /*compute matrix inverse*/
        double xi_x[3][3];
        inverse(x_xi,xi_x);

        /*evaluate na_x*/
        for (int a=0;a<4;a++)
            for (int d=0;d<3;d++) {
                NX[e][a][d]=0;
                for (int k=0;k<3;k++)
                    NX[e][a][d]+=na_xi[a][k]*xi_x[k][d];
            }
    }
}

/*compute inverse of a 3x3 matrix using the adjugate method*/
void FESolver::inverse(double M[3][3], double V[3][3]) { TRACE_ME;
    double a=M[0][0];
    double b=M[0][1];
    double c=M[0][2];
    double d=M[1][0];
    double e=M[1][1];
    double f=M[1][2];
    double g=M[2][0];
    double h=M[2][1];
    double i=M[2][2];

    V[0][0]=(e*i-f*h);
    V[1][0]=-(d*i-f*g);
    V[2][0]=(d*h-e*g);
    V[0][1]=-(b*i-c*h);
    V[1][1]=(a*i-c*g);
    V[2][1]=-(a*h-b*g);
    V[0][2]=(b*f-c*e);
    V[1][2]=-(a*f-c*d);
    V[2][2]=(a*e-b*d);
    double det = a*V[0][0]+b*V[1][0]+c*V[2][0];

    double idet=0;
    if (fabs(det)<1e-12) {
        std::cerr<<"Matrix is not invertible, setting to [0]!"<<std::endl;}
    else idet=1/det;

    /*1/det*/
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++)
            V[i][j]*=idet;
}

/*Newton Rhapson solver, input is the ion density*/
void FESolver::solveNonLinear(double *ion_den) { TRACE_ME;

    /*allocate memory for y*/
    double *y = new double[neq];
    double *G = new double[neq];

    /*clear y values*/
    for (int i=0;i<neq;i++) y[i]=0;

    bool converged = false;
    double L2;
    for (int it=0;it<10;it++) {
        /*builds the "ff" part of the force vector*/
        buildF1Vector(ion_den);

        /*form G=K*d-F*/
        matVecMultiply(G,K,d, neq);    //G=K*d
        vecVecSubtract(G,G,F0, neq); //G=G-F giving us G=K*d-F
        vecVecSubtract(G,G,F1, neq);

        buildJmatrix();

        solveLinear(J,y,G);

        /*now that we have y, update solution */
        for (int n=0;n<neq;n++) d[n]-=y[n];

        /*compute residue*/
        double sum=0;
        for (int u=0;u<neq;u++) {
            sum+=y[u]*y[u];
        }
        L2 = sqrt(sum)/neq;
        if (L2<1e-2) {
            std::cout<<" NR converged in "<<it+1<<" iterations with L2="<<std::setprecision(3)<<L2<<std::endl;
            converged=true;
            break;
        }

    }

    if (!converged) std::cerr<<"NR failed to converge, L2 = "<<L2<<std::endl;
    delete[] y;
    delete[] G;

}

/*builds J matrix for NR solver*/
void FESolver::buildJmatrix() { TRACE_ME;
    /*first compute exponential term*/
    double *fp_term = new double[neq];
    double *FP = new double[neq];

    for (int n=0;n<neq;n++) FP[n] = 0;

    for (int n=0;n<neq;n++) {
        fp_term[n] = -QE/EPS0*n0*exp((d[n]-phi0)/kTe)*(1/kTe);
    }

    /*now set J=K*/
    for (int i=0;i<neq;i++)
        for (int j=0;j<neq;j++)
            J[i][j] = K[i][j];

    /*build fprime vector*/
    double fe[4];

    for (int e=0;e<n_elements;e++) {
        for (int a=0;a<4;a++) {
            double ff=0;
            int A = LM[e][a];
            if (A>=0)    /*if unknown node*/ {
                /*perform quadrature*/
                for (int k=0;k<n_int;k++)
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) {
                            /*change of limits*/
                            double xi = 0.5*(l[i]+1);
                            double eta = 0.5*(l[j]+1);
                            double zeta = 0.5*(l[k]+1);

                            double Na=evalNa(a,xi,eta,zeta);
                            ff += fp_term[A]*Na*detJ[e]*W[i]*W[j]*W[k];
                        }
                ff*=(1.0/8.0);    /*change of limits*/
            }
            fe[a] = ff;
        }

        /*assembly*/
        for (int a=0;a<4;a++)    /*tetrahedra*/ {
            int P = LM[e][a];
            if (P<0) continue;    /*skip g nodes*/

            FP[P] += fe[a];
        }
    }


    /*subtract diagonal term*/
    for (int u=0;u<neq;u++) {
        J[u][u]-=FP[u];
    }

    delete[] fp_term;
    delete[] FP;
}

/*preassembles the K matrix and "h" and "g" parts of the force vector*/
void FESolver::preAssembly() { TRACE_ME;
    /*loop over elements*/
    for (int e=0;e<n_elements;e++) {
        Tetra &tet = volume.elements[e];
        double ke[4][4];

        for (int a=0;a<4;a++)
            for (int b=0;b<4;b++) {
                ke[a][b] = 0;    /*reset*/

                /*perform quadrature*/
                for (int k=0;k<n_int;k++)
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) {
                            double nax[3],nbx[3];

                            getNax(nax,e,a);
                            getNax(nbx,e,b);

                            /*dot product*/
                            double dot=0;
                            for (int d=0;d<3;d++) dot+=nax[d]*nbx[d];
                            ke[a][b] += dot*detJ[e]*W[i]*W[j]*W[k];
                        }
            }

        /*we now have the ke matrix*/
        addKe(e,ke);

        /*force vector*/
        double fe[4];

        for (int a=0;a<4;a++) {
            /*second term int(na*h), always zero since support only h=0*/
            double fh=0;

            /*third term, -sum(kab*qb)*/
            double fg = 0;
            for (int b=0;b<4;b++) {
                int n = tet.con[b];
                double gb = g[n];
                fg-=ke[a][b]*gb;
            }

            /*combine*/
            fe[a] = fh + fg;
        }
        addFe(F0, e,fe);
    }  /*end of element*/


}

/*computes "ff" part of F*/
void FESolver::buildF1Vector(double *ion_den) { TRACE_ME;
    double *f = new double[neq];
    /*start by computing the RHS term on all unknown nodes*/
    for (int n=0;n<n_nodes;n++) {
        int A = ID[n];
        if (A<0) continue;    /*skip known nodes*/
        f[A] = (QE/EPS0)*(ion_den[n]+n0*exp((d[A]-phi0)/kTe));
    }

    /*loop over elements*/
    for (int e=0;e<n_elements;e++) {
        double fe[4];
        for (int a=0;a<4;a++) {
            /*first term is int(na*f), set to zero for now*/
            double ff=0;
            int A = LM[e][a];
            if (A>=0)    /*if unknown node*/ {
                /*perform quadrature*/
                for (int k=0;k<n_int;k++)
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) {
                            /*change of limits*/
                            double xi = 0.5*(l[i]+1);
                            double eta = 0.5*(l[j]+1);
                            double zeta = 0.5*(l[k]+1);

                            double Na=evalNa(a,xi,eta,zeta);
                            ff += f[A]*Na*detJ[e]*W[i]*W[j]*W[k];
                        }
                ff*=(1.0/8.0);    /*change of limits*/
                fe[a] = ff;
            }
        }

        addFe(F1, e,fe);
    }

    delete[] f;
}

/*simple Gauss-Seidel solver for A*x=b*/
void FESolver::solveLinear(double **A, double *x, double *b) { TRACE_ME;
    int it;
    const double tol=1e-4;
    double L2;

    for (int u=0;u<neq;u++)
        if (fabs(A[u][u])<1e-12)
            std::cout<<"Zero diagonal on "<<u<<std::endl;

    bool converged=false;
    for (it=0;it<10000;it++) {
        for (int u=0;u<neq;u++) {
            /*skip over unused nodes*/
            if (fabs(A[u][u])<1e-12) continue;

            double sum=0;
            for (int v=0;v<neq;v++) {
                if (u==v) continue;
                sum+=A[u][v]*x[v];
            }
            x[u] = (b[u]-sum)/A[u][u];
        }

        /*periodically compute residue*/
        if (it%25==0) {
            double L=0;
            for (int u=0;u<neq;u++) {
                double sum=0;
                for (int v=0;v<neq;v++)
                    sum+=A[u][v]*x[v];
                double r=b[u]-sum;
                L+=r*r;
            }
            L2 = sqrt(L)/neq;
            if (L2<tol) {converged=true; break;}
        }
    }

    if (!converged) std::cerr<<" GS failed to converge in "<<it<<" iterations, "<<std::setprecision(3)<<": L2="<<L2<<std::endl;
}

/*wrapper for solving the non-linear Poisson's equation*/
void FESolver::computePhi(double *ion_den) { TRACE_ME;

    /*solve the system*/
    solveNonLinear(ion_den);

    /*combine d and g to phi*/
    for (int n=0;n<n_nodes;n++) {
        /*zero on non-g nodes*/
        uh[n] = g[n];

        /*is this a non-g node?*/
        int A=ID[n];
        if (A>=0)
            uh[n] += d[A];
    }
}

/*updates electric field*/
void FESolver::updateEf() { TRACE_ME;
    /*interpolate electric field*/
    for (int e=0;e<n_elements;e++) {
        Tetra &tet = volume.elements[e];
        for (int d=0;d<3;d++) ef[e][d]=0;

        for (int a=0;a<4;a++) {
            int A = tet.con[a];
            double nx[3];
            getNax(nx,e,a);
            /*minus sign since negative gradient*/
            for (int d=0;d<3;d++) ef[e][d]-=nx[d]*uh[A];
        }
    }
}

