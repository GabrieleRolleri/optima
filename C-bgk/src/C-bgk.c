/*  Lattice Boltzmann sample, written in C
 *
 *  Copyright (C) 2006 Jonas Latt
 *  Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland 
 *  E-mail: Jonas.Latt@cui.unige.ch
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public 
 *  License along with this program; if not, write to the Free 
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

/* unsteady.c:
 * This example examines an unsteady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, where as the outlet implements an outflow condition:
 * grad_x u = 0. At Reynolds numbers around 100, an unstable periodic pattern is
 * created, the Karman vortex street.
 */

#include "lb.h"
#include "boundaries.h"
#include <stdio.h>
#include <stdlib.h>

  // These constants define the flow geometry and are commented in
  //   the function setConstants()
int lx, ly;
int obst_x, obst_y, obst_r;
float uMax, Re, nu, omega;
int maxT, tSave;

  // The dynamics that are to be specified on different regions of 
  //   the flow: LBGK in the bulk, bounce-back on the cylinder, 
  //   regularized boundary condition on the four domain boundaries.
  //   Alternatively, the Zou/He boundary condition can be used:
  //   replace upperRegularized by upperZouHe etc.
Dynamics       bulkDynamics         = { &bgk, (void*)&omega };
Dynamics       bounceBackDynamics   = { &bounceBack, 0 };

  // The velocity to be imposed on the upper/lower boundary: u=0
VelocityBCData zeroVelocityBoundary = { &bulkDynamics, 0., 0. };

Dynamics upperBoundary = { &upperRegularized,
                           (void*)&zeroVelocityBoundary };
Dynamics lowerBoundary = { &lowerRegularized,
                           (void*)&zeroVelocityBoundary };
Dynamics* leftBoundary;   // Those two objects are initialized in the
Dynamics* rightBoundary;  //   function iniData()

  // These arrays contain the velocities that are to be imposed on 
  //   the inlet (poiseuilleBoundary) and the outlet 
  //   (pressureBoundary) of the channel
VelocityBCData* poiseuilleBoundary;
PressureBCData* pressureBoundary;

  // The main object containing the simulation
Simulation sim;

void setConstants() {
    lx       = 250;     // channel length
    ly       = 50;      // channel height
    obst_x = lx/5;      // position of the cylinder; the cylinder is
    obst_y = ly/2;      // offset from the center to break symmetry
    obst_r = ly/10+1;   // radius of the cylinder

    uMax  = 0.02;       // maximum velocity of the Poiseuille inflow
    Re    = 100;        // Reynolds number
    nu    = uMax * 2.*obst_r / Re;  // kinematic fluid viscosity
    omega = 1. / (3*nu+1./2.);      // relaxation parameter

    maxT   = 100000;       // total number of iterations
    tSave  = 100;          // frequency of periodic saves to disk
}

  // Memory allocation and default initialisation of the simulation 
  //   and the left/right boundaries
void iniData() {
    int iX, iY;
    poiseuilleBoundary =
        (VelocityBCData*) calloc(ly+2, sizeof(VelocityBCData));
    pressureBoundary =
        (PressureBCData*) calloc(ly+2, sizeof(PressureBCData));

    leftBoundary  = (Dynamics*) calloc(ly+2, sizeof(Dynamics));
    rightBoundary = (Dynamics*) calloc(ly+2, sizeof(Dynamics));

    for (iY=2; iY<=ly-1; ++iY) {
        poiseuilleBoundary[iY].bulkDynamics   = &bulkDynamics;
        poiseuilleBoundary[iY].uy             = 0.;
        pressureBoundary[iY].bulkDynamics = &bulkDynamics;
        pressureBoundary[iY].rho          = 1.;
        pressureBoundary[iY].uPar         = 0.;

        leftBoundary[iY].dynamicsFun = &leftRegularized;
        leftBoundary[iY].selfData
            = (void*)&poiseuilleBoundary[iY];
        rightBoundary[iY].dynamicsFun = &rightPressureRegularized;
        rightBoundary[iY].selfData
            = (void*)&pressureBoundary[iY];
    }
}

  // De-allocation of the memory
void freeData() {
    free(rightBoundary);
    free(leftBoundary);
    free(poiseuilleBoundary);
    free(pressureBoundary);
}

  // compute parabolic Poiseuille profile
float computePoiseuille(int iY) {
    float y = (float)(iY-1);
    float L = (float)(ly-1);
    return 4.*uMax / (L*L) * (L*y-y*y);
}

  // Specify the geometry of the simulation
void iniGeometry() {
    int iX, iY;
    for(iX=1; iX<=lx; ++iX) {
        for(iY=1; iY<=ly; ++iY) {
              // All bulk nodes are initialized at equilibrium with constant
              // density and a velocity determined by a y-dependend Poiseuille
              // profile.
            float uPoiseuille = computePoiseuille(iY);
            iniEquilibrium(&sim.lattice[iX][iY], 1., uPoiseuille, 0.);
              // on the obstacle, set bounce back dynamics
            if ( (iX-obst_x)*(iX-obst_x) +
                 (iY-obst_y)*(iY-obst_y) <= obst_r*obst_r )
            {
                setDynamics(&sim, iX, iY, &bounceBackDynamics);
            }
              // elsewhere, use lbgk dynamics
            else {
                setDynamics(&sim, iX, iY, &bulkDynamics);
            }
        }
    }

      // upper and lower boundary: u=0
    for (iX=1; iX<=lx; ++iX) {
        setDynamics(&sim, iX, 1, &lowerBoundary);
        setDynamics(&sim, iX, ly, &upperBoundary);
    }

      // left boundary: Poiseuille profile, constant through time
      // right boundary: initially Poiseuille profile, then outlet
      //   condition grad_x u = 0
    for (iY=2; iY<=ly-1; ++iY) {
        float uPoiseuille = computePoiseuille(iY);
        poiseuilleBoundary[iY].ux   = uPoiseuille;
        setDynamics(&sim, 1, iY, &leftBoundary[iY]);
        setDynamics(&sim, lx, iY, &rightBoundary[iY]);
    }
}

  // Compute a second order extrapolation on the right boundary to
  // ensure a zero-gradient boundary condition on the pressure.
  // This must be recomputed at every time step. The velocity is
  // constrained to be perpendicular to the outflow surface.
void updateZeroGradientBoundary() {
    int iY;
    float rho1, ux1, uy1, rho2, ux2, uy2;
    for (iY=2; iY<=ly-1; ++iY) {
        computeMacros(sim.lattice[lx-1][iY].fPop, &rho1, &ux1,&uy1);
        computeMacros(sim.lattice[lx-2][iY].fPop, &rho2, &ux2,&uy2);
        pressureBoundary[iY].rho = 4./3.*rho1 - 1./3.*rho2;
        pressureBoundary[iY].uPar = 0.;
    }
}

static const float t[9] = { 4./9., 1./9., 1./9., 1./9., 1./9.,
                             1./36., 1./36., 1./36., 1./36. };
  // lattice velocities
static const int c[9][2] = {
    {0,0},
    {1,0}, {0,1}, {-1,0}, {0,-1},
    {1,1}, {-1,1}, {-1,-1}, {1,-1}
};


int main(int argc, char *argv[]) {
      // initialisation of a lx * ly simulation
    setConstants();
    printf("\nlx=%d, ly=%d, omega=%f\n\n", lx, ly, omega);

    iniData();
    constructSim(&sim, lx, ly);
    iniGeometry();

      // the main loop over time steps
    int iT;
    for (iT=0; iT<maxT; ++iT) {
        if (iT%1000==0) {
            printf("t=%d\n", iT);
        }

          // on the right boundary, outlet condition grad_x u = 0
        updateZeroGradientBoundary();

        collide(&sim);

          // the data are written to disk after collision, to be
          //   that the macroscopic variables are computed 
          //   correctly on the boundaries
        if (iT%tSave==0 && iT>0) {
            printf("t=%d\n", iT);
            saveVel(&sim, "vel.dat");
        }

        propagate(&sim);

          // By default: periodic boundary conditions. In this case,
          //   this is important, because the upper and lower
          //   boundaries are horizontally periodic, so that no
          //   special corner nodes are needed.
        makePeriodic(&sim);
    }

    destructSim(&sim);
    freeData();
}
