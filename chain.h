#ifndef _CHAIN_H_
#define _CHAIN_H_

#include <vector>
#include "bead.h"
using namespace std;

class Chain
{
public:
    Chain();
    Chain(int length);
    ~Chain();
    bool CrankshaftMove();
    bool ReptationMove();
    double GetEnergy();
	double ModifyEnergy(int modType, int modPos, Bead oldBead, Bead newBead, Bead *oldBead2 = NULL, Bead *newBead2 = NULL);
	double GetVectorAngle(Bead &prev, Bead &curr, Bead &next);
    bool IsChainStable(double energyArray[5]);
    bool IsChangeAcceptable(double oldEnergy, double newEnergy);
    bool RotateBead(int i, double bondlength, double rotateAngle);
    void PrintChain();
    void MonteCarloChain();
    double GetKnotLength();

private:
    vector<Bead> chain;
    int length;
    double energy;
    bool stable;
    double knotLength;
	double ex2, ey2, ez2, ex, ey, ez;
	double bendEnergyStorage;
};

#endif