#include <iostream>
#include <cmath>
#include <random>
#include <time.h>
#include "const.h"
#include "parameter.h"
#include "simplefunc.h"
#include "chain.h"

using namespace std;

Chain::Chain()
{
	length = chain.size();
	stable = false;
	energy = 0;
	knotLength = 0;
}

// SUCCESS!!!!!
Chain::Chain(int length) {
	Bead temp = {0, 0, 0, 0, 0};
	Bead temp2;
	for (int i=0; i<length; i++) {
		chain.push_back(temp);
		temp2 = GenerateRandomBead(BONDLENGTH, temp, true);
		BeadCopy(&temp, temp2);
	}
	this->length = length;
	this->energy = GetEnergy();
	knotLength = 0;
}

Chain::~Chain()
{
}

bool Chain::CrankshaftMove() {
	if(length >= 3) {
		int i = random(length - 3); // here i is start point of three points
		Bead tempa, tempb, tempc;
		double oldEnergy = energy;
		tempa = this->chain[i];
		tempb = this->chain[i+1];
		tempc = this->chain[i+2];
		RotateBead(i, BONDLENGTH, randomf(0, 2*PI));
		if( IsChangeAcceptable(oldEnergy, ModifyEnergy(0, i+1, tempb, (this->chain[i+1]), &tempc, &(this->chain[i+2]) ))) {
			return true;
		}
		else {
			BeadCopy(&(this->chain[i]), tempa);
			BeadCopy(&(this->chain[i+1]), tempb);
			BeadCopy(&(this->chain[i+2]), tempc);
			this->energy = oldEnergy;
			return false;
		}
	}
	return false;
}

bool Chain::ReptationMove() {
	int i = random(1);
	double oldEnergy = energy;
	Bead temp;
	if(i == 0) {
		// head move to tail
		temp = chain.front();
		chain.push_back(GenerateRandomBead(BONDLENGTH,chain.back(),true));
		chain.erase(chain.begin());
	} else {
		// tail move to head
		temp = chain.back();
		chain.insert(chain.begin(), GenerateRandomBead(BONDLENGTH, chain.front(), false));
		chain.pop_back();
	}
	if(IsChangeAcceptable(oldEnergy, GetEnergy())) {
		return true;
	} else {
		if (i == 0) {
			chain.pop_back();
			chain.insert(chain.begin(), temp);
		} else {
			chain.erase(chain.begin());
			chain.push_back(temp);
		}
		energy = oldEnergy;
		return false;
	}
}

double Chain::GetEnergy() {
	double res = 0, freeEnergy = 0, bendEnergy = 0, rParallel2 = 0;
	double contourLength = (length - 1) * BONDLENGTH;
	double ex = 0, ex2 = 0, ey = 0, ey2 = 0, ez = 0, ez2 = 0, theta = 0;
	for(int i=0; i<length; i++){
		ex += chain[i].x, ey += chain[i].y, ez += chain[i].z;
		ex2 += chain[i].x*chain[i].x, ey2 += chain[i].y*chain[i].y, ez2 += chain[i].z*chain[i].z;
		if (i < length - 1 && i > 0) {
			theta = GetVectorAngle(chain[i-1], chain[i], chain[i+1]);
			bendEnergy += theta * theta;
		}
	}
	bendEnergyStorage = bendEnergy;
	rParallel2 = (ex2 + ey2 + ez2 - ex*ex - ey*ey - ez*ez)/length;
	this->ex2 = ex2, this->ey2 = ey2, this->ez2 = ez2, this->ex = ex, this->ey = ey, this->ez = ez;
	freeEnergy = KBT * ((rParallel2) / (contourLength * PERSISTENCE_LENGTH) 
		+ (OMEGA * contourLength * contourLength) / (rParallel2 * SLIT_HEIGHT));
	bendEnergy = 0.5 * bendEnergy * PERSISTENCE_LENGTH * KBT / BONDLENGTH;
	res = freeEnergy + bendEnergy;
	energy = res;
	return res;
}

double Chain::ModifyEnergy(int modType, int modPos, Bead oldBead, Bead newBead, Bead *oldBead2, Bead *newBead2) {
	double res = 0, freeEnergy = 0, bendEnergy = 0, rParallel2 = 0, theta = 0, bendEnergyTemp = 0;
	double contourLength = (length - 1) * BONDLENGTH;
    // minus old variance
    this->ex2 -= (oldBead.x * oldBead.x), this->ey2 -= (oldBead.y * oldBead.y), this->ez2 -= (oldBead.z * oldBead.z);
    this->ex -= oldBead.x, this->ey -= oldBead.y, this->ez -= oldBead.z;
    // plus new variance
    this->ex2 += (newBead.x * newBead.x), this->ey2 += (newBead.y * newBead.y), this->ez2 += (newBead.z * newBead.z);
    this->ex += newBead.x, this->ey += newBead.y, this->ez += newBead.z;
    // calculate free energy
    rParallel2 = (this->ex2 + this->ey2 + this->ez2 - this->ex*this->ex - this->ey*this->ey - this->ez*this->ez) / length;
    freeEnergy = KBT * ((rParallel2) / (contourLength * PERSISTENCE_LENGTH) 
                + (OMEGA * contourLength * contourLength) / (rParallel2 * SLIT_HEIGHT));

	switch (modType) {
	case 0: {
		// crankshaft move
		// old1: [modPos], old2: [modPos+1]
		if (modPos - 1 > 0) {
			theta = GetVectorAngle(chain[modPos-2], chain[modPos-1], oldBead);
			bendEnergyTemp -= theta * theta; // old
			theta = GetVectorAngle(chain[modPos-2], chain[modPos-1], newBead);
			bendEnergyTemp += theta * theta; // new
		}
		if (modPos + 1 < length - 1) {
			theta = GetVectorAngle(oldBead, *oldBead2, chain[modPos+2]);
			bendEnergyTemp -= theta * theta;
			theta = GetVectorAngle(newBead, *newBead2, chain[modPos+2]);
			bendEnergyTemp += theta * theta;
		}
		bendEnergyStorage += bendEnergyTemp, bendEnergyTemp = 0;
		bendEnergy = 0.5 * bendEnergyStorage * PERSISTENCE_LENGTH * KBT / BONDLENGTH;
		energy = bendEnergy + freeEnergy;

		return energy;
			}

	case 1:{
		// reptation move, head move to tail. -head + new to tail. old is old head, new is new tail. angle - (0<->1), + (length-2<->length-1)()
		if (length >= 3) {
			theta = GetVectorAngle(oldBead, chain[0], chain[1]);
			bendEnergyTemp -= theta * theta;
			theta = GetVectorAngle(chain[length-3], chain[length-2], chain[length-1]);
			bendEnergyTemp += theta * theta;
		}
		bendEnergyStorage += bendEnergyTemp;
		bendEnergy = 0.5 * bendEnergyStorage * PERSISTENCE_LENGTH * KBT / BONDLENGTH;
		energy = bendEnergy + freeEnergy;

		return energy;
		   }

	case 2: {
		// reptation move, tail move to head. -tail + new to chain[1]. old is old tail. new is chain[0]. angle - (old <-> length-1), + (1<->2)
		if (length >= 3) {
			theta = GetVectorAngle(chain[length-2], chain[length-1], oldBead);
			bendEnergyTemp -= theta * theta;
			theta = GetVectorAngle(newBead, chain[1], chain[2]);
			bendEnergyTemp += theta * theta;
		}
		bendEnergyStorage += bendEnergyTemp;
		bendEnergy = 0.5 * bendEnergyStorage * PERSISTENCE_LENGTH * KBT / BONDLENGTH;
		energy = bendEnergy + freeEnergy;

		return energy;
			}

	default: return energy;
	}
}

double Chain::GetVectorAngle(Bead &prev, Bead &curr, Bead &next) {
	// get angle around i
	double theta = 0;
	double x1 = curr.x - prev.x,
		y1 = curr.y - prev.y,
		z1 = curr.z - prev.z,
		x2 = next.x - curr.x,
		y2 = next.y - curr.y,
		z2 = next.z - curr.z;
	theta = acosl((x1 * x2 + y1 * y2 + z1 * z2) / (BONDLENGTH * BONDLENGTH));
	return theta;
}

bool Chain::IsChainStable(double energyArray[ENERGY_ARRAY_LENGTH]) {
	for (int i = 1; i<ENERGY_ARRAY_LENGTH; i++) {
		if(fabs(energyArray[i] - energyArray[i-1]) > STABLE_EPSILON)
			return false;
	}
	return true;
}

// ====== Parameters wait to be tuned ======
bool Chain::IsChangeAcceptable (double oldEnergy, double newEnergy) {
	double dE = newEnergy - oldEnergy;
	if(dE < 0) {
		return true;
	} else {
		double probability = expl((0 - dE)/KBT);
		double i = randomf(0-PURANGE, PURANGE);
		if(i<probability) return true;
		else return false;
	}
}

// SUCCESS!!!!!
bool Chain::RotateBead(int i, double bondlength, double rotateAngle) {
    Bead *a = &(this->chain[i]), *b = &(this->chain[i+1]), *c = &(this->chain[i+2]);
    if (rotateAngle >= 2 * PI) return false;
	double x1 = b->x - a->x, 
           y1 = b->y - a->y, 
           z1 = b->z - a->z, 
           x2 = c->x - b->x, 
           y2 = c->y - b->y, 
           z2 = c->z - b->z;
    x2 += x1, y2 += y1, z2 += z1;
    double nx = x2 / 2.0, ny = y2 / 2.0, nz = z2 / 2.0; // get vector n
    double len = sqrtf(nx * nx + ny * ny + nz * nz); 
    double x, y, z;
    nx /= len, ny /= len, nz /= len;   
    x =  x1 * (cosl(rotateAngle) + nx * nx * (1 - cosl(rotateAngle))) +    
     y1 * (nx * ny * (1 - cosl(rotateAngle)) - nz * sinl(rotateAngle)) + 
    z1 * (nx * nz * (1 - cosl(rotateAngle)) + ny * sinl(rotateAngle));
    y = x1 * (nx * ny * (1 - cosl(rotateAngle)) + nz * sinl(rotateAngle)) +  
    y1 * (ny * ny * (1 - cosl(rotateAngle)) + cosl(rotateAngle)) + 
    z1 * (ny * nz * (1 - cosl(rotateAngle)) - nx * sinl(rotateAngle));
    z = x1 * (nx * nz * (1 - cosl(rotateAngle)) - ny * sinl(rotateAngle)) + 
    y1 * (ny * nz * (1 -cosl(rotateAngle)) + nx * sinl(rotateAngle)) + 
    z1 * (nz * nz * (1 - cosl(rotateAngle)) + cosl(rotateAngle));//transform by matrix

	BeadCopy(b, RectangularToSpherical(x, y, z, bondlength, chain[i], true)), 
		BeadCopy(c, RectangularToSpherical((x2-x), (y2-y), (z2-z), bondlength,chain[i+1], true));

    /*cout<<x<<"  "<<y<<"  "<<z<<"\n";
    cout<<bondlength*cosl(b->theta)*cosl(b->phi)<<"  "<<bondlength*cosl(b->theta)*sinl(b->phi)<<"  "<<bondlength*sinl(b->theta)<<"\n\n";
    cout<<x2<<"  "<<y2<<"  "<<z2<<"\n";
    cout<<x+bondlength*cosl(c->theta)*cosl(c->phi)<<"  "<<y+bondlength*cosl(c->theta)*sinl(c->phi)<<"  "<<z+bondlength*sinl(c->theta)<<"\n";*/
    return true;
}

// SUCCESS!!!!!
void Chain::PrintChain() {
	cout<<"Chain Info:\n"<<"Length: "<<length<<"\nFree Energy: "<<energy<<"\n";
	for (int i=0; i<length; i++) {
		cout<<"NO."<<i<<",    Phi: "<<chain[i].phi<<",    Theta: "<<chain[i].theta<<"\n";
		cout<<"         "<<"x: "<<chain[i].x<<",    y: "<<chain[i].y<<",    z: "<<chain[i].z<<"\n";
	}
}

// ====== Unfinished ======
void Chain::MonteCarloChain() {
	double oldEnergy, newEnergy;
	double energyArray[ENERGY_ARRAY_LENGTH];
	int i=0;
	do {
		i++;
		energyArray[i%ENERGY_ARRAY_LENGTH] = energy;
		CrankshaftMove();
		ReptationMove();
	}while(!IsChainStable(energyArray));

	cout<<i;
}

// ====== Unfinished ======
double Chain::GetKnotLength() {
	double res = 0;

	knotLength = res;
	return res;
}