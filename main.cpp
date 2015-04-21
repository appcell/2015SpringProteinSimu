#include <iostream>
#include <cmath>
#include <random>
#include <time.h>
#include <stdlib.h>
#include "const.h"
#include "parameter.h"
#include "simplefunc.h"
#include "chain.h"

using namespace std;

void BeadCopy(Bead *_this, Bead b)  
{  
	_this->phi = b.phi;
	_this->theta = b.theta;
	_this->x = b.x;
	_this->y = b.y;
	_this->z = b.z;
} 

// SUCCESS!!!!!
Bead GenerateRandomBead(double bondlength, Bead &prev, bool flag) {
	// flag == 1, prev is prev; flag != 1, prev is next
	Bead newBead = {0, 0, 0, 0, 0};
	int i = (flag == true)? 1 :(-1);
	Bead *pointer = (flag == true)? &newBead:&prev;
	pointer->phi = randomf(-PI, PI);
	pointer->theta = randomf(-PI/2, PI/2);
	newBead.x = i * bondlength*cosl(pointer->theta)*cosl(pointer->phi) + prev.x;
	newBead.y = i * bondlength*cosl(pointer->theta)*sinl(pointer->phi) + prev.y;
	newBead.z = i * bondlength*sinl(pointer->theta) + prev.z;
	return newBead;
}

// SUCCESS!!!!!
Bead RectangularToSpherical(double deltaX, double deltaY, double deltaZ, double bondlength, Bead &prev, bool flag) {
	Bead res = {0, 0, 0, 0, 0}, * pointer = (flag == true)? &res:&prev;
	int i = (flag == true)? 1 : (-1);
	pointer->theta = asinf(deltaZ/bondlength);
	if(equal(deltaX,0)) 
		if(deltaY>0)
			pointer->phi = PI/2;
		else
			pointer->phi = -PI/2;
	else if(equal(deltaY,0))
		if(deltaX>0) 
			pointer->phi = 0;
		else 
			pointer->phi = PI;
	else {
		pointer->phi = atanf(deltaY/deltaX);
		if(deltaY<0 && deltaX<0)
			pointer->phi += PI;
	}
	res.x = i * deltaX + prev.x;
	res.y = i * deltaY + prev.y;
	res.z = i * deltaZ + prev.z;

	return res;
}

int main() {
	Chain tmp(400);
	tmp.MonteCarloChain();
	cout<<"\n"<<tmp.GetEnergy();

	return 0;
}