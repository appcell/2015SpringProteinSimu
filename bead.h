#ifndef _BEAD_H_
#define _BEAD_H_

typedef struct
{
	double phi;
	double theta;
	double x;
	double y;
	double z;
}Bead;
void BeadCopy(Bead *_this, Bead b);
Bead GenerateRandomBead(double bondlength, Bead &prev, bool flag);
Bead RectangularToSpherical(double x, double y, double z, double bondlength, Bead &prev, bool flag);

#endif