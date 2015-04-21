#ifndef _PARAMETER_H_
#define _PARAMETER_H_

/* =============================================
             Physical Parameters
============================================= */
#define BONDLENGTH 5 // bond length (nm)
#define PERSISTENCE_LENGTH 3
#define OMEGA 3
#define SLIT_HEIGHT 5

/* =============================================
         Program Related Parameters
============================================= */
#define MAXLENGTH 200 
#define MINLENGTH 10
#define PURANGE 10 // Metropolis probability random number range
#define STABLE_EPSILON 0.03 // range of Energy to decide if result converges
#define ENERGY_ARRAY_LENGTH 5 // length of energyArray. If all nums in this array are restricted to +/-STABLE_EPSILON range, then chain is marked as stable.

#endif