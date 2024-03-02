#include "chipmunk/cpRobust.h"


cpBool
cpCheckPointGreater(const cpVect a, const cpVect b, const cpVect c)
{
	return fix14_mul((b.y - a.y), (a.x + b.x - (c.x << 1))) > fix14_mul((b.x - a.x), (a.y + b.y - (c.y << 1)));
}

cpBool
cpCheckAxis(cpVect v0, cpVect v1, cpVect p, cpVect n){
	return cpvdot(p, n) <= cpfmax(cpvdot(v0, n), cpvdot(v1, n));
}
