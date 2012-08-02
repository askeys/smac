/*=================================================================
 *
 * ICP algorithm implemented in c++
 *
 * written by
 *
 * Per Bergström 2007-10-09
 *
 * email: per.bergstrom@ltu.se
 *
 * Uses kd-tree written by Guy Shechter
 * http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4586&objectType=file
 *
 *=================================================================*/

#include "kdtree_common.h"
//#include "kdtree_common.cc"

void icp(
double *trpr,
double *ttpr,
double *modelz,
unsigned int nmodelz,
double *dataz,
double *qltyz,
unsigned int ndataz,
unsigned int *randvecz,
unsigned int nrandvecz,
unsigned int nrandz,
unsigned int iimax,
Tree* tree
);
