/*
 * Registration.cpp
 *
 *  Created on: 08.12.2014
 *      Author: Markus Kühn
 */
#include "Registration.h"

#include "obcore/base/tools.h"
#include "obcore/base/Timer.h"
#include "obcore/math/mathbase.h"

namespace obvious {

const char* g_states[] = { "NDT_IDLE", "NDT_PROCESSING", "NDT_NOTMATCHABLE",
		"NDT_MAXITERATIONS", "NDT_TIMEELAPSED", "NDT_SUCCESS", "NDT_ERROR" };

const char* Registration::state2char(EnumState eState) {
	return g_states[eState];
}
;

virtual void Registration::reset() {
	_Tfinal4x4->setIdentity();
	if (_sceneTmp)
		System<double>::copy(_sizeScene, _dim, _scene, _sceneTmp);
}

void Registration::setMaxIterations(unsigned int iterations)
{
  Registration::_maxIterations = iterations;
}

unsigned int Registration::getMaxIterations()
{
  return Registration::_maxIterations;
}

Matrix Registration::getFinalTransformation4x4()
{
  Matrix T = *_Tfinal4x4;
  return T;
}

Matrix Registration::getFinalTransformation()
{
  Matrix T(_dim+1, _dim+1);
   for(int r=0; r<_dim; r++)
   {
      for(int c=0; c<_dim; c++)
      {
         T(r,c) = (*_Tfinal4x4)(r,c);
      }
      T(r,_dim) = (*_Tfinal4x4)(r,3);
   }

   for(int c=0; c<_dim; c++)
      T(_dim,c) = 0;

   T(_dim,_dim) = 1;

  return T;
}

Matrix Registration::getLastTransformation()
{
  Matrix T = *_Tlast;
  return T;
}


void Registration::applyTransformation(double** data, unsigned int size, unsigned int dim, Matrix* T)
{
  // Apply rotation
  Matrix R(*T, 0, 0, dim, dim);
  Matrix::multiply(R, *data, size, dim);

  // Apply translation
  if(_dim < 3)
  {

#pragma omp parallel
{
#pragma omp for
  for(unsigned int i=0; i<size; i++)
  {
    data[i][0] += (*T)(0,3);
    data[i][1] += (*T)(1,3);
  }
}

  }
  else
  {

#pragma omp parallel
{
#pragma omp for
  for(unsigned int i=0; i<size; i++)
  {
    data[i][0] += (*T)(0,3);
    data[i][1] += (*T)(1,3);
    data[i][2] += (*T)(2,3);
  }
}

  } // end if

}


void Registration::checkMemory(unsigned int rows, unsigned int cols, unsigned int &memsize, double** &mem)
{
  // first instantiation of buffer
  if(mem == NULL)
  {
    memsize = rows;
    System<double>::allocate(rows, cols, mem);
  }
  // resize buffer, if needed
  if(rows > memsize)
  {
    memsize = rows;
    if(mem != NULL) System<double>::deallocate(mem);
    System<double>::allocate(rows, cols, mem);
  }
}



}

