#ifndef NDT_H_
#define NDT_H_

#include <iostream>
#include <string.h>
using namespace std;

#include "obcore/base/CartesianCloud.h"
#include "obcore/base/System.h"
#include "obvision/registration/Registration.h"
#include "obcore/math/linalg/linalg.h"
#include "obvision/registration/Registration.h"

using namespace obvious;

namespace obvious
{

/**
 * NDT return states
 */
enum EnumNdtState { NDT_IDLE 			= 0,
  NDT_PROCESSING 		= 1,
  NDT_NOTMATCHABLE 	= 2,
  NDT_MAXITERATIONS 	= 3,
  NDT_TIMEELAPSED 	= 4,
  NDT_SUCCESS 		= 5,
  NDT_CONVERGED   = 6,
  NDT_ERROR			= 7 };


struct NdtCell
{
  vector<double*> coords;
  double* centroid;
  Matrix* cov;
  Matrix* cov_inv;
  bool isOccupied(){ return (coords.size()>=5); };
};

/**
 * @class Ndt
 * @brief Represents the normal distribution transform
 * @author Stefan May
 **/
class Ndt: public obvious::Registration
{
public:
  /**
   * Standard constructor
   */
  Ndt(int minX, int maxX, int minY, int maxY);

  /**
   * Destructor
   */
  ~Ndt();

  /**
   * Sample model point cloud to NDT space
   * @param coords model coordinates
   * @param probability probability of coordinates of being sampled (range [0.0 1.0])
   */
  void setModel(Matrix* coords, double probability=1.0);

  /**
   * Copy scene to internal buffer
   * @param coords scene coordinates
   * @param probability probability of coordinates of being sampled (range [0.0 1.0])
   */
  void setScene(Matrix* coords, double probability=1.0);

  /**
   * Reset state of NDT algorithm
   */
  void reset();

  /**
   * Start iteration
   * @param rms return value of RMS error
   * @param iterations return value of performed iterations
   * @param Tinit apply initial transformation before iteration
   * @return  processing state
   */
  EnumState iterate(double* rms, unsigned int* iterations, Matrix* Tinit=NULL);

private:

  int _minX;
  int _maxX;
  int _minY;
  int _maxY;

  /**
   * size of internal scene buffer
   */
  unsigned int _sizeSceneBuf;

  NdtCell** _model;

  /**
   * the scene
   */
  double** _scene;

  /**
   * local copy of scene
   */
  double** _sceneTmp;

  /**
   * size of scene
   */
  unsigned int _sizeScene;

  double _d1;

  double _d2;

};

}

#endif /*NDT_H_*/
