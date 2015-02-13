#ifndef RANSACMATCHING_H_
#define RANSACMATCHING_H_

#include <flann/flann.hpp>
#include "obcore/math/linalg/linalg.h"
#include "omp.h"

namespace obvious
{

/**
 * @class RansacMatching
 * @brief Matching algorithm with RANSAC scheme
 * @author Stefan May
 **/
class RansacMatching
{
public:
  /**
   * Constructor
   * @param trials number of trials / matching guesses
   * @param epsThresh threshold for rating good matches
   * @param phiMax maximum rotation
   * @param sizeControlSet approximate set of control set
   */
  RansacMatching(unsigned int trials = 50, double epsThresh = 0.03, unsigned int sizeControlSet = 180, bool clipPeripheralArea = false);

  /**
   * Destructor
   */
  virtual ~RansacMatching();

  /**
   * Matching method
   * @param M model
   * @param S scene
   * @param maskM Mask with fields = true for valid points in the scene S
   * @param maskS Mas with fields = true for valid points in the model M
   * @param validPointsM Number of valid points in the model M
   * @param validPointsS Number of valid points in the scene S
   * @param phiMax The algorithm will not obtain rotations between two scans bigger than this.
   * @param resolution Resolution of the laser scan (angular distance between obtained measurements)
   * @return 3x3 registration matrix
   */
  obvious::Matrix match(obvious::Matrix* M, obvious::Matrix* S, bool* maskM, bool* maskS, unsigned int validPointsM, unsigned int validPointsS, double phiMax, double resolution);

private:

  // squared distance threshold
  double _epsSqr;

  // number of trials
  unsigned int _trials;

  // approximate control set
  unsigned int _sizeControlSet;

  // clip peripheral area of laser (non-overlapping area)
  bool _clipPeripheralArea;

  // tree for accelerating NN search
  flann::Index<flann::L2<double> >* _index;
  flann::Matrix<double>* _model;

};

}

#endif /* RANSACMATCHING_H_ */
