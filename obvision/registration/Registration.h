/*
 * Registration.h
 *
 *  Created on: 08.12.2014
 *      Author: mkuehn
 */

#ifndef REGISTRATION_H_
#define REGISTRATION_H_

#include <iostream>
#include <string.h>
using namespace std;

#include "obcore/base/CartesianCloud.h"
#include "obcore/base/System.h"

#include "obcore/math/linalg/linalg.h"

using namespace obvious;

namespace obvious {
/**
 * Registration return states
 */
enum EnumState {
	IDLE = 0,
	PROCESSING = 1,
	NOTMATCHABLE = 2,
	MAXITERATIONS = 3,
	TIMEELAPSED = 4,
	SUCCESS = 5,
	CONVERGED = 6,
	ERROR = 7
};

/**
 * @class Registration
 * @brief This class represents a super class for registration algorithms
 * @author Markus Kühn
 **/
class Registration {
public:
	/**
	 * Standard constructor
	 * @param assigner pair assignment strategy for point clouds
	 * @param estimator transformation estimator
	 * @param pDistanceThreshold thresholding strategy
	 */
	Registration(){};

	/**
	 * Destructor
	 */
	virtual ~Registration(){};

	/**
	 * convert enumeration to char*
	 * @param eState state enumeration
	 * @return state string
	 */
	static const char* state2char(EnumState eState);

	/**
	 * Copy model to internal buffer
	 * @param coords model coordinates, as tuples or triples
	 * @param normals model normals, as tuples or triples, may be NULL
	 * @param size number of points, i.e. coordinate triples
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
//	virtual void setModel(double* coords, const unsigned int size, double probability = 1.0) = 0;

	/**
	 * Convenience method extracting all points from cloud to double array
	 * @param coords model coordinatess
	 * @param normals model normals, may be NULL
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
virtual void setModel(Matrix* coords,Matrix* normals = NULL, double probability = 1.0) = 0;

	/**
	 * Copy scene to internal buffer
	 * @param coords scene coordinates, as tuples or triples
	 * @param normals scene normals, as tuples or triples, may be NULL
	 * @param size number of points, i.e. coordinate triples
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
//	virtual void setScene(double* coords, const unsigned int size, double probability = 1.0) = 0;

	/**
	 * Convenience method extracting data from cloud to double array
	 * @param coords scene coordinates
	 * @param normals scene normals, may be NULL
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
	virtual void setScene(Matrix* coords,Matrix* normals = NULL, double probability = 1.0) = 0;


	/**
	 * Abstract method implemented by derived classes for loading algorithm specific parameters
	 * @param file Path to file.
	 * @return -1 if loading the parameters was not successful.
	 */
	virtual int loadParametersFromXML(string filePath) = 0;


	/**
	 * Reset state of Registration algorithm (resets iteration count, estimated transformation and point pairs)
	 */
	virtual void reset() = 0;

	virtual void activateTrace() = 0;

	virtual void serializeTrace(char* folder, unsigned int delay) = 0;


	/**
	 * Set maximum number of iteration steps
	 * @param iterations maximum number of iteration steps
	 */
	void setMaxIterations(unsigned int iterations);

	/**
	 * Get maximum number of iteration steps
	 * @return maximum number of iteration steps
	 */
	unsigned int getMaxIterations();

//	/**
//	 * Perform one iteration
//	 * @param rms return value of RMS error
//	 * @param pairs return value of pair assignments, i.e. number of pairs
//	 * @return processing state
//	 */
//	virtual EnumState step() = 0;

	/**
	 * Start registration process
	 * @param Tinit Apply initial transformation before iteration.
	 * @param verbose Set true if information about the registration are required in the console.
	 * @return  processing state
	 */
	virtual EnumState align(Matrix* Tinit=NULL, bool verbose=false) = 0;

	/**
	 * Get final 4x4 rotation matrix determined through iteration
	 * @return final transformation matrix
	 */
	Matrix getFinalTransformation4x4();

	/**
	 * Get final rotation matrix determined through iteration
	 * @return final transformation matrix
	 */
	Matrix getFinalTransformation();

	/**
	 * Get last rotation matrix determined within the last iteration step
	 * @return last transformation matrix
	 */
	Matrix getLastTransformation();

protected:


	/**
	 * apply transformation to data array
	 * @param data 2D or 3D coordinates
	 * @param size number of points
	 * @param dim dimensionality
	 * @param T transformation matrix
	 */
	void applyTransformation(double** data, unsigned int size, unsigned int dim,
			Matrix* T);

	/**
	 * internal memory check routine
	 * @param rows row size of needed memory
	 * @param memsize row size of target memory
	 * @param mem target memory
	 */
	void checkMemory(unsigned int rows, unsigned int cols,
			unsigned int &memsize, double** &mem);

	/**
	 * Creates an one-dimensional bool array to randomly use points or not
	 * @param size Number of points.
	 * @param probability Probability a point is used after masking a pointcloud
	 * @return mask
	 */
	static bool* createSubsamplingMask(unsigned int* size, double probability);


	/**
	 * maximum number of iterations
	 */
	unsigned int _maxIterations;

	/**
	 * final transformation matrix, found after iteration (fixed dimensions for 2D and 3D case)
	 */
	Matrix* _Tfinal4x4;

	/**
	 * transformation matrix of last step
	 */
	Matrix* _Tlast;

	/**
	 * dimension of space
	 */
	int _dim;
};

}

#endif /* REGISTRATION_H_ */
