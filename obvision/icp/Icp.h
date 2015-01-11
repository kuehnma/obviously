#ifndef ICP_H_
#define ICP_H_

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
using namespace std;

#include "obvision/registration/Registration.h"
#include "obvision/icp/assign/PairAssignment.h"
#include "obvision/icp/assign/filter/IPostAssignmentFilter.h"
#include "obvision/icp/IRigidEstimator.h"
#include "obcore/base/CartesianCloud.h"
#include "obcore/base/System.h"

#include "obcore/math/linalg/linalg.h"

#include "obvision/icp/IcpTrace.h"

using namespace obvious;

namespace obvious {

/**
 * @class Icp
 * @brief Represents the iterative closest point algorithm
 * @author Stefan May and Dirk Holz
 **/
class Icp: public obvious::Registration {
public:
	/**
	 * Standard constructor
	 * @param assigner pair assignment strategy for point clouds
	 * @param estimator transformation estimator
	 * @param pDistanceThreshold thresholding strategy
	 */
	Icp(PairAssignment* assigner, IRigidEstimator* estimator);

	/**
	 * Destructor
	 */
	~Icp();

	/**
	 * Activate internal trace writing. While the method iterate is executed, all states are traced.
	 */
	void activateTrace();

	/**
	 * Deactivate internal trace writing.
	 */
	void deactivateTrace();

	/**
	 * Accessor to pair assigner
	 * @return assigner
	 */
	PairAssignment* getPairAssigner();

	/**
	 * Accessor to rigid estimator
	 * @return estimator
	 */
	IRigidEstimator* getRigidEstimator();

	/**
	 * Copy model to internal buffer
	 * @param coords model coordinates, as tuples or triples
	 * @param normals model normals, as tuples or triples, may be NULL
	 * @param size number of points, i.e. coordinate triples
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
	void setModel(double* coords, double* normals, const unsigned int size,
			double probability = 1.0);

	/**
	 * Convenience method extracting all points from cloud to double array
	 * @param coords model coordinates
	 * @param normals model normals, may be NULL
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
	void setModel(Matrix* coords, Matrix* normals = NULL, double probability =
			1.0);

	/**
	 * Copy scene to internal buffer
	 * @param coords scene coordinates, as tuples or triples
	 * @param normals scene normals, as tuples or triples, may be NULL
	 * @param size number of points, i.e. coordinate triples
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
	void setScene(double* coords, double* normals, const unsigned int size,
			double probability = 1.0);

	/**
	 * Convenience method extracting data from cloud to double array
	 * @param coords scene coordinates
	 * @param normals scene normals, may be NULL
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
	void setScene(Matrix* coords, Matrix* normals = NULL, double probability =
			1.0);

	/**
	 * Reset state of ICP algorithm (resets iteration count, estimated transformation and point pairs)
	 */
	void reset();

	/**
	 * @deprecated
	 * Set maximal RMS error interrupting iteration
	 * @param rms RMS error
	 */
	void setMaxRMS(double rms);

	/**
	 * @deprecated
	 * Get maximal RMS error interrupting iteration
	 * @return RMS error
	 */
	double getMaxRMS();

	/**
	 * @deprecated
	 * If the RMS error is not reducing within this number, the matching process is considered as successful.
	 * @param convCnt convergence counter
	 */
	void setConvergenceCounter(unsigned int convCnt);

	/**
	 * @deprecated
	 * Access convergence counter
	 * @param convCnt convergence counter
	 */
	unsigned int getConvergenceCounter();




	/**
	 * Perform one iteration
	 * @param rms return value of RMS error
	 * @param pairs return value of pair assignments, i.e. number of pairs
	 * @return processing state
	 */
	EnumState step(double* rms, unsigned int* pairs);

	/**
	 * Start iteration
	 * @param rms return value of RMS error
	 * @param pairs return value of pair assignments, i.e. number of pairs
	 * @param iterations return value of performed iterations
	 * @param Tinit apply initial transformation before iteration
	 * @return  processing state
	 */
	EnumState iterate(double* rms, unsigned int* pairs,
			unsigned int* iterations, Matrix* Tinit = NULL);

	/**
	 * Start registration process
	 * @param Tinit Apply initial transformation before iteration.
	 * @param verbose Set true if information about the registration are required in the console.
	 * @return  processing state
	 */
	EnumState align(Matrix* Tinit = NULL, bool verbose = false);

	/**
	 * Abstract method implemented by derived classes for loading algorithm specific parameters
	 * @param file Path to file.
	 */
	int loadParametersFromXML(string filepath);

	/**
	 * Serialize assignment to trace folder
	 * @param folder trace folder (must not be existent)
	 * @param delay animation delay (specified in delay*1/100s)
	 */
	void serializeTrace(char* folder, unsigned int delay = 10);



private:

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
	 * maximal RMS error interrupting iteration
	 */
	double _maxRMS;

	/**
	 * the model (coordinates)
	 */
	double** _model;

	/**
	 * normals of model
	 */
	double** _normalsM;

	/**
	 * normals of scene
	 */
	double** _normalsS;

	/**
	 * local copy of scene normals
	 */
	double** _normalsSTmp;

	/**
	 * size of model
	 */
	unsigned int _sizeModel;

	/**
	 * size of internal model buffer
	 */
	unsigned int _sizeModelBuf;

	/**
	 * size of internal scene buffer
	 */
	unsigned int _sizeSceneBuf;

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

	/**
	 * point assigner
	 */
	PairAssignment* _assigner;

	/**
	 * transformation estimator
	 */
	IRigidEstimator* _estimator;

	/**
	 * convergence counter
	 */
	unsigned int _convCnt;

	/**
	 * tracing instance (applied while iterating)
	 */
	IcpTrace* _trace;
};

}

#endif /*ICP_H_*/
