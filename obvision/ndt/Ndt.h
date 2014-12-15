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
#include <eigen3/Eigen/Dense>
#include "obvision/ndt/NdtTrace.h"
#include "obcore/base/CartesianCloud.h"

using namespace obvious;

namespace obvious {

struct NdtCell {
	vector<double*> coords;
	double* centroid;
	Matrix* cov;
	Matrix* cov_inv;
	bool isOccupied() {
		return (coords.size() >= 3);
	}
	;
};

/**
 * @class Ndt
 * @brief Represents the normal distribution transform
 * @author Stefan May
 **/
class Ndt: public obvious::Registration {
public:
	/**
	 * Standard constructor
	 */
	Ndt(int minX, int maxX, int minY, int maxY, double cellSize);

	/**
	 * Destructor
	 */
	~Ndt();

	/**
	 * Sample model point cloud to NDT space
	 * @param coords model coordinates
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
	void setModel(Matrix* coords, double probability = 1.0);

	/**
	 * Copy scene to internal buffer
	 * @param coords scene coordinates
	 * @param probability probability of coordinates of being sampled (range [0.0 1.0])
	 */
	void setScene(Matrix* coords, double probability = 1.0);

	/**
	 * Reset state of NDT algorithm
	 */
	void reset();

	/**
	 * Encapsulates one step of calculating the score, hessian matrix and score gradient for each iteration
	 * @param hessian Hessian matrix
	 * @param score_gradient Gradient of the score
	 * @param score Score for this iteration(Negativ is good)
	 * @return processing state
	 */
	EnumState step(Eigen::Matrix3d &hessian, Eigen::Vector3d &score_gradient,
			double &score);

	/**
	 * Start iteration
	 * @param rms return value of RMS error
	 * @param iterations return value of performed iterations
	 * @param Tinit apply initial transformation before iteration
	 * @return  processing state
	 */
	EnumState iterate(double* rms, unsigned int* iterations, Matrix* Tinit =
	NULL);

	/**
	 * Activate internal trace writing. While the method iterate is executed, all states are traced.
	 */
	void activateTrace();

	/**
	 * Deactivate internal trace writing.
	 */
	void deactivateTrace();

	/**
	 * Serialize assignment to trace folder
	 * @param folder trace folder (must not be existent)
	 * @param delay animation delay (specified in delay*1/100s)
	 */
	void serializeTrace(char* folder, unsigned int delay = 10);

private:

	/**
	 * Compose a 2D Transformation in a 4x4 Matrix.
	 * The Matrix is composed by the Rotation around the Z - and X,Y translation
	 * @param tf Empty Eigen Matrix 4d as an in/out parameter
	 * @param angle The rotation around the Z-Axis
	 * @param t_x The translation in x-direction
	 * @param t_y The translation in y-direction
	 */
	static void composeTransformation(Eigen::Matrix4d& tf, double angle,
			double t_x, double t_y);

	/**
	 * Convert between different Matrix Classes
	 */
	static void EigenMatrix4dToObviouslyMatrix(obvious::Matrix* obviousMat,
			Eigen::Matrix4d eigenMat);
	/**
	 * caculate score for a given parametervektor offset and scene
	 * Used for line search
	 */

	/**
	 * Finds a Corresponding cell for coordinats x,y
	 * @param pointX X coordinate of a point
	 * @param pointY Y coordinate of a point
	 * @return A corresponding cell or NULL for nothing
	 */
	NdtCell getCorrespondingCell(double pointX, double pointY);


	double calculateScore(double** scene, Eigen::Vector3d &poseIncrement);

	/**
	 * Search for a good stepsize with backtracking
	 */
	double lineSearch(Eigen::Vector3d &gradientInit, Eigen::Vector3d &poseIncrement, double scoreInit);

	int _minX;
	int _maxX;
	int _minY;
	int _maxY;

	/**
	 * Size of a squared cell. (cellSize * cellSize)
	 */
	double _cellSize;

	/**
	 * size of internal scene buffer
	 */
	unsigned int _sizeSceneBuf;

	NdtCell** _model;

	/**
	 * Number of columns in the model
	 */
	int _numCellsX;

	/**
	 * Number of rows in the model
	 */
	int _numCellsY;

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
	 * Trace instance for debugging
	 */
	NdtTrace* _trace;

	double** _searchTmp;


};

}

#endif /*NDT_H_*/
