#include "Ndt.h"

#include "obcore/base/tools.h"
#include "obcore/base/Timer.h"
#include "obcore/math/mathbase.h"

namespace obvious {

Ndt::Ndt(int minX, int maxX, int minY, int maxY, double cellSize) :
		Registration() {

	//Initialize Cell Structure
	_minX = minX;
	_maxX = maxX;
	_minY = minY;
	_maxY = maxY;
	_cellSize = cellSize;
	_numCellsX = floor(_maxX - _minX) / _cellSize;
	_numCellsY = floor(_maxY - _minY) / _cellSize;
	System<NdtCell>::allocate(_numCellsY, _numCellsX, _model);
	for (int y = 0; y < _numCellsY; y++) {
		for (int x = 0; x < _numCellsX; x++) {
			_model[y][x].centroid = new double[2];
			_model[y][x].cov = new Matrix(2, 2);
			_model[y][x].cov_inv = new Matrix(2, 2);
		}
	}

	//Initialize Members;
	Registration::_maxIterations = 3;
	Registration::_dim = 2;
	_scene = NULL;
	_sceneTmp = NULL;
	_searchTmp = NULL;
	_sizeSceneBuf = 0;
	_sizeScene = 0;

	Registration::_Tfinal4x4 = new Matrix(4, 4);
	Registration::_Tlast = new Matrix(4, 4);
	Registration::_Tfinal4x4->setIdentity();
	Registration::_Tlast->setIdentity();

	this->reset();
}

Ndt::~Ndt() {
	System<NdtCell>::deallocate(_model);
	if (_scene != NULL)
		System<double>::deallocate(_scene);
	if (_sceneTmp != NULL)
		System<double>::deallocate(_sceneTmp);
	delete Registration::_Tlast;
}

void Ndt::setModel(Matrix* coords, double probability) {

	if (_trace) {
		_trace->reset();
		_trace->setModel(coords);
	}

	if (coords->getCols() != (size_t) _dim) {
		cout << "WARNING: Model is not of correct dimensionality. Needed: "
				<< _dim << endl;
		return;
	}

	unsigned int sizeSource = coords->getRows();
	unsigned int sizeModel = sizeSource;
	bool* mask = Registration::createSubsamplingMask(&sizeModel, probability);

	for (unsigned int i = 0; i < sizeSource; i++) {
		if (mask[i]) //fixme take all points for now
		{
			double* coord = new double[2];
			coord[0] = (*coords)(i, 0);
			coord[1] = (*coords)(i, 1);
			int y = floor((coord[1] - _minY) / _cellSize);
			int x = floor((coord[0] - _minX) / _cellSize);
			_model[y][x].coords.push_back(coord);
		}
	}

	//Iterate cells and  calculate covariance and centroid if a cell is occupied Eq. 6.2 & 6.3 (Magnusson,2009)
	for (int y = 0; y < _numCellsY; y++) {
		for (int x = 0; x < _numCellsX; x++) {
			NdtCell cell = _model[y][x];
			if (!cell.isOccupied())
				continue;

			vector<double*> v = cell.coords;

			// Centroid
			double sumX = 0.0;
			double sumY = 0.0;
			for (unsigned int i = 0; i < v.size(); i++) {
				cout << x << " " << y << " " << v[i][0] << " " << v[i][1]
						<< endl;
				sumX += v[i][0];
				sumY += v[i][1];
			}
			cell.centroid[0] = sumX / (double) v.size();
			cell.centroid[1] = sumY / (double) v.size();
			cout << "centroid: " << cell.centroid[0] << " " << cell.centroid[1]
					<< endl << endl;
			;

			//Covariance
			Matrix* cov = cell.cov;
			(*cov)(0, 0) = 0.0;
			(*cov)(0, 1) = 0.0;
			(*cov)(1, 0) = 0.0;
			(*cov)(1, 1) = 0.0;
			for (unsigned int i = 0; i < v.size(); i++) {
				double c[2];
				c[0] = v[i][0] - cell.centroid[0];
				c[1] = v[i][1] - cell.centroid[1];
				(*cov)(0, 0) += c[0] * c[0];
				(*cov)(0, 1) += c[0] * c[1];
				(*cov)(1, 0) += c[1] * c[0];
				(*cov)(1, 1) += c[1] * c[1];
			}
			(*cov)(0, 0) /= v.size() - 1;
			(*cov)(0, 1) /= v.size() - 1;
			(*cov)(1, 0) /= v.size() - 1;
			(*cov)(1, 1) /= v.size() - 1;
			*(cell.cov_inv) = cov->getInverse();
			cov->print();
			cell.cov_inv->print();
		}
	}

	delete[] mask;
}

void Ndt::setScene(Matrix* coords, double probability) {

	if (coords->getCols() != (size_t) _dim) {
		cout << "WARNING: Scene is not of correct dimensionality " << _dim
				<< endl;
		return;
	}

	unsigned int sizeSource = coords->getRows();
	_sizeScene = sizeSource;
	bool* mask = Registration::createSubsamplingMask(&_sizeScene, probability);

	Registration::checkMemory(_sizeScene, _dim, _sizeSceneBuf, _scene);
	unsigned int idx = 0;
	for (unsigned int i = 0; i < sizeSource; i++) {
		if (mask[i]) {
			for (unsigned int j = 0; j < (unsigned int) _dim; j++)
				_scene[idx][j] = (*coords)(i, j);
			idx++;
		}
	}
	if (_sceneTmp)
		System<double>::deallocate(_sceneTmp);
	System<double>::allocate(_sizeScene, _dim, _sceneTmp);
	System<double>::copy(_sizeScene, _dim, _scene, _sceneTmp);

	delete[] mask;
}

void Ndt::reset() {
	Registration::_Tfinal4x4->setIdentity();
	if (_sceneTmp)
		System<double>::copy(_sizeScene, _dim, _scene, _sceneTmp);
}

NdtCell Ndt::getCorrespondingCell(double pointX, double pointY) {
	int x = floor((pointX - _minX) / _cellSize);
	int y = floor((pointY - _minY) / _cellSize);
	if (x < 0 || x > _numCellsY || y < 0 || y > _numCellsY) {
		assert(false);
	}

	/** get information from related cell**/
	NdtCell cell = _model[y][x];
	return _model[y][x];
}

//only used for lineSearch!!!!!!!!
double Ndt::calculateScore(double** scene, Eigen::Vector3d &poseIncrement) {

	Eigen::Matrix4d tf = Eigen::Matrix4d::Zero();
	composeTransformation(tf, poseIncrement(2), poseIncrement(0),
			poseIncrement(1));
	Matrix* obTF = new Matrix(4, 4);
	EigenMatrix4dToObviouslyMatrix(obTF, tf);
	Registration::applyTransformation(scene, _sizeScene, _dim, obTF);

	double score;
	double c1, c2, d3, d1, d2;
	double integral, outlier_ratio, support_size;
	integral = 0.1;
	outlier_ratio = 0.35;
	support_size = 1.0; //=current resolution
	c1 = (1 - outlier_ratio) / integral;
	c2 = outlier_ratio / pow(support_size, 3);
	d3 = -log(c2);
	d1 = -(-log(c1 + c2) - d3);
	d2 = -log((-log(c1 * exp(-0.5) + c2) - d3) / -d1);

	for (unsigned int j = 0; j < _sizeScene; j++) {
		/** project transformed scene point to cell **/
		NdtCell cell = getCorrespondingCell(scene[j][0], scene[j][1]);
		if (!cell.isOccupied())
			continue; //Jump is there are not enough points in the cell or no cell was found //fixme use neighbors

		//fixme known inconsistency because of obviously->eigen conversion
		Matrix *cov_inv = new Matrix(*cell.cov_inv);
		Eigen::Vector2d eMappedPoint = Eigen::Vector2d::Zero(); //relative point to cell centroid
		eMappedPoint(0) = scene[j][0] - cell.centroid[0];
		eMappedPoint(1) = scene[j][1] - cell.centroid[1];
		Eigen::Matrix2d eCovInv = Eigen::Matrix2d::Zero();
		eCovInv(0, 0) = (*cov_inv)(0, 0);
		eCovInv(0, 1) = (*cov_inv)(0, 1);
		eCovInv(1, 0) = (*cov_inv)(1, 0);
		eCovInv(1, 1) = (*cov_inv)(1, 1);

		/** Calcuate Point Score **/
		double l = eMappedPoint.transpose() * eCovInv * eMappedPoint; //likelyhood
		double px = -d1 * exp(-d2 * 0.5 * l);  //point score
		score -= px; //overall score

	}
	return score;
}

EnumState Ndt::step(Eigen::Matrix3d &hessian, Eigen::Vector3d &score_gradient,
		double &score) {

	double** transformedScene;
	System<double>::allocate(_sizeScene, _dim, transformedScene);
	transformedScene = _scene;
	Registration::applyTransformation(transformedScene, _sizeScene, _dim,
			Registration::_Tfinal4x4);

	//NOTE: SceneTMP already contains the transformed points according to the latest transformation parameters

	//Parameterization of score function; Eq. 6.8 (Magnusson,2009)
	double c1, c2, d3, d1, d2;
	double integral, outlier_ratio, support_size;
	integral = 0.1;
	outlier_ratio = 0.35;
	support_size = 1.0; //=current resolution
	c1 = (1 - outlier_ratio) / integral;
	c2 = outlier_ratio / pow(support_size, 3);
	d3 = -log(c2);
	d1 = -(-log(c1 + c2) - d3);
	d2 = -log((-log(c1 * exp(-0.5) + c2) - d3) / -d1);

	EnumState retState = NOTMATCHABLE;
	for (unsigned int j = 0; j < _sizeScene; j++) {
		/** project transformed scene point to cell **/
		int x = floor((transformedScene[j][0] - _minX) / _cellSize);
		int y = floor((transformedScene[j][1] - _minY) / _cellSize);
		if (x < 0 || x > _numCellsY || y < 0 || y > _numCellsY) {
			assert(false);
			continue;
		}

		/** get information from related cell**/
		NdtCell cell = _model[y][x];
		if (!cell.isOccupied())
			continue; //Jump is there are not enough points in the cell //fixme use neighbors
		retState = PROCESSING; //necessary to known if there was even one iteration through the code

		//fixme known inconsistency because of obviously->eigen conversion
		Matrix *cov_inv = new Matrix(*cell.cov_inv);
		Eigen::Vector2d eTransformedPt(transformedScene[j][0],
				transformedScene[j][1]);
		Eigen::Vector2d eMappedPoint = Eigen::Vector2d::Zero(); //relative point to cell centroid
		eMappedPoint(0) = eTransformedPt(0) - cell.centroid[0];
		eMappedPoint(1) = eTransformedPt(1) - cell.centroid[1];
		Eigen::Matrix2d eCovInv = Eigen::Matrix2d::Zero();
		eCovInv(0, 0) = (*cov_inv)(0, 0);
		eCovInv(0, 1) = (*cov_inv)(0, 1);
		eCovInv(1, 0) = (*cov_inv)(1, 0);
		eCovInv(1, 1) = (*cov_inv)(1, 1);

		/** Calcuate Point Score Eq. 6.9 & 6.10 (Magnusson,2009)**/
		double l = eMappedPoint.transpose() * eCovInv * eMappedPoint; //likelyhood
		double px = -d1 * exp(-d2 * 0.5 * l); //point score
		score -= px; //overall score
		//cout<<"likelihood "<<l<<"\n point score "<< px<<endl;

		/** Jacobian 2D NDT Eq. 6.15 (Magnusson,2009)**/
		double angleZ = atan2((*Registration::_Tfinal4x4)(1, 0),
				(*Registration::_Tfinal4x4)(0, 0));
		Eigen::MatrixXd eJacobian = Eigen::MatrixXd::Zero(2, 3);
		eJacobian(0, 0) = 1;		//first column
		eJacobian(1, 0) = 0;
		eJacobian(0, 1) = 0;		//second column
		eJacobian(1, 1) = 1;
		eJacobian(0, 2) = -eTransformedPt(0) * sin(angleZ)
				- eTransformedPt(1) * cos(angleZ);
		eJacobian(1, 2) = eTransformedPt(0) * cos(angleZ)
				- eTransformedPt(1) * sin(angleZ);

		/** Score Gradient Eq. 6.12 (Magnusson,2009)**/
		//cout<<"point "<<eMappedPoint<<"\n Cov"<<eCovInv.matrix() <<"\n jacobi "<<eJacobian.matrix()<<endl;
		double factor = exp(
				-d2 * 0.5
						* (double) (eMappedPoint.transpose() * eCovInv
								* eMappedPoint));
		for (int g = 0; g < 3; g++) {
			double xCJg = eMappedPoint.transpose() * eCovInv * eJacobian.col(g);
			score_gradient[g] += d1 * d2 * xCJg * factor;
		}

		/** Hessian Caculation Eq. 6.13 (Magnusson,2009) **/
		for (int g = 0; g < 3; g++) {
			for (int h = 0; h < 3; h++) {
				//second order partial derivatives Eq. 6.8 (Magnusson,2009)
				// use different values if g & h == 2
				Eigen::Vector2d secondDerivativeVec(0, 0);
				if (g == 2 && h == 2) {
					secondDerivativeVec(0) = -eTransformedPt(0) * cos(angleZ)
							+ eTransformedPt(0) * sin(angleZ);
					secondDerivativeVec(1) = -eTransformedPt(0) * sin(angleZ)
							- eTransformedPt(0) * cos(angleZ);
				}
				double xCJg = eMappedPoint.transpose() * eCovInv
						* eJacobian.col(g);
				double xCJh = eMappedPoint.transpose() * eCovInv
						* eJacobian.col(h);
				double xCSec = eMappedPoint.transpose() * eCovInv
						* secondDerivativeVec;
				double JhCJg = eJacobian.col(h).transpose() * eCovInv
						* eJacobian.col(g);
				hessian(g, h) += d1 * d2 * factor * (-d2 * xCJg * xCJh) + xCSec
						+ JhCJg;
			}

		}
	}
	return retState;
}

EnumState Ndt::iterate(double* rms, unsigned int* iterations, Matrix* Tinit) {

	unsigned int iter = 0;

	if (Tinit) {
		applyTransformation(_sceneTmp, _sizeScene, _dim, Tinit);
		(*Registration::_Tfinal4x4) = (*Tinit) * (*Registration::_Tfinal4x4);
	}

	EnumState eRetval = PROCESSING;
	while (eRetval == PROCESSING) {
		cout << "Iteration: " << iter << endl;
		/** Reset iteration parameters **/
		double score = 0;
		Eigen::Vector3d score_gradient = Eigen::Vector3d::Zero();
		Eigen::Matrix3d hessian = Eigen::Matrix3d::Zero();

		//calculate one iteration and check return for errors
		eRetval = step(hessian, score_gradient, score);
		if (eRetval == NOTMATCHABLE) {
			cout << "ERROR: There were no points in occupied cells" << endl;
			break;
		} else if (score == 0) {
			cout << "ERROR:Score is zero" << endl;
			eRetval = ERROR;
			break;
		}

		/**Solve registration equation Hessian * deltaParam = - gradient **/
		Eigen::Vector3d deltaParam;
		//deltaParam = hessian.inverse() * -score_gradient;
		deltaParam = hessian.fullPivLu().solve(-score_gradient);
		deltaParam.normalize();	//fixeme determine a stepsize here e.g. with More & Thuente line search
// 		//Line Search
//		double lineSearchResult = lineSearch(score_gradient, deltaParam, score);
//		cout<<"Line Search"<<lineSearchResult<<endl;
//		deltaParam *= lineSearchResult;

		deltaParam *= 0.1; //todo use line search here e.g. More & Thuente, 1992
		cout << "Iteration finished with:\nScore: " << score << "\n Hessian: "
				<< hessian.matrix() << "\n gradient " << score_gradient
				<< "\n Delta Param: " << deltaParam << "\n\n\n" << endl;

		/** Apply Transformation **/
		Eigen::Matrix4d deltaTF = Eigen::Matrix4d::Zero();
		composeTransformation(deltaTF, deltaParam(2), deltaParam(0),
				deltaParam(1));
		Matrix* obDeltaTF = new Matrix(4, 4);
		EigenMatrix4dToObviouslyMatrix(obDeltaTF, deltaTF);
		Registration::_Tlast = obDeltaTF;
		Registration::applyTransformation(_sceneTmp, _sizeScene, _dim,
				Registration::_Tlast);

		(*_Tfinal4x4) = (*_Tlast) * (*_Tfinal4x4);
		iter++;

		if (iter >= Registration::_maxIterations)
			eRetval = MAXITERATIONS;

		if (_trace) {
			_trace->addAssignment(_sceneTmp, _sizeScene);
		}
	}
	*iterations = iter;
	return eRetval;
}

void Ndt::EigenMatrix4dToObviouslyMatrix(obvious::Matrix* obviousMat,
		Eigen::Matrix4d eigenMat) {
//Rotation
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			(*obviousMat)(i, j) = eigenMat(i, j);
		}
	}
}

void Ndt::composeTransformation(Eigen::Matrix4d& tf, double angle, double t_x,
		double t_y) {
	Eigen::AngleAxisd appliedRotation;
	Eigen::Translation3d appliedTranslation;

	appliedRotation = Eigen::AngleAxisd(angle, Eigen::Vector3d::UnitZ());
	appliedTranslation = Eigen::Translation3d(t_x, t_y, 0);
	tf = Eigen::Matrix4d((appliedTranslation * appliedRotation).matrix());

}

void Ndt::activateTrace() {
	_trace = new NdtTrace(Registration::_dim);
}

void Ndt::deactivateTrace() {
	delete _trace;
	_trace = NULL;
}

void Ndt::serializeTrace(char* folder, unsigned int delay) {
	if (_trace)
		_trace->serialize(folder, delay);
}

}

/////////////////////////////////////////////////////////////////////////////////////////////////
// LINE SEARCH Backtracking
////////////////////////////////////////////////////////////////////////////////////////////////
double Ndt::lineSearch(Eigen::Vector3d &gradientInit,
		Eigen::Vector3d &poseIncrement, double scoreInit) {

	double stepSizeMax = 0.5;
	double stepSize = stepSizeMax;
	double shrinkF = 0.5; //T
	double controlF = 0.5; //c
	double m = poseIncrement.transpose() * gradientInit;

	double t = -shrinkF * m;

	//from magnusson code
//	if(m >= 0.0) {
//		poseIncrement = -poseIncrement;
//		m = -m;
//		if (m >= 0.0){return 0.05;} //recovery step
//	}

	if (_searchTmp)
		System<double>::deallocate(_searchTmp);
	System<double>::allocate(_sizeScene, _dim, _searchTmp);

	int iter = 0;
	double scoreInc = 0;
	while ((scoreInit - scoreInc) < (stepSize * t)) { //Armijo–Goldstein condition
		stepSize *= controlF;

		Eigen::Vector3d reducedIncrement = poseIncrement * stepSize;
		//calculate score for f(x + steSize * poseInc)
		System<double>::copy(_sizeScene, _dim, _sceneTmp, _searchTmp);
		scoreInc = calculateScore(_searchTmp, reducedIncrement);

		double test1 = scoreInit - scoreInc;
		double test2 = stepSize * t;

		if (test2 == 0) {
			cout << "Recovery" << endl; //recovery step
			return 0.05;
		}

		//cout<<"StepSizeIT "<<stepSize<<"\n t="<<t<< "\nfirst part "<<test1<<"\n second "<< test2 <<"\nScoreInc "<< scoreInc<<endl;
		iter++;
	}

	//cout<<"Stepsize: "<<stepSize<<endl;
	return stepSize;
}

