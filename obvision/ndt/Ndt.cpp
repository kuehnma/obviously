#include "Ndt.h"

#include "obcore/base/tools.h"
#include "obcore/base/Timer.h"
#include "obcore/math/mathbase.h"

namespace obvious {

const char* g_ndt_states[] = { "NDT_IDLE", "NDT_PROCESSING", "NDT_NOTMATCHABLE",
		"NDT_MAXITERATIONS", "NDT_TIMEELAPSED", "NDT_SUCCESS", "NDT_ERROR" };

Ndt::Ndt(int minX, int maxX, int minY, int maxY) {

	_minX = minX;
	_maxX = maxX;
	_minY = minY;
	_maxY = maxY;

	System<NdtCell>::allocate(_maxY - _minY, _maxX - _minX, _model);
	for (int y = 0; y < _maxY - _minY; y++) {
		for (int x = 0; x < _maxX - _minX; x++) {
			_model[y][x].centroid = new double[2];
			_model[y][x].cov = new Matrix(2, 2);
			_model[y][x].cov_inv = new Matrix(2, 2);
		}
	}

	Registration::_maxIterations = 3;
	Registration::_dim = 2;
	_scene = NULL;
	_sceneTmp = NULL;
	_sizeSceneBuf = 0;
	_sizeScene = 0;

	Registration::_Tfinal4x4 = new Matrix(4, 4);
	Registration::_Tlast = new Matrix(4, 4);
	Registration::_Tfinal4x4->setIdentity();
	Registration::_Tlast->setIdentity();

	// magic numbers taken from ROS implementation
	_d1 = 1.0;
	_d2 = 0.05;

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
	if (coords->getCols() != (size_t) _dim) {
		cout << "WARNING: Model is not of correct dimensionality. Needed: "
				<< _dim << endl;
		return;
	}

	unsigned int sizeSource = coords->getRows();
	unsigned int sizeModel = sizeSource;
	bool* mask = Registration::createSubsamplingMask(&sizeModel, probability);

	for (unsigned int i = 0; i < sizeSource; i++) {
		//if(mask[i])
		if (true) //fixme take all points
		{
			double* coord = new double[2];
			coord[0] = (*coords)(i, 0);
			coord[1] = (*coords)(i, 1);
			int x = floor(coord[0]) - _minX;
			int y = floor(coord[1]) - _minY;
			_model[y][x].coords.push_back(coord);
		}
	}

	//Iterate cells and  calculate covariance and centroid if a cell is occupied
	for (int y = 0; y < _maxY - _minY; y++) {
		for (int x = 0; x < _maxX - _minX; x++) {
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
			(*cov)(0, 0) /= v.size();
			(*cov)(0, 1) /= v.size();
			(*cov)(1, 0) /= v.size();
			(*cov)(1, 1) /= v.size();
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

EnumState Ndt::step(Matrix &hessian, double &score_gradient[3], double &score) {

	for (unsigned int j = 0; j < _sizeScene; j++) {

		double transformedScene[_sizeScene][_dim];
		applyTransformation(transformedScene, _sizeScene, _dim, Tfinal4x4);

		// project transformed scene point to cell
		int x = floor(transformedScene[j][0] - _minX);
		int y = floor(transformedScene[j][1] - _minY);

		//Check if projection was successful
		if (x < 0 && x > _maxX && y < 0 && y > _maxY) {
			continue;
		}

		NdtCell cell = _model[y][x];
		// coord zero mean
		Vector c_zm(2);
		c_zm(0) = transformedScene[j][0] - cell.centroid[0];
		c_zm(1) = transformedScene[j][1] - cell.centroid[1];
		Matrix *cov_inv = cell.cov_inv;

		Vector tmp = Matrix::multiply(*cov_inv, c_zm, false);

		// likelihood
		double l = c_zm(0) * tmp(0) + c_zm(1) * tmp(1);

		//score
		double px = exp(-0.5 * l);
		score -= px;

		Matrix jacobian1 = new Matrix(2, 1);
		Matrix jacobian2 = new Matrix(2, 1);
		Matrix jacobian3 = new Matrix(2, 1);
		jacobian1(0, 0) = 1;
		jacobian1(1, 0) = 0;
		jacobian2(0, 0) = 0;
		jacobian2(1, 0) = 1;
		jacobian3(0, 0) = -_sceneTmp[j][0] * sin(angle)
				- _sceneTmp[j][1] * cos(angle);
		jacobian3(1, 2) = _sceneTmp[j][0] * cos(angle)
				- _sceneTmp[j][1] * sin(angle);

		//gradient
		Vector grad(3);
		double factor = exp(-0.5 * (c_zm(0) * tmp(0) + c_zm(1) * tmp(1)));
		double xInvC1 = c_zm(0) * *cov_inv(0, 0) + c_zm(1) * *cov_inv(1, 0);
		double xInvC2 = c_zm(0) * *cov_inv(0, 1) + c_zm(1) * *cov_inv(1, 1);
		double xInvCJ3 = (xInvC1 * j30 + xInvC2 * j31);

		double j30 = -_sceneTmp[j][0] * sin(angle)
				- _sceneTmp[j][1] * cos(angle);
		double j31 = _sceneTmp[j][0] * cos(angle)
				- _sceneTmp[j][1] * sin(angle);
		score_gradient[0] += xInvC1 * factor;
		score_gradient[1] += xInvC2 * factor;
		score_gradient[2] += xInvCJ3 * factor;

		//hessian
		double j3InvCj3 =
		hessian(0,0) += factor * (xInvC1 * -xInvC1  + 0 + *cov_inv(0,0) );
		hessian(0,1) += factor * (xInvC1 * -xInvC2  + 0 + *cov_inv(1,0) );
		hessian(0,2) += factor * (xInvC1 * -xInvCJ3 + 0 + *cov_inv() )
	}
}

EnumState Ndt::iterate(double* rms, unsigned int* iterations, Matrix* Tinit) {
	Registration::_Tfinal4x4->setIdentity();

	if (Tinit) {
		applyTransformation(_sceneTmp, _sizeScene, _dim, Tinit);
		(*Registration::_Tfinal4x4) = (*Tinit) * (*Registration::_Tfinal4x4);
	}

	EnumState eRetval = PROCESSING;
	unsigned int iter = 0;
	//double rms_prev = 10e12;

	for (unsigned int i = 0; i < Registration::_maxIterations; i++) {
		double score = 0;
		double score_gradient[3] = { 0, 0, 0 };
		Matrix hessian = new Matrix(3, 3);

		eRetval = step(hessian, score_gradient, score)

		Registration::applyTransformation(_sceneTmp, _sizeScene, _dim,
				Registration::_Tlast);
		iter++;
	}

	*iterations = iter;

	return eRetval;
}

}
