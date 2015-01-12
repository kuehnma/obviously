#include "Icp.h"
#include "obcore/base/tools.h"
#include "obcore/base/Timer.h"
#include "obcore/math/mathbase.h"

namespace obvious {

Icp::Icp(PairAssignment* assigner, IRigidEstimator* estimator) :
		Registration() {
	_assigner = assigner;
	_estimator = estimator;

	Registration::_maxIterations = 3;
	Registration::_dim = _assigner->getDimension();
	_maxRMS = 0.1;
	_model = NULL;
	_scene = NULL;
	_sceneTmp = NULL;
	_normalsM = NULL;
	_normalsS = NULL;
	_normalsSTmp = NULL;
	_sizeModelBuf = 0;
	_sizeSceneBuf = 0;
	_sizeModel = 0;
	_sizeScene = 0;

	Registration::_Tfinal4x4 = new Matrix(4, 4);
	Registration::_Tlast = new Matrix(4, 4);
	Registration::_Tfinal4x4->setIdentity();
	Registration::_Tlast->setIdentity();
	_convCnt = 5;

	this->reset();

	_trace = NULL;
}

Icp::~Icp() {
	if (_model != NULL)
		System<double>::deallocate(_model);
	if (_scene != NULL)
		System<double>::deallocate(_scene);
	if (_sceneTmp != NULL)
		System<double>::deallocate(_sceneTmp);
	if (_normalsM != NULL)
		System<double>::deallocate(_normalsM);
	if (_normalsS != NULL)
		System<double>::deallocate(_normalsS);
	if (_normalsSTmp != NULL)
		System<double>::deallocate(_normalsSTmp);
	delete Registration::_Tlast;
	if (_trace) {
		delete _trace;
		_trace = NULL;
	}
}

void Icp::activateTrace() {
	_trace = new IcpTrace(Registration::_dim);
}

void Icp::deactivateTrace() {
	delete _trace;
	_trace = NULL;
}

PairAssignment* Icp::getPairAssigner() {
	return _assigner;
}

IRigidEstimator* Icp::getRigidEstimator() {
	return _estimator;
}

void Icp::setModel(double* coords, double* normals, const unsigned int size,
		double probability) {
	_sizeModel = size;
	bool* mask = Registration::createSubsamplingMask(&_sizeModel, probability);

	unsigned int sizeNormalsBuf = _sizeModelBuf;
	Registration::checkMemory(_sizeModel, Registration::_dim, _sizeModelBuf,
			_model);
	unsigned int idx = 0;
	for (unsigned int i = 0; i < size; i++) {
		if (mask[i]) {
			for (unsigned int j = 0; j < (unsigned int) Registration::_dim; j++)
				_model[idx][j] = coords[Registration::_dim * i + j];
			idx++;
		}
	}

	if (normals) {
		Registration::checkMemory(_sizeModel, Registration::_dim,
				sizeNormalsBuf, _normalsM);
		idx = 0;
		for (unsigned int i = 0; i < size; i++) {
			if (mask[i]) {
				for (unsigned int j = 0; j < (unsigned int) Registration::_dim;
						j++)
					_normalsM[idx][j] = normals[Registration::_dim * i + j];
				idx++;
			}
		}
	}

	_assigner->setModel(_model, idx);
	_estimator->setModel(_model, idx, _normalsM);

	delete[] mask;
}

void Icp::setModel(Matrix* coords, Matrix* normals, double probability) {
	if (coords->getCols() != (size_t) Registration::_dim) {
		cout << "WARNING: Model is not of correct dimensionality. Needed: "
				<< Registration::_dim << endl;
		return;
	}

	unsigned int sizeSource = coords->getRows();
	_sizeModel = sizeSource;
	bool* mask = Registration::createSubsamplingMask(&_sizeModel, probability);

	unsigned int sizeNormals = _sizeModelBuf;

	Registration::checkMemory(_sizeModel, Registration::_dim, _sizeModelBuf,
			_model);

	unsigned int idx = 0;
	for (unsigned int i = 0; i < sizeSource; i++) {
		if (mask[i]) {
			for (unsigned int j = 0; j < (unsigned int) Registration::_dim; j++)
				_model[idx][j] = (*coords)(i, j);
			idx++;
		}
	}

	if (normals) {
		Registration::checkMemory(_sizeModel, Registration::_dim, sizeNormals,
				_normalsM);
		idx = 0;
		for (unsigned int i = 0; i < sizeSource; i++) {
			if (mask[i]) {
				for (unsigned int j = 0; j < (unsigned int) Registration::_dim;
						j++)
					_normalsM[idx][j] = (*normals)(i, j);
				idx++;
			}
		}
	}

	_assigner->setModel(_model, idx);
	_estimator->setModel(_model, idx, _normalsM);

	delete[] mask;
}

void Icp::setScene(double* coords, double* normals, const unsigned int size,
		double probability) {
	if (size == 0) {
		cout << "Scene of size 0 passed ... ignoring" << endl;
		return;
	}

	_sizeScene = size;
	bool* mask = Registration::createSubsamplingMask(&_sizeScene, probability);

	unsigned int sizeNormalsBuf = _sizeSceneBuf;
	Registration::checkMemory(_sizeScene, Registration::_dim, _sizeSceneBuf,
			_scene);
	unsigned int idx = 0;
	for (unsigned int i = 0; i < size; i++) {
		if (mask[i]) {
			for (unsigned int j = 0; j < (unsigned int) Registration::_dim; j++)
				_scene[idx][j] = coords[Registration::_dim * i + j];
			idx++;
		}
	}
	if (_sceneTmp)
		System<double>::deallocate(_sceneTmp);
	System<double>::allocate(_sizeScene, Registration::_dim, _sceneTmp);
	System<double>::copy(_sizeScene, Registration::_dim, _scene, _sceneTmp);

	if (normals) {
		Registration::checkMemory(_sizeScene, Registration::_dim,
				sizeNormalsBuf, _normalsS);
		idx = 0;
		for (unsigned int i = 0; i < size; i++) {
			if (mask[i]) {
				for (unsigned int j = 0; j < (unsigned int) Registration::_dim;
						j++)
					_normalsS[idx][j] = normals[Registration::_dim * i + j];
				idx++;
			}
		}
		if (_normalsSTmp)
			System<double>::deallocate(_normalsSTmp);
		System<double>::allocate(_sizeScene, Registration::_dim, _normalsSTmp);
		System<double>::copy(_sizeScene, Registration::_dim, _normalsS,
				_normalsSTmp);
	}

	delete[] mask;
}

void Icp::setScene(Matrix* coords, Matrix* normals, double probability) {
	if (coords->getCols() != (size_t) Registration::_dim) {
		cout << "WARNING: Scene is not of correct dimensionality "
				<< Registration::_dim << endl;
		return;
	}

	unsigned int sizeSource = coords->getRows();
	_sizeScene = sizeSource;
	bool* mask = Registration::createSubsamplingMask(&_sizeScene, probability);

	unsigned int sizeNormals = _sizeSceneBuf;

	Registration::checkMemory(_sizeScene, Registration::_dim, _sizeSceneBuf,
			_scene);
	unsigned int idx = 0;
	for (unsigned int i = 0; i < sizeSource; i++) {
		if (mask[i]) {
			for (unsigned int j = 0; j < (unsigned int) Registration::_dim; j++)
				_scene[idx][j] = (*coords)(i, j);
			idx++;
		}
	}
	if (_sceneTmp)
		System<double>::deallocate(_sceneTmp);
	System<double>::allocate(_sizeScene, Registration::_dim, _sceneTmp);
	System<double>::copy(_sizeScene, Registration::_dim, _scene, _sceneTmp);

	if (normals) {
		Registration::checkMemory(_sizeScene, Registration::_dim, sizeNormals,
				_normalsS);
		idx = 0;
		for (unsigned int i = 0; i < sizeSource; i++) {
			if (mask[i]) {
				for (unsigned int j = 0; j < (unsigned int) Registration::_dim;
						j++)
					_normalsS[idx][j] = (*normals)(i, j);
				idx++;
			}
		}
		if (_normalsSTmp)
			System<double>::deallocate(_normalsSTmp);
		System<double>::allocate(_sizeScene, Registration::_dim, _normalsSTmp);
		System<double>::copy(_sizeScene, Registration::_dim, _normalsS,
				_normalsSTmp);
	}

	delete[] mask;
}

void Icp::reset() {
	Registration::_Tfinal4x4->setIdentity();
	_assigner->reset();
	if (_sceneTmp)
		System<double>::copy(_sizeScene, Registration::_dim, _scene, _sceneTmp);
	if (_normalsSTmp)
		System<double>::copy(_sizeScene, Registration::_dim, _normalsS,
				_normalsSTmp);
}
/**
 * @deprecated
 */
void Icp::setMaxRMS(double rms) {
	_maxRMS = rms;
}

/**
 * @deprecated
 */
double Icp::getMaxRMS() {
	return _maxRMS;
}

/**
 * @deprecated
 */
void Icp::setConvergenceCounter(unsigned int convCnt) {
	_convCnt = convCnt;
}

/**
 * @deprecated
 */
unsigned int Icp::getConvergenceCounter() {
	return _convCnt;
}

EnumState Icp::step(double* rms, unsigned int* pairs) {
	Timer t;
	if (_model == NULL || _sceneTmp == NULL)
		return ERROR;

	EnumState retval = PROCESSING;

	vector<StrCartesianIndexPair>* pvPairs;
	_estimator->setScene(_sceneTmp, _sizeScene, _normalsSTmp);
	_assigner->determinePairs(_sceneTmp, _sizeScene);
	pvPairs = _assigner->getPairs();
	*pairs = pvPairs->size();

	if (_trace) {
		_trace->addAssignment(_sceneTmp, _sizeScene, *pvPairs);
	}

	if (pvPairs->size() > 2) {
		// Estimate transformation
		_estimator->setPairs(pvPairs);

		// get mapping error
		*rms = _estimator->getRMS();

		// estimate transformation
		_estimator->estimateTransformation(Registration::_Tlast);

		Registration::applyTransformation(_sceneTmp, _sizeScene,
				Registration::_dim, Registration::_Tlast);
		if (_normalsS)
			Registration::applyTransformation(_normalsSTmp, _sizeScene,
					Registration::_dim, Registration::_Tlast);

		// update overall transformation
		(*Registration::_Tfinal4x4) = (*Registration::_Tlast)
				* (*Registration::_Tfinal4x4);
	} else {
		retval = NOTMATCHABLE;
	}

	return retval;
}

EnumState Icp::iterate(double* rms, unsigned int* pairs,
		unsigned int* iterations, Matrix* Tinit) {
	if (_trace) {
		_trace->reset();
		_trace->setModel(_model, _sizeModel);
	}

	Registration::_Tfinal4x4->setIdentity();

	if (Tinit) {
		Registration::applyTransformation(_sceneTmp, _sizeScene,
				Registration::_dim, Tinit);
		if (_normalsSTmp)
			Registration::applyTransformation(_normalsSTmp, _sizeScene,
					Registration::_dim, Tinit);
		(*Registration::_Tfinal4x4) = (*Tinit) * (*Registration::_Tfinal4x4);
	}

	EnumState eRetval = PROCESSING;
	unsigned int iter = 0;
	double rms_prev = 10e12;
	unsigned int conv_cnt = 0;
	while (eRetval == PROCESSING) {
		eRetval = step(rms, pairs);
		iter++;

		if (fabs(*rms - rms_prev) < 10e-10)
			conv_cnt++;
		else
			conv_cnt = 0;
		if ((*rms <= _maxRMS || conv_cnt >= _convCnt))
			eRetval = SUCCESS;
		else if (iter >= Registration::_maxIterations)
			eRetval = MAXITERATIONS;

		rms_prev = *rms;
	}
	*iterations = iter;

	return eRetval;
}

EnumState Icp::align(Matrix* Tinit, bool verbose) {
	double rms;
	unsigned int pairs;
	unsigned int it;
	iterate(&rms, &pairs, &it, Tinit);
	if (verbose) {
		cout << endl << "Error: " << _estimator->getRMS() << endl;
		cout << "Iterations: " << _estimator->getIterations() << endl;
	}
}

string getValue(string line) {
	unsigned int start = line.find(": ") + strlen(": ");
	unsigned int end = line.find(",");
	if (end == string::npos) {
		end = line.find("\n");
	}

	return line.substr(start, end - start);

}

int Icp::loadParametersFromXML(string filePath) {

	string algorithm;
	string line;
	ifstream file(filePath.c_str());

	if (file.is_open()) {
		while (getline(file, line)) {
			if (line.find("algorithm") != string::npos) {
				algorithm = getValue(line);
			} else if (line.find("iterations") != string::npos) {
				Registration::_maxIterations = std::atoi(
						getValue(line).c_str());
			} else if (line.find("maxRMS") != string::npos) {
				_maxRMS = atof(getValue(line).c_str());
			} else if (line.find("ConvCnt") != string::npos) {
				_convCnt = atoi(getValue(line).c_str());
			}
		}
	} else {
		cout << "Parsing Registration Parameters: Unable to open file" << endl;
		return -1;
	}
	file.close();

	//Check parameters
	if (algorithm != "ICP" || Registration::_maxIterations <= 0 || _maxRMS < 0
			|| _convCnt < 0) {
		cout << "Parsing Registration Parameters: Parameters invalid." << endl
				<< "ICP intended to use." << endl << "Algorithm is "
				<< algorithm << endl << "Iterations: "
				<< Registration::_maxIterations << " MaxRMS: " << _maxRMS
				<< " ConvCnt: " << _convCnt << endl << endl;
		return -1;
	}

	return 0;
}

void Icp::serializeTrace(char* folder, unsigned int delay) {
	if (_trace)
		_trace->serialize(folder, delay);
}

}
