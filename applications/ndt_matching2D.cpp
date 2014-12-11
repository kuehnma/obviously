/**
 * Sample application showing the usage of the OBVIOUS NDT implementation.
 * @author Stefan May
 * @date 23.08.2014
 */

#include <string.h>
#include <iostream>

#include "obcore/math/mathbase.h"
#include "obcore/math/linalg/linalg.h"
#include "obcore/base/Timer.h"
#include "obvision/ndt/Ndt.h"
#include "obvision/registration/Registration.h"

using namespace std;
using namespace obvious;

int main(int argc, char** argv) {
	// Model coordinates
	obvious::Matrix* M = new obvious::Matrix(1000, 2);

	for (int i = 0; i < 200; i++) {
		double di = (double) i;
		(*M)(i, 0) = sin(di / 100.0);
		(*M)(i, 1) = sin(di / 20.0);
	}

	obvious::Matrix T = MatrixFactory::TransformationMatrix33(deg2rad(4.0), 0.2,
			0.2);
	obvious::Matrix S = M->createTransform(T);

	//fixme dynamic boarder computation
	Registration* ndt = new Ndt(-10, 10, -10, 10, 0.5);
	ndt->activateTrace();
	ndt->setModel(M);
	ndt->setScene(&S);
	ndt->setMaxIterations(100);

	double rms;
	unsigned int it;
	EnumState state = ndt->iterate(&rms, &it);
	obvious::Matrix F = ndt->getFinalTransformation();
	F.invert();
	cout << "NDT state: " << state << ", with " << it << " iterations): "
				<< endl;

	cout << "Iterations: " << it << endl;

	cout << "Applied transformation:" << endl;
	T.print();

	char folder[6] = "trace";
	ndt->serializeTrace(folder, 10);

	cout << endl << "Estimated transformation:" << endl;
	F.print();

	delete ndt;

	return 0;
}
