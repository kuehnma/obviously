#ifndef NDTTRACE_H_
#define NDTTRACE_H_

#include <iostream>
#include <vector>
#include "obvision/icp/assign/assignbase.h"
#include "obcore/math/linalg/linalg.h"

using namespace std;

namespace obvious
{

/**
 * @class NdtTrace
 * @brief Represents a trace module for the ndt
 * @author Stefan May (IcpTrace), changed by Markus Kühn
 **/
class NdtTrace
{
public:
	/**
	 * Default constructor
	 */
	NdtTrace(unsigned int dim);
		 
	/**
	 * Destructor
	 */
	~NdtTrace();
	
	/**
	 * Reset trace record
	 */
	void reset();

	/**
   * Set model of trace record
   * @param model model data
   */
	void setModel(Matrix* model);

	/**
	 * Add scene assignment to trace record
	 * @param scene scene data
	 * @param sizeS size of scene
	 * @param pairs tuple of assigned indices
	 */
	void addAssignment(double** scene, unsigned int sizeS);
	
	/**
	 * Serialize assignment to trace folder
	 * @param folder trace folder (must not be existent)
	 * @param delay animation delay (specified in delay*1/100s)
	 */
	void serialize(char* folder, unsigned int delay);

private:
	
	unsigned int _dim;

	Matrix* _M;
	
	vector<Matrix*> _scenes;
	
};

}

#endif /*NDTTRACE_H_*/
