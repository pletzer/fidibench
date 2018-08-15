/**
 * @brief Write data on a uniform grid to a VTK file
 * @author Alexander Pletzer
 *
 * This software is provided with the hope that it will be 
 * useful. It comes with no warranty whatsoever. 
 * Please send bug reports to alexander@gokliya.net.
 */

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "MultiArrayIter.h"

#ifndef WRITE_VTK_H
#define WRITE_VTK_H

/**
 * Write data to a VTK file
 * @param filename file name
 * @param globalDims global dimsensions
 * @param xmins low corner point of the domain
 * @param xmaxs high corner point of the domain
 * @param field flat array of data (assumed to be in column major order)
 */
void writeVTK(const std::string& filename, 
	const std::vector<size_t>& globalDims,
	const std::vector<double>& xmins, 
	const std::vector<double>& xmaxs, 
	const std::vector<double>& field,
	const std::string& name);

#endif // WRITE_VTK_H
