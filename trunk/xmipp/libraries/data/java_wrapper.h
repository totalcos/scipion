/***************************************************************************
 *
 * Author: Jesus Cuenca (jcuenca@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef JAVA_WRAPPER_H
#define JAVA_WRAPPER_H

#include <string>
#include <iostream>
#include <unistd.h>
#include <error.h>
#include "image.h"
#include "filename.h"

using namespace std;

// abstraction of ImageJ's ImagePlus simplified for C++
class ImagePlusC {
public:
	// filename is a string as returned by metadata.getStrFromValue - see MetaData for details
	std::string filename;
	// array of pixels
	double * image;
	int width, height;
	bool readHeaderOnly;
	int nImages;
	int slice; // which slice to read

	// empty constructor
	ImagePlusC() {
		filename = "";
		image = NULL;
		readHeaderOnly = false;
		width = height = 0;
		nImages = 0;
		slice = -1;
	}

	~ImagePlusC() {
		if (image != NULL)
			delete image;
	}

	// x=0..width-1, y=0..height-1
	double getPixel(int x, int y) {
		return image[(y * width) + x];
	}
};

/**
 * Read the image file specified in ip and save all the relevant info into ip too
 */
int readImage(ImagePlusC &ip) {
	Image<float> in;
	int error = 0;
	try {
		//cout<<ip.filename;
		/*char cCurrentPath[FILENAME_MAX];
		 getcwd(cCurrentPath, sizeof(cCurrentPath));
		 cout << cCurrentPath;*/

		FileName fn = ip.filename;
		// cout << fn;
		// read header / whole file (depending on ip.readHeaderOnly)
		error = in.read(fn, !ip.readHeaderOnly, ip.slice);

		ip.width = in.data.xdim;
		ip.height = in.data.ydim;
		ip.nImages = NSIZE(in());

		/* cout << "readdata:" << !ip.readHeaderOnly << ", w:" << ip.width << ", h:"
				<< ip.height << ", slice:" << ip.slice << "\n";
		cout << "#" << ip.nImages << "\n"; */

		if (ip.readHeaderOnly == false) {
			ip.image = new double[ip.width * ip.height];
			for (int i = 0; i < ip.height; i++) {
				for (int j = 0; j < ip.width; j++) {
					ip.image[(ip.height - 1 - i) * ip.width + ip.width - 1 - j]
							= in.getPixel(i, j);
				}
			}
		}

	} catch (XmippError xe) {
		std::cerr << xe;
		error = 1;
	} catch (std::exception& e) {
		std::cerr << e.what();
	} catch (...) {
		std::cerr << "Unhandled exception";
	}

	return error;
}

Image<double> readFullImage(const std::string &filename) {
	try {
		FileName fn(filename);
		Image<double> img;

		img.read(fn);

		return img;
	} catch (XmippError xe) {
		std::cerr << xe;
	} catch (std::exception& e) {
		std::cerr << e.what();
	} catch (...) {
		std::cerr << "Unhandled exception";
	}

	return Image<double>();
}

Image<double> readImageHeader(const std::string &filename) {
	try {
		FileName fn(filename);
		Image<double> img;

		cout << "FileName: " << fn << endl;

		img.read(fn, false);

		return img;
	} catch (XmippError xe) {
		std::cerr << xe;
	} catch (std::exception& e) {
		std::cerr << e.what();
	} catch (...) {
		std::cerr << "Unhandled exception";
	}

	return Image<double>();
}


#endif
