/***************************************************************************
 * Authors:     AUTHOR_NAME (jvargas@cnb.csic.es)
 *
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

#ifndef VALIDATION_NONTILT_H_
#define VALIDATION_NONTILT_H_
#define PI 3.14159265

#include <data/xmipp_program.h>
#include "reconstruct_significant.h"
#include <math.h>

/**@defgroup Validation without tilt
   @ingroup ReconsLibrary */
//@{
class ProgValidationNonTilt: public XmippProgram
{


public:
    /** Filenames */
    FileName fnDir, fnSym, fnInit;

    /** Sampling rate of the volume and projections */
    double sampling_rate;


public:

    void readParams();

    void defineParams();

    void run();

public:
    void obtainSumU(const MetaData & tempMd,std::vector<double> & sum_u,std::vector<double> & H0);

    void obtainSumW(const MetaData & tempMd,std::vector<double> & sum_W,std::vector<double> & sum_u,std::vector<double> & H, const double factor);


};
#endif /* VALIDATION_NONTILT_H_ */
