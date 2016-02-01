# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import pyworkflow.protocol.params as params
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtReconstruct3D
from pyworkflow.em.packages.xmipp3.convert import setOfParticlesToMd, getImageLocation
import xmipp


class XmippProtHelixInitial(ProtReconstruct3D):
    """    
    Reconstruct an helical initial tube from 2d classification.
    """
    _label = 'helical initial'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputClasses2D', params.PointerParam,
                      pointerClass='SetOfClasses2D',
                      label="Input classes", important=True,
                      help='Select the input images from the project.')

        group = form.addGroup('Helical parameters')
        group.addParam('heightFraction', params.FloatParam,
                       label='Height fraction',
                       help="")
        group.addParam('helixParams', params.StringParam,
                       label='Helix params: ',
                       help="")

        form.addParallelSection(threads=2, mpi=0)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_xmd': self._getExtraPath('input_particles.xmd'),
            'initial_volume': self._getPath('initial_volume.vol'),
            'helical_volume': self._getPath('helical_volume.vol')
            }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
        inputClasses = self.inputClasses2D.get()
        inputImages = inputClasses.getImages()
        sampling = inputImages.getSamplingRate()

        self._createFilenameTemplates()
        self._insertFunctionStep('convertInputStep', inputClasses.getObjId())
        self._insertFunctionStep('reconstructStep', sampling)
        self._insertFunctionStep('symmetrizeStep', sampling)
        self._insertFunctionStep('createOutputStep')
        
    #--------------------------- STEPS functions --------------------------------------------
    def createHorizontalMask(self, inputClasses):
        """ Create an horizontal line mask to be used for alignment. """
        firstAvg = inputClasses.getFirstItem().getRepresentative()
        maskFile = self._getExtraPath('horizontal_mask.spi')

        params =  '-i %s' % getImageLocation(firstAvg)
        params += ' --mask rectangular -%d -10' % firstAvg.getXDim()
        params += ' --create_mask %s' % maskFile

        self.runJob('xmipp_transform_mask', params)
        #xmipp_transform_mask "aligned2_ref.xmp" "--create_mask" "mask.spi" "--mask" "rectangular" "-80" "-10"

        mdClasses = xmipp.MetaData()
        for cls2D in inputClasses:
            objId = mdClasses.addObject()
            clsLoc = getImageLocation(cls2D.getRepresentative())
            mdClasses.setValue(xmipp.MDL_IMAGE, clsLoc, objId)
            mdClasses.setValue(xmipp.MDL_ITEM_ID, long(cls2D.getObjId()), objId)

        mdFile = self._getExtraPath('averages.xmd')
        mdClasses.write(mdFile)
        params =  ' -i %s' % mdFile
        params += ' --ref %s' % maskFile
        params += ' --do_not_check_mirrors --iter 1'
        params += ' --oroot %s' % self._getExtraPath('align')

        self.runJob('xmipp_image_align', params)

    def convertInputStep(self, classesId):
        inputClasses = self.inputClasses2D.get()
        self.createHorizontalMask(inputClasses)

        mdAll = xmipp.MetaData()

        for cls2D in inputClasses:
            md = xmipp.MetaData()
            setOfParticlesToMd(cls2D, md)
            # TODO: Determine which class is vertical or horizontal
            # for horizontal classes we need to drop x-shift
            # for vertical classes we need to drop y-shift
            mdAll.unionAll(md)

        mdAll.fillConstant(xmipp.MDL_ANGLE_TILT, 90)
        mdAll.fillRandom(xmipp.MDL_ANGLE_ROT, "uniform", 0, 360)

        #mdAll.fillConstant(xmipp.MDL_ANGLE_ROT, 0)
        #mdAll.fillRandom(xmipp.MDL_ANGLE_TILT, "uniform", 0, 180)

        mdAll.write(self._getFileName('input_xmd'))

    def reconstructStep(self, sampling):
        """ Create the input file in STAR format as expected by Xmipp.
        If the input particles comes from Xmipp, just link the file. 
        """
        params =  '  -i %s' % self._getFileName('input_xmd')
        params += '  -o %s' % self._getFileName('initial_volume')
        params += ' --thr %d' % 1 # TODO: self.numberOfThreads.get()
        params += ' --sampling %f' % sampling

        self.runJob('xmipp_reconstruct_fourier', params)

    def symmetrizeStep(self, sampling):
        params =  '  -i %s' % self._getFileName('initial_volume')
        params += '  -o %s' % self._getFileName('helical_volume')
        params += ' --sampling %f' % sampling
        params += ' --sym helicalDihedral'
        params += ' --heightFraction %f' % self.heightFraction
        params += ' --helixParams %s' % self.helixParams

        self.runJob('xmipp_transform_symmetrize', params)
        # xmipp_transform_symmetrize "-i" "rec_fourier.vol" "--sampling" "4.52" "--sym" "helicalDihedral"
        #  "--heightFraction" "0.63" "--helixParams" "29.535039" "-67.153007" "-o" "rec_fourier_sym.vol"
            
    def createOutputStep(self):
        inputImages = self.inputClasses2D.get().getImages()
        volume = Volume()
        volume.setFileName(self._getFileName('helical_volume'))
        volume.setSamplingRate(inputImages.getSamplingRate())
        
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(self.inputClasses2D, volume)
    
    #--------------------------- INFO functions -------------------------------------------- 
    def _validate(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
