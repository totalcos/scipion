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
import pyworkflow.em.metadata as md
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
        group.addParam('zetaRise', params.FloatParam,
                       label='Zeta rise (A)',
                       help="Zeta rise in Angstroms.")
        group.addParam('phiAngle', params.FloatParam,
                       label='Phi angle (deg)',
                       help="Phi angle in degrees. Use positive values "
                            "for right handed helixes. Negative left handed. ")
        group.addParam('dihedral', params.BooleanParam,
                      label='Is a dihedral helix?')

        form.addParallelSection(threads=2, mpi=0)

    #--------------------------- INSERT steps functions --------------------------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
            'input_xmd': self._getExtraPath('input_particles.xmd'),
            'initial_volume': self._getExtraPath('initial_volume.vol'),
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

    def convertInputStep(self, classesId):
        inputClasses = self.inputClasses2D.get()
        clsDict = self.detectOrientation(inputClasses)

        mdAll = xmipp.MetaData()

        for cls2D in inputClasses:
            md = xmipp.MetaData()
            setOfParticlesToMd(cls2D, md)
            rot = clsDict[cls2D.getObjId()]

            if rot is None:
                print "WARNING: Ignoring class %s" % cls2D.getObjId()
            else:
                if rot == 90:
                    print "Rotating by 90 images of class %s" % cls2D.getObjId()
                    md.operate("anglePsi=anglePsi+90")
                mdAll.unionAll(md)

        mdAll.fillConstant(xmipp.MDL_ANGLE_TILT, 90)
        mdAll.fillRandom(xmipp.MDL_ANGLE_ROT, "uniform", 0, 360)

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
        params += ' --sym helical%s' % ('Dihedral' if self.dihedral else '')
        params += ' --heightFraction %f' % self.heightFraction
        params += ' --helixParams %f %f' % (self.zetaRise, self.phiAngle)

        self.runJob('xmipp_transform_symmetrize', params)

    def createOutputStep(self):
        inputImages = self.inputClasses2D.get().getImages()
        volume = Volume()
        volume.setFileName(self._getFileName('helical_volume'))
        volume.setSamplingRate(inputImages.getSamplingRate())
        
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(self.inputClasses2D, volume)

        # Let use an special iterator to skip missing particles
        # discarded by classification (in case of cl2d)
        setMdIter = md.SetMdIterator(self._getFileName('input_xmd'),
                                     sortByLabel=md.MDL_ITEM_ID,
                                     updateItemCallback=self._createItemMatrix)
        outputImages = self._createSetOfParticles()
        outputImages.copyInfo(inputImages)
        outputImages.copyItems(inputImages,
                               updateItemCallback=setMdIter.updateItem)
        outputImages.setAlignmentProj()
        self._defineOutputs(outputParticles=outputImages)
        self._defineSourceRelation(self.inputClasses2D, outputImages)
    
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
    def _createItemMatrix(self, item, row):
        from pyworkflow.em.packages.xmipp3.convert import createItemMatrix
        import pyworkflow.em as em
        createItemMatrix(item, row, align=em.ALIGN_PROJ)

    def createHorizontalMask(self, inputClasses):
        """ Create an horizontal line mask to align with
        class averages and detect if there need to be
        retated by 90 degrees.
        """
        firstAvg = inputClasses.getFirstItem().getRepresentative()
        maskFile = self._getExtraPath('horizontal_mask.spi')

        params =  '-i %s' % getImageLocation(firstAvg)
        params += ' --mask rectangular -%d -10' % firstAvg.getXDim()
        params += ' --create_mask %s' % maskFile

        self.runJob('xmipp_transform_mask', params)
        #xmipp_transform_mask "aligned2_ref.xmp" "--create_mask" "mask.spi" "--mask" "rectangular" "-80" "-10"
        return maskFile

    def detectOrientation(self, inputClasses):
        """ Check if the class averages are vertical or horizontal.
        If they are vertical, we need to rotate images by 90 degrees
        to have all of them horizontal before reconstruction.
        """
        maskFile = self.createHorizontalMask(inputClasses)

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

        clsDict = {}
        mdAlign = xmipp.MetaData(self._getExtraPath('align_alignment.xmd'))
        for objId in mdAlign:
            clsId = mdAlign.getValue(xmipp.MDL_ITEM_ID, objId)
            clsImg = mdAlign.getValue(xmipp.MDL_IMAGE, objId)
            absAngle = abs(mdAlign.getValue(xmipp.MDL_ANGLE_PSI, objId))

            if absAngle > 180:
                rot = None
            else:
                rot = 0 if (absAngle < 20 or absAngle > 160) else 90
            clsDict[clsId] = rot
            print " class %d: Rotate %s" % (clsId, rot)

        return clsDict
