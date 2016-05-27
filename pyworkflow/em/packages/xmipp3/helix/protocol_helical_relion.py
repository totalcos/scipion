# **************************************************************************
# *
# * Authors:    Jaime Martin-Benito
#               J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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

from os.path import join
from glob import glob

import pyworkflow.protocol.params as params
import pyworkflow.em.metadata as md
import pyworkflow.utils as pwutils

from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtReconstruct3D
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.packages.xmipp3.convert import setOfParticlesToMd, getImageLocation
import pyworkflow.em.packages.relion as relion
from pyworkflow.em.packages.relion.constants import ANGULAR_SAMPLING_LIST
import xmipp


class XmippProtHelixRelion(ProtReconstruct3D):
    """    
    This protocol generates the first iteration to perform a helical
    reconstruction using RELION
    PURPOSE: 3D reconstruction and classification of
    helical structures using RELION.
    """
    _label = 'helical relion'
    
    #--------------------------- DEFINE param functions --------------------------------------------   
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('doContinue', params.BooleanParam, default=False,
                      label='Continue from previous run?')

        form.addParam('previousRun', params.PointerParam,
                      pointerClass='XmippProtHelixRelion',
                      condition='doContinue',
                      label='Previous run')

        form.addParam('inputParticles', params.PointerParam,
                      pointerClass='SetOfParticles',
                      condition='not doContinue',
                      label="Input particles", important=True,
                      help='Select the input images from the project.')

        form.addParam('referenceVolume', params.PointerParam,
                      pointerClass='Volume',
                      condition='not doContinue',
                      label="Input volume", important=True,
                      help='Initial reference 3D map, it should have the same '
                           'dimensions and the same pixel size as your input particles.')

        form.addParam('isMapAbsoluteGreyScale', params.BooleanParam,
                      default=False,
                      condition='not doContinue',
                      label="Is initial 3D map on absolute greyscale?",
                      help='The probabilities are based on squared differences, '
                           'so that the absolute grey scale is important. \n'
                           'Probabilities are calculated based on a Gaussian noise model,'
                           'which contains a squared difference term between the reference and the experimental image. '
                           'This has a consequence that the reference needs to be on the same absolute intensity '
                           'grey-scale as the experimental images. RELION and XMIPP reconstruct maps at their absolute '
                           'intensity grey-scale. Other packages may perform internal normalisations of the reference '
                           'density, which will result in incorrect grey-scales. Therefore: if the map was reconstructed '
                           'in RELION or in XMIPP, set this option to Yes, otherwise set it to No. If set to No, RELION '
                           'will use a (grey-scale invariant) cross-correlation criterion in the first iteration, and '
                           'prior to the second iteration the map will be filtered again using the initial low-pass filter. '
                           'This procedure is relatively quick and typically does not negatively affect the outcome of the '
                           'subsequent MAP refinement. Therefore, if in doubt it is recommended to set this option to No.')

        form.addParam('numberOfReferences', params.IntParam, default=1,
                      condition='not doContinue',
                      label='Number of references',
                      help='')

        form.addParam('totalIters', params.IntParam, default=1,
                      condition='doContinue',
                      label='Total number of iterations',
                      help='The number of iterations will be from last iteration '
                           'of previous protocol until this number.')

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

        line = group.addLine('Symmetry search radius (px)')
        line.addParam('radMin', params.IntParam, default=2, label='min')
        line.addParam('radMax', params.IntParam, default=40, label='min')

        form.addSection(label='Relion Parameters')

        form.addParam('regularisationParamT', params.IntParam, default=2,
                      label='Regularisation parameter T',
                      help='Regularisation parameter (Tau value)')

        form.addParam('referenceMask', params.PointerParam,
                      pointerClass='VolumeMask', allowsNull=True,
                      label='Reference mask (optional)',
                      help='A volume mask containing a (soft) mask with the same dimensions '
                           'as the reference(s), and values between 0 and 1, with 1 being 100% protein '
                           'and 0 being 100% solvent. The reconstructed reference map will be multiplied '
                           'by this mask. If no mask is given, a soft spherical mask based on the <radius> '
                           'of the mask for the experimental images will be applied.\n\n'
                           'In some cases, for example for non-empty icosahedral viruses, it is also useful '
                           'to use a second mask. Check _Advaced_ level and select another volume mask')

        form.addParam('maskDiameterA', params.IntParam, default=-1,
                      condition='not doContinue',
                      label='Particle mask diameter (A)',
                      help='The experimental images will be masked with a soft circular mask '
                           'with this <diameter>. '
                           'Make sure this diameter is not set too small because that may mask '
                           'away part of the signal! If set to a value larger than the image '
                           'size no masking will be performed.\n\n'
                           'The same diameter will also be used for a spherical mask of the '
                           'reference structures if no user-provided mask is specified.')

        form.addParam('doCTF', params.BooleanParam, default=True,
                      label='Do CTF-correction?', condition='not doContinue',
                      help='If set to Yes, CTFs will be corrected inside the MAP refinement. '
                           'The resulting algorithm intrinsically implements the optimal linear, '
                           'or Wiener filter. Note that input particles should contains CTF parameters.')

        form.addParam('haveDataBeenPhaseFlipped', params.LabelParam,
                      label='Have data been phase-flipped?      (Don\'t answer, see help)',
                      help='The phase-flip status is recorded and managed by Scipion. \n'
                           'In other words, when you import or extract particles, \n'
                           'Scipion will record whether or not phase flipping has been done.\n\n'
                           'Note that CTF-phase flipping is NOT a necessary pre-processing step \n'
                           'for MAP-refinement in RELION, as this can be done inside the internal\n'
                           'CTF-correction. However, if the phases have been flipped, the program\n'
                           'will handle it.')

        form.addParam('doNorm', params.BooleanParam, default=False,
                      label='Do normalisation?',
                      help='')

        form.addParam('iniHigh', params.IntParam, default=60,
                      condition='not doContinue',
                      label='Resolution limit (A)',
                      help='Resolution (in Angstroms) to which to limit '
                           'refinement in the first iteration')

        group = form.addGroup('Sampling')
        group.addParam('angularSamplingDeg', params.EnumParam, default=2,
                      choices=ANGULAR_SAMPLING_LIST,
                      label='Angular sampling interval (deg)',
                      help='There are only a few discrete angular samplings possible because '
                       'we use the HealPix library to generate the sampling of the first '
                       'two Euler angles on the sphere. The samplings are approximate numbers '
                       'and vary slightly over the sphere.')
        group.addParam('oversampling', params.IntParam, default=1,
                       label='Oversampling',
                       help='Adaptive oversampling order to speed-up calculations. \n'
                            'Default 1 (0=no oversampling, 1=2x, 2=4x, etc...)')

        group.addParam('offsetSearchRangePix', params.FloatParam, default=5,
                      label='Offset search range (pix)',
                      help='Probabilities will be calculated only for translations in a circle '
                           'with this radius (in pixels). The center of this circle changes at '
                           'every iteration and is placed at the optimal translation for each '
                           'image in the previous iteration.')
        group.addParam('offsetSearchStepPix', params.FloatParam, default=1.0,
                      label='Offset search step (pix)',
                      help='Translations will be sampled with this step-size (in pixels). '
                           'Translational sampling is also done using the adaptive approach. '
                           'Therefore, if adaptive=1, the translations will first be evaluated'
                           'on a 2x coarser grid.')

        group.addParam('limitTilt', params.FloatParam, default=89,
                       label='Limit tilt angle',
                       help='Limited tilt angle: positive for keeping side view. Default -91s'
                            '(i.e. 85 will keep aprox. tilt angles between 85 and 90)'
                            'negative for keeping top viewsLimit tilt angle.'
                            '(i.e. -10 will keep aprox. tilt angles between 0 and 10)')

        self.rootName = 'rnp'

        form.addParallelSection(threads=2, mpi=4)

    #--------------------------- INSERT steps functions ------------------------

    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        iterDir = self._getExtraPath('iter_%(iter)03d')
        iterRoot = join(iterDir, self.rootName)

        myDict = {
            'solventMask': self._getExtraPath('solvent_mask.mrc'),
            'iterDir': iterDir,
            'iterRoot': iterRoot,
            'optimiser': iterRoot + '_it%(iter)03d_optimiser.star',
            'relionVol': iterRoot + '_it%(iter)03d_class%(ref)03d.mrc',
            'c1Vol': iterRoot + '_it%(iter)03d_class%(ref)03d_C1.mrc',
            'symParams': iterRoot + '_it%(iter)03d_class%(ref)03d_params.xmd'
            }
        self._updateFilenamesDict(myDict)

    def _insertAllSteps(self):
        # Update the filenames dictionary
        self._createFilenameTemplates()
        # Store the sampling rate for use in commands
        self.sampling = self.getInputParticles().getSamplingRate()

        # This protocol could run a first iteration (not doContinue)
        # or continue from a previous run with an iterative refinement
        if not self.doContinue:
            self._insertStepsForIter(1)
        else:
            prevIter = self.previousRun.get().getLastIter()
            for i in range(prevIter+1, self.totalIters.get() + 1):
                self._insertStepsForIter(i)

        self._insertFunctionStep('createOutputStep')
        
    #--------------------------- STEPS functions -------------------------------

    def convertInputStep(self, classesId):
        ih = ImageHandler()
        if self.doContinue:
            pass
            # TODO: Copy the last iteration folder from previous run?
        else:
            ih.convert(self.inputReference.get(),
                       self._getFileName('solventMask'))

    def _insertStepsForIter(self, iteration):
        self._insertFunctionStep('prepareIterStep', iteration)
        self._insertFunctionStep('relionStep', iteration)
        self._insertFunctionStep('findSymmetryStep', iteration)
        self._insertFunctionStep('symmetrizeStep', iteration)


    def prepareIterStep(self, iteration):
        """ Prepare folder and files for a given iteration. """
        # Create the folder for this iteration
        pwutils.makePath(self._getFileName('iterDir', iter=iteration))

    def printJob(self, program, args, **kwargs):
        print "program: %s" % program
        print "   args: %s" % args

    def relionStep(self, iteration):
        """ Run relion classify with only one iteration. """
        self.info("Getting optimiser file from the previous iteration: %s" %
                  self._getFileName('optimiser', iter=iteration-1))

        args = " --o %s" % self._getFileName('iterRoot', iter=iteration)
        args += " --continue %s" % self._getFileName('optimiser',
                                                     iter=iteration-1)
        args += " --iter %s" % iteration
        args += " --tau2_fudge %d" % self.regularisationParamT
        args += " --solvent_mask %s" % self._getFileName('solventMask')
        args += " --oversampling %d" % self.oversampling
        args += " --healpix_order %f" % self.angularSamplingDeg
        args += " --offset_range %d" % self.offsetSearchRangePix.get()
        args += " --offset_step %d" % self.offsetSearchStepPix.get() * 2

        self.printJob('relion_refine_mpi', args, env=relion.getEnviron())

    def findSymmetryStep(self, iteration):
        self.info("Searching helical symmetry using XMIPP")

        for ref in range(1, self.numberOfReferences.get()+1):
            args =  '  -i %s' % self._getFileName('relionVol',
                                                  iter=iteration, ref=ref)
            args += '  -o %s' % self._getFileName('symParams',
                                                  iter=iteration, ref=ref)
            args += ' --sampling %f' % self.sampling
            args += ' --sym helical%s' % ('Dihedral' if self.dihedral else '')
            args += ' --heightFraction %f' % self.heightFraction
            args += ' --localHelical %f %f' % (self.zetaRise, self.phiAngle)
            args += '  -z 1.000000 500.000000 1 --rotHelical -360.000000 360.000000 1  --mask tube -<rad_min> -<rad_max> -<image_size> '

            self.printJob('xmipp_transform_symmetrize', args)

    def symmetrizeStep(self, iteration):
        self.info("Imposing helical symmetry using XMIPP")

        for ref in range(1, self.numberOfReferences.get()+1):
            relionVol = self._getFileName('relionVol',
                                          iter=iteration, ref=ref)
            c1Vol = self._getFileName('c1Vol',
                                      iter=iteration, ref=ref)
            # We store the original relion output volume as C1_*
            # Then we symmetrize that volume with Xmipp and
            # replace the original one, so for next iteration
            # we trick Relion input for further refinement
            #pwutils.moveFile(relionVol, c1Vol)
            print "move file: %s -> %s" % (relionVol, c1Vol)
            args =  '  -i %s' % c1Vol
            args += '  -o %s' % relionVol
            args += ' --sampling %f' % self.sampling
            args += ' --sym helical%s' % ('Dihedral' if self.dihedral else '')
            args += ' --heightFraction %f' % self.heightFraction
            args += ' --helixParams %f %f' % (self.zetaRise, self.phiAngle)

            self.printJob('xmipp_transform_symmetrize', args)

    def createOutputStep(self):
        return

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
        errors = []

        if self.doContinue:
            lastIter = self.previousRun.get().getLastIter()
            if self.totalIters <= lastIter:
                errors.append('The total number of iterations should be greater'
                              ' than %s \n(last iteration of previous run.)' %
                              lastIter)

        return errors

    def _summary(self):
        """ Should be overriden in subclasses to 
        return summary message for NORMAL EXECUTION. 
        """
        return []
    
    #--------------------------- UTILS functions --------------------------------------------
    def getInputParticles(self):
        """ Return the input particles, either directly if not continue
         or the ones coming from previous protocol.
        """
        if self.doContinue:
            return self.previousRun.get().getInputParticles()
        else:
            return self.inputParticles.get()

    def getLastIter(self):
        """ Find last completed iteration. """
        iters = glob(self._getExtraPath('iter_???'))
        iters.sort()

        if iters:
            return int(iters[-1].split('iter_')[1])

        return None
