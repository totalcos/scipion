# **************************************************************************
# *
# * Authors:     Josue Gomez Blanco     (josue.gomez-blanco@mcgill.ca) [1]
# *              J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# * [1] McGill University
# * [2] SciLifeLab, Stockholm University
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.protocol.params import (PointerParam, FloatParam, FileParam,
                                        BooleanParam, IntParam, LEVEL_ADVANCED)
from pyworkflow.em.data import Volume
from pyworkflow.em.protocol import ProtAnalysis3D, ImageHandler
import pyworkflow.em.metadata as md
from pyworkflow.utils import exists
from convert import isVersion3


class ProtRelionPostprocess(ProtAnalysis3D):
    """
    Relion post-processing protocol for automated masking,
    overfitting estimation, MTF-correction and B-factor sharpening.
    """
    _label = 'post-processing'
    
    def _createFilenameTemplates(self):
        """ Centralize how files are called for iterations and references. """
        myDict = {
                  'finalVolume': self._getTmpPath("relion_class001.mrc"),
                  'half1': self._getTmpPath("relion_half1_class001_unfil.mrc"),
                  'half2': self._getTmpPath("relion_half2_class001_unfil.mrc"),
                  'mask': self._getTmpPath("input_mask.mrc")
                  }

        self._updateFilenamesDict(myDict)

    # -------------------------- DEFINE param functions -----------------------
    def _defineParams(self, form):
        
        form.addSection(label='Input')
        form.addParam('protRefine', PointerParam, pointerClass="ProtRefine3D",
                      label='Select a previous refinement protocol',
                      help='Select any previous refinement protocol to get the '
                           '3D half maps. Note that it is recommended that the '
                           'refinement protocol uses a gold-standard method.')
        form.addParam('solventMask', PointerParam, pointerClass="VolumeMask",
                      label='Solvent mask',
                      help="Provide a soft mask where the protein is white "
                           "(1) and the solvent is black (0). Often, the "
                           "softer the mask the higher resolution estimates "
                           "you will get. A soft edge of 5-10 pixels is often "
                           "a good edge width.")
        form.addParam('calibratedPixelSize', FloatParam, default=0,
                      label='Calibrated pixel size (A)',
                      help="Provide the final, calibrated pixel size in "
                           "Angstroms. If 0, the input pixel size will be used. "
                           "This value may be different from the pixel-size "
                           "used thus far, e.g. when you have recalibrated "
                           "the pixel size using the fit to a PDB model. "
                           "The X-axis of the output FSC plot will use this "
                           "calibrated value.")

        form.addSection(label='Sharpening')
        form.addParam('mtf', FileParam,
                      label='MTF-curve file',
                      help='User-provided STAR-file with the MTF-curve '
                           'of the detector.'
                           'Relion param: <--mtf>')
        form.addParam('doAutoBfactor', BooleanParam, default=True,
                      label='Estimate B-factor automatically?',
                      help='If set to Yes, then the program will use the '
                           'automated procedure described by Rosenthal and '
                           'Henderson (2003, JMB) to estimate an overall '
                           'B-factor for your map, and sharpen it accordingly.')
        line = form.addLine('B-factor resolution (A): ',
                            condition='doAutoBfactor',
                            help='There are the frequency (in Angstroms), '
                                 'lowest and highest, that will be included in '
                                 'the linear fit of the Guinier plot as '
                                 'described in Rosenthal and Henderson '
                                 '(2003, JMB).')
        line.addParam('bfactorLowRes', FloatParam, default='10.0', label='low')
        line.addParam('bfactorHighRes', FloatParam, default='0.0', label='high')
        form.addParam('bfactor', FloatParam, default=-350,
                      condition='not doAutoBfactor',
                      label='Provide B-factor:',
                      help= 'User-provided B-factor (in A^2) for map '
                            'sharpening, e.g. -400. Use negative values for '
                            'sharpening. Be careful: if you over-sharpen\n'
                            'your map, you may end up interpreting noise for '
                            'signal!\n'
                            'Relion param: *--adhoc_bfac*')
        
        form.addSection(label='Filtering')
        form.addParam('skipFscWeighting', BooleanParam, default=False,
                      label='Skip FSC-weighting for sharpening?',
                      help='If set to No (the default), then the output map '
                           'will be low-pass filtered according to the '
                           'mask-corrected, gold-standard FSC-curve. '
                           'Sometimes, it is also useful to provide an ad-hoc '
                           'low-pass filter (option below), as due to local '
                           'resolution variations some parts of the map may '
                           'be better and other parts may be worse than the '
                           'overall resolution as measured by the FSC. In '
                           'such  cases, set this option to Yes and provide '
                           'an ad-hoc filter as described below.')
        form.addParam('lowRes', FloatParam, default=5,
                      condition='skipFscWeighting',
                      label='Low-pass filter (A):',
                      help='This option allows one to low-pass filter the map '
                           'at a user-provided frequency (in Angstroms). When '
                           'using a resolution that is higher than the '
                           'gold-standard FSC-reported resolution, take care '
                           'not to interpret noise in the map for signal...')
        form.addParam('filterEdgeWidth', IntParam, default=2,
                      expertLevel=LEVEL_ADVANCED,
                      label='Low-pass filter edge width:',
                      help='Width of the raised cosine on the low-pass filter '
                           'edge (in resolution shells)\n'
                           'Relion param: *--filter_edge_width*')
        form.addParam('randomizeAtFsc', FloatParam, default=0.8,
                      expertLevel=LEVEL_ADVANCED,
                      label='Randomize phases threshold',
                      help='Randomize phases from the resolution where FSC '
                           'drops below this value\n'
                           'Relion param: *--randomize_at_fsc*')
        
        form.addParallelSection(threads=0, mpi=0)
    
    # -------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        objId = self.protRefine.get().getObjId()
        self._createFilenameTemplates()
        self._defineParamDict()
        self._insertFunctionStep('convertInputStep', objId)
        self._insertFunctionStep('postProcessStep', self.paramDict)
        self._insertFunctionStep('createOutputStep')
    
    # -------------------------- STEPS functions -------------------------------
    def convertInputStep(self, protId):
        protRef = self.protRefine.get()
        protClassName = protRef.getClassName()

        # FIXME: JMRT I think this a dirty hack to plugin other refinement output here
        # I would suggest better that we standardize the outputs: finalMap, half1 and half2
        # in the base class ProtRefine3D and use the general approach here, instead of
        # dealing with it here as it is now, that is not extensible.

        if protClassName.startswith('ProtRelionRefine3D'):
            protRef._createFilenameTemplates()
            final = protRef._getFileName("finalvolume", ref3d=1)
            half1 = protRef._getFileName("final_half1_volume", ref3d=1)
            half2 = protRef._getFileName("final_half2_volume", ref3d=1)
            
        elif protClassName.startswith('ProtFrealign'):
            protRef._createFilenameTemplates()
            lastIter = protRef._getLastIter()
            final = protRef._getFileName('iter_vol', iter=lastIter)
            half1 = protRef._getFileName('iter_vol1', iter=lastIter)
            half2 = protRef._getFileName('iter_vol2', iter=lastIter)
            
        elif protClassName.startswith('XmippProtProjMatch'):
            iterN = protRef.getLastIter()
            protRef._initialize()
            final = protRef._getFileName('reconstructedFileNamesIters',
                                            iter=iterN, ref=1)
            half1 = protRef._getFileName('reconstructedFileNamesItersSplit1',
                                         iter=iterN, ref=1)
            half2 = protRef._getFileName('reconstructedFileNamesItersSplit2',
                                         iter=iterN, ref=1)

        elif protClassName.startswith('EmanProtRefine'):
            protRef._createFilenameTemplates()
            numRun = protRef._getRun()
            protRef._createIterTemplates(numRun)
            iterN = protRef._lastIter()
            final = protRef._getFileName("mapFull", run=numRun, iter=iterN)
            half1 = protRef._getFileName("mapEvenUnmasked", run=numRun)
            half2 = protRef._getFileName("mapOddUnmasked", run=numRun)
            
        ih = ImageHandler()
        ih.convert(final, self._getFileName("finalVolume"))
        ih.convert(half1, self._getFileName("half1"))
        ih.convert(half2, self._getFileName("half2"))
        ih.convert(self.solventMask.get(), self._getFileName('mask'))
    
    def postProcessStep(self, paramDict):
        params = ' '.join(['%s %s' % (k, str(v))
                           for k, v in self.paramDict.iteritems()])
        self.runJob('relion_postprocess', params)
    
    def createOutputStep(self):
        volume = Volume()
        volume.setFileName(self._getExtraPath('postprocess.mrc'))
        vol = self.protRefine.get().outputVolume
        pxSize = vol.getSamplingRate()
        volume.setSamplingRate(pxSize)
        self._defineOutputs(outputVolume=volume)
        self._defineSourceRelation(vol, volume)

    # -------------------------- INFO functions --------------------------------
    def _validate(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION.
        """
        errors = []
        mtfFile = self.mtf.get()
        
        if mtfFile and not exists(mtfFile):
            errors.append("Missing MTF-file '%s'" % mtfFile)

        protClassName = self.protRefine.get().getClassName()
        if protClassName.startswith('SpiderProtRefinement'):
            errors.append("Relion post-process protocol not implemented for "
                          "Spider - refinement.")
        
        if protClassName.startswith('XmippProtReconstructHighRes'):
            errors.append("Relion post-process protocol not implemented for "
                          "Xmipp - highres.")
        return errors
    
    def _summary(self):
        """ Should be overwritten in subclasses to
        return summary message for NORMAL EXECUTION. 
        """
        summary = []
        postStarFn = self._getExtraPath("postprocess.star")
        if exists(postStarFn):
            mdResol = md.RowMetaData(postStarFn)
            resol = mdResol.getValue(md.RLN_POSTPROCESS_FINAL_RESOLUTION)
            summary.append("Final resolution: *%0.2f A*" % resol)
        
        return summary
        
    # -------------------------- UTILS functions -------------------------------
    def _defineParamDict(self):
        """ Define all parameters to run relion_postprocess"""
        volume = self.protRefine.get().outputVolume
        cps = self.calibratedPixelSize.get()
        angpix = cps if cps > 0 else volume.getSamplingRate()

        # It seems that in Relion3 now the input should be the map
        # filename and not the prefix as before
        if isVersion3():
            inputFn = self._getFileName('half1')
        else:
            inputFn = self._getTmpPath("relion")

        self.paramDict = {'--i': inputFn,
                          '--o': self._getExtraPath('postprocess'),
                          '--angpix': angpix,
                          # Expert params
                          '--filter_edge_width': self.filterEdgeWidth.get(),
                          '--randomize_at_fsc': self.randomizeAtFsc.get()
                          }
        self.paramDict['--mask'] = self._getFileName('mask')

        mtfFile = self.mtf.get()
        if mtfFile:
            self.paramDict['--mtf'] = mtfFile

        if self.doAutoBfactor:
            self.paramDict['--auto_bfac'] = ''
            self.paramDict['--autob_lowres'] = self.bfactorLowRes.get()
            self.paramDict['--autob_highres'] = self.bfactorHighRes.get()
        else:
            self.paramDict['--adhoc_bfac'] = self.bfactor.get()

        if self.skipFscWeighting:
            self.paramDict['--skip_fsc_weighting'] = ''
            self.paramDict['--low_pass'] = self.lowRes.get()
    
    def _getRelionMapFn(self, fn):
        return fn.split(':')[0]
