# ******************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [1]
# *
# * [1] SciLifeLab, Stockholm University
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
# ******************************************************************************

import pyworkflow.protocol.params as params
import pyworkflow.em as em
import pyworkflow.em.metadata as md

import convert


class ProtRelionCtfRefinement(em.ProtParticles):
    """
    Wrapper protocol for the Relion's per-particle CTF refinement.
    """

    _label = 'ctf refinement'

    def _defineParams(self, form):
        form.addSection(label='Input')
        #TODO: conditions on particles?
        form.addParam('inputParticles', params.PointerParam,
                      important=True,
                      label='Input particles',
                      pointerClass='SetOfParticles',
                      help='Provide a set of particles for local CTF refinement.')
        form.addParam('inputPostprocess', params.PointerParam,
                      important=True,
                      label='Input Postprocess',
                      pointerClass='ProtRelionPostprocess',
                      help='Select a PostProcess job. The mask used for this '
                           'postprocessing will be applied to the unfiltered '
                           'half-maps and should encompass the entire complex. '
                           'The resulting FSC curve will be used for weighting '
                           'the different frequencies. \n\n'
                           'Note that for helices it is common practice to use '
                           'a mask only encompassing the central 30% or so of '
                           'the box. This gives higher resolution estimates, '
                           'as it disregards ill-defined regions near the box '
                           'edges. However, for ctf_refine it is better to '
                           'use a mask encompassing (almost) the entire box, '
                           'as otherwise there may not be enough signal. ')

        form.addSection(label='Fit')
        form.addParam('minResolution', params.FloatParam, default=30,
                      label='Minimum resolution for fits (A)',
                      help="The minimum spatial frequency (in Angstrom) used "
                           "in the beamtilt fit.")

        form.addParam('doCtfFitting', params.BooleanParam, default=True,
                      label='Perform CTF parameter fitting?',
                      help='If set to Yes, then relion_ctf_refine will be '
                           'used to estimate the selected parameters below.')
        form.addParam('fitPartDefocus', params.BooleanParam, default=True,
                      condition='doCtfFitting',
                      label='Fit per-particle defocus?',
                      help='If set to Yes, then relion_ctf_refine will '
                           'estimate a per-particle defocus.')
        form.addParam('rangeDefocusFit', params.IntParam, default=2000,
                      condition='doCtfFitting and fitPartDefocus',
                      label='Range for defocus fit (A)',
                      help='The range in (Angstrom) for the defocus fit of '
                           'each particle.')
        form.addParam('fitAstig', params.EnumParam, default=0,
                      choices=['no', 'per-micrograph', 'per-particle'],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Fit astigmatism? ',
                      help="If *per-micrograph*, ctf_refine will try to refine "
                           "astigmatism on a per-micrograph basis. This will "
                           "require many particles and good signal-to-noise "
                           "ratios per micrograph.\n\n"
                           "If *per-particle*, astigmatism will be estimated "
                           "on a per-particle basis. This requires very strong "
                           "data, i.e. very large particles with excellent "
                           "signal-to-noise ratios.")
        form.addParam('fitMicPhaseShift', params.BooleanParam, default=False,
                      condition='doCtfFitting',
                      label='Fit per-micrograph phase-shift?',
                      help='If set to Yes, ctf_refine will try to refine a '
                           'phase-shift (amplitude contrast) on a per-'
                           'micrograph basis. This may be useful for Volta-'
                           'phase plate data, but will require many particles '
                           'and good signal-to-noise ratios per micrograph.')
        form.addParam('doBeamtiltEstimation', params.BooleanParam, default=False,
                      label='Perform beamtilt estimation?',
                      help='If set to Yes, then relion_ctf_refine will also '
                           'estimate the beamtilt over the entire data set. '
                           'This option is only recommended for high-resolution '
                           'data sets, i.e. significantly beyond 3 Angstrom '
                           'resolution.')

        form.addParallelSection(threads=1, mpi=1)

    # -------------------------- STEPS functions -------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('refineCtfStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        inputParts = self.inputParticles.get()
        imgStar = self._getPath('input_particles.star')

        self.info("Converting set from '%s' into '%s'" %
                  (inputParts.getFileName(), imgStar))

        convert.writeSetOfParticles(inputParts, imgStar, self._getExtraPath(),
                                    alignType=em.ALIGN_PROJ,
                                    fillRandomSubset=True)

    def _getInputVolumes(self, postStar):
        """ Parse the input volumes: halfs and mask
        from the postprocess.star file. """
        half1 = half2 = mask = ''
        with open(postStar) as f:
            for line in f:
                if line.startswith('_rlnUnfilteredMapHalf1'):
                    half1 = line.split()[1]
                if line.startswith('_rlnUnfilteredMapHalf2'):
                    half2 = line.split()[1]
                if line.startswith('_rlnMaskName'):
                    mask = line.split()[1]
                    break
        return half1, half2, mask

    def refineCtfStep(self):
        args = "--i %s " % self._getPath('input_particles.star')
        args += "--o %s " % self._getExtraPath()
        postStar = self.inputPostprocess.get()._getExtraPath('postprocess.star')
        args += "--f %s " % postStar
        args += "--m1 %s --m2 %s --mask %s " % self._getInputVolumes(postStar)
        args += "--kmin_tilt %0.3f " % self.minResolution

        if self.doCtfFitting:
            args += "--fit_defocus "

        fitAstig = self.fitAstig.get()
        if fitAstig == 1:
            args += "--glob_astig "
        elif fitAstig == 2:
            args += "--astig "

        if self.fitMicPhaseShift:
            args += "--fit_phase "

        if self.doBeamtiltEstimation:
            args += "--fit_beamtilt "

        args += "--j %d " % self.numberOfThreads

        self.runJob("relion_ctf_refine", args)

    def createOutputStep(self):
        imgSet = self.inputParticles.get()
        outImgSet = self._createSetOfParticles()
        outImgSet.copyInfo(imgSet)

        outImgsFn = self._getExtraPath('particles_ctf_refine.star')
        imgSet.setAlignmentProj()
        rowIterator = md.iterRows(outImgsFn, sortByLabel=md.RLN_IMAGE_ID)
        outImgSet.copyItems(imgSet,
                            updateItemCallback=self._updateItemCtf,
                            itemDataIterator=rowIterator)

        self._defineOutputs(outputParticles=outImgSet)
        self._defineTransformRelation(self.inputParticles, outImgSet)

    def _updateItemCtf(self, particle, row):
        particle.setCTF(convert.rowToCtfModel(row))
        #TODO: Add other field from the .star file when other options?

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        errors = []
        return errors

