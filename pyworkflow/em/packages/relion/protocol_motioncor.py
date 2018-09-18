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

import os
from itertools import izip
from math import ceil

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
import pyworkflow.protocol.constants as cons
import pyworkflow.utils as pwutils
import pyworkflow.em as em
from pyworkflow.em.protocol import ProtAlignMovies
from pyworkflow.gui.plotter import Plotter
from pyworkflow.protocol import STEPS_SERIAL


class ProtRelionMotioncor(ProtAlignMovies):
    """
    Wrapper protocol for the Relion's implementation of motioncor algorithm.
    """

    _label = 'motioncor'

    def __init__(self, **kwargs):
        ProtAlignMovies.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif'] else 'mrc'

    # -------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):

        line = form.addLine('Frames for corrected SUM',
                             help='First and last frames to use in corrected '
                                  'average (starts counting at 1 and 0 as last '
                                  'means util the last frame in the movie). ')
        line.addParam('sumFrame0', params.IntParam, default=1,
                      label='from')
        line.addParam('sumFrameN', params.IntParam, default=0,
                      label='to')

        form.addParam('doDW', params.BooleanParam, default=False,
                       label='Do dose-weighting?',
                       help='If set to Yes, the averaged micrographs will be '
                            'dose-weighted. \n\n'
                            'NOTE: In Scipion the Voltage and and Dose '
                            'information is provided during import, so you '
                            'do not need to provide them anymore. ')

        form.addParam('saveNonDW', params.BooleanParam, default=False,
                       condition='doDW',
                       label='Save non-dose weighted as well?',
                       help='Aligned but non-dose weighted images are '
                            'sometimes useful in CTF estimation, although '
                            'there is no difference in most cases. Whichever '
                            'the choice, CTF refinement job is always done on '
                            'dose-weighted particles.')

        group = form.addGroup("Motion")

        group.addParam('bfactor', params.IntParam, default=150,
                      label='Bfactor',
                      help="The B-factor that will be applied to the "
                           "micrographs.")

        line = group.addLine('Number of patches',
                            help='Number of patches to be used for patch based '
                                 'alignment. Set to *0 0* to do only global motion '
                                 'correction. \n')
        line.addParam('patchX', params.IntParam, default=1, label='X')
        line.addParam('patchY', params.IntParam, default=1, label='Y')

        group.addParam('groupFrames', params.IntParam, default=1,
                      label='Group frames',
                      help="The B-factor that will be applied to the "
                           "micrographs.")

        group.addParam('binFactor', params.FloatParam, default=1.,
                       label='Binning factor',
                       help='Bin the micrographs this much by a windowing '
                            'operation in the Fourier Tranform. Binning at '
                            'this level is hard to un-do later on, but may be '
                            'useful to down-scale super-resolution images. '
                            'Float-values may be used. Do make sure though '
                            'that the resulting micrograph size is even.')

        group.addParam('gainRotation', params.EnumParam, default=0,
                       choices=['No rotation (0)',
                                ' 90 degrees (1)',
                                '180 degrees (2)',
                                '270 degrees (3)'],
                      label='Gain rotation',
                      help="Rotate the gain reference by this number times 90 "
                           "degrees clockwise in relion_display. This is the "
                           "same as -RotGain in MotionCor2. \n"
                           "Note that MotionCor2 uses a different convention "
                           "for rotation so it says 'counter-clockwise'.")

        group.addParam('gainFlip', params.EnumParam, default=0,
                       choices=['No flipping        (0)',
                                'Flip upside down   (1)',
                                'Flip left to right (2)'],
                      label='Gain flip',
                      help="Flip the gain reference after rotation. "
                           "This is the same as -FlipGain in MotionCor2. "
                           "0 means do nothing, 1 means flip Y (upside down) "
                           "and 2 means flip X (left to right).")

        form.addParam('extraParams', params.StringParam, default='',
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Additional parameters',
                      help="Extra parameters for motioncorr (NOT motioncor2)")

        form.addParam('doComputePSD', params.BooleanParam, default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label="Compute PSD (before/after)?",
                      help="If Yes, the protocol will compute for each movie "
                           "the average PSD before and after alignment, "
                           "for comparison")

        form.addParam('doComputeMicThumbnail', params.BooleanParam,
                      default=False,
                      expertLevel=cons.LEVEL_ADVANCED,
                      label='Compute micrograph thumbnail?',
                      help='When using this option, we will compute a '
                           'micrograph thumbnail and keep it with the '
                           'micrograph object for visualization purposes. ')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions -------------------------------
    def _processMovie(self, movie):
        """ `which relion_run_motioncorr` --i Import/job001/movies.star
        --o MotionCorr/job002/ --first_frame_sum 1 --last_frame_sum 0
        --bin_factor 2 --bfactor 150
        --angpix 0.98 --patch_x 1 --patch_y 1 --gainref gain.mrc
        --gpu "0" --dose_weighting --voltage 300 --dose_per_frame 1
        --preexposure 0
        """
        movieFolder = self._getOutputMovieFolder(movie)
        inputStar = os.path.join(movieFolder,
                                 '%s_input.star' % self._getMovieRoot(movie))
        self.writeInputStar(inputStar, movie)

        pwutils.makePath(movieFolder, 'output')
        # The program will run in the movie folder, so let's put
        # the input files relative to that
        args = "--i %s --o output/ " % os.path.basename(inputStar)
        args += "--use_motioncor2 --use_own --motioncor2_exe fake_mc2 "
        f0, fN = self._getRange(movie)
        args += "--first_frame_sum %d --last_frame_sum %d " % (f0, fN)
        args += "--bin_factor %f --bfactor %d " % (self.binFactor, self.bfactor)
        args += "--angpix %f " % (movie.getSamplingRate())
        args += "--patch_x %d --patch_y %d " % (self.patchX, self.patchY)
        args += "--j %d " % self.numberOfThreads

        if self.doDW:
            preExp, dose = self._getCorrectedDose(self.inputMovies.get())
            voltage = movie.getAcquisition().getVoltage()
            args += "--dose_weighting "
            if self.saveNonDW:
                args += " --save_noDW "
            args += "--voltage %d " % voltage
            args += "--dose_per_frame %f " % dose
            args += "--preexposure %f " % preExp

        self.runJob("relion_run_motioncorr", args, cwd=movieFolder)

        self._computeExtra(movie)

        self._moveFiles(movie)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        # Check base validation before the specific ones for Motioncorr
        errors = ProtAlignMovies._validate(self)
        return errors

    # ------------------------ Extra BASE functions ---------------------------
    def _getRelPath(self, baseName, refPath):
        return os.path.relpath(self._getExtraPath(baseName), refPath)

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _createOutputMovies(self):
        return False

    def _createOutputMicrographs(self):
        return not self.doDW or bool(self.saveNonDW)

    def _createOutputWeightedMicrographs(self):
        return bool(self.doDW)

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setPlotInfo(movie, mic)
        self._setMotionValues(movie, mic)

    def _setMotionValues(self, movie, mic):
        """ Parse motion values from the 'corrected_micrographs.star' file
        generated for each movie. """
        fn = os.path.join(self._getOutputMovieFolder(movie), 'output',
                          'corrected_micrographs.star')
        # FIXME: There are a few quick and dirty solutions here:
        # 1) we are reading the star file by "hand"
        # 2) assuming the data line always starts with 'output/',
        #   since we use that folder as output
        # 3) current implementation only works for a single item
        #   so, a re-write is needed when processing in batch
        #   and corrected_micrographs.star file can contain more rows.
        # Assume file format is the following:
        #
        # data_
        #
        # loop_
        # _rlnMicrographName  # 1
        # _rlnMicrographMetadata  # 2
        # _rlnAccumMotionTotal  # 3
        # _rlnAccumMotionEarly  # 4
        # _rlnAccumMotionLate  # 5
        # output / cct_1.mrc output / cct_1.star   5.347620 4.311143 1.036477
        with open(fn) as f:
            for line in f:
                l = line.strip()
                if l.startswith('output/'):
                    parts = l.split()
                    mic._rlnAccumMotionTotal = pwobj.Float(parts[2])
                    mic._rlnAccumMotionEarly = pwobj.Float(parts[3])
                    mic._rlnAccumMotionLate = pwobj.Float(parts[4])
                    break

    def _getMovieShifts(self, movie):
        outStar = self._getMovieOutFn(movie, '.star')
        first, last = self._getRange(movie)
        n = last - first + 1
        # JMRT: I will parse the shift manually here from the .star file
        # to avoid the need to include new labels and depend on binding code
        # The following structure of the block is assumed:
        """
        data_global_shift

        loop_
        _rlnMicrographFrameNumber #1
        _rlnMicrographShiftX #2
        _rlnMicrographShiftY #3
                   1     0.000000     0.000000
                   2     -0.91811     1.351012
                   ...
        """
        xShifts, yShifts = [], []

        with open(outStar) as f:
            # Locate the desired data block
            found_data = False
            found_loop = False

            for line in f:
                l = line.strip()
                if found_data:
                    if found_loop:
                        if not l:
                            break
                        if not l.startswith('_rln'):
                            parts = l.split()
                            xShifts.append(float(parts[1]))
                            yShifts.append(float(parts[2]))
                            # Avoid reading values from non aligned frames
                            if len(xShifts) == n:
                                break
                    else:
                        if l == 'loop_':
                            found_loop = True
                else:
                    if l == 'data_global_shift':
                        found_data = True

        return xShifts, yShifts

    # --------------------------- UTILS functions -----------------------------
    def writeInputStar(self, starFn, *images):
        """ Easy way to write a simple star file with a single micrographs.
        Used by the relion implementation of motioncor.
        """
        with open(starFn, 'w') as f:
            f.write("data_\nloop_\n_rlnMicrographMovieName\n")
            for img in images:
                f.write("%s\n" % os.path.basename(img.getFileName()))

    def _getMovieOutFn(self, movie, suffix):
        movieBase = pwutils.removeBaseExt(movie.getFileName()).replace('.', '_')
        return os.path.join(self._getOutputMovieFolder(movie), 'output',
                            '%s%s' % (movieBase, suffix))

    def _getAbsPath(self, baseName):
        return os.path.abspath(self._getExtraPath(baseName))

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getPsdCorr(self, movie):
        return self._getNameExt(movie, '_psd_comparison', 'psd', extra=True)

    def _getPsdJpeg(self, movie):
        return self._getNameExt(movie, '_psd', 'jpeg', extra=True)

    def _setPlotInfo(self, movie, mic):
        mic.plotGlobal = em.Image(location=self._getPlotGlobal(movie))
        if self.doComputePSD:
            mic.psdCorr = em.Image(location=self._getPsdCorr(movie))
            mic.psdJpeg = em.Image(location=self._getPsdJpeg(movie))
        if self.doComputeMicThumbnail:
            mic.thumbnail = em.Image(
                location=self._getOutputMicThumbnail(movie))

    def _computeExtra(self, movie):
        """ Compute thumbnail, PSD and plots. """
        inputMovies = self.inputMovies.get()
        movieFolder = self._getOutputMovieFolder(movie)
        outMicFn = self._getMovieOutFn(movie, '.mrc')

        if self.doComputeMicThumbnail:
            self.computeThumbnail(outMicFn,
                                  outputFn=self._getOutputMicThumbnail(movie))

        if self.doComputePSD:
            #fakeShiftsFn = self.writeZeroShifts(movie)
            movieFn = movie.getFileName()
            aveMicFn = os.path.join(movieFolder,
                                    pwutils.removeBaseExt(movieFn) + "_tmp.mrc")
            self.averageMovie(movie, movieFn, aveMicFn,
                              binFactor=self.binFactor.get(),
                              dark=inputMovies.getDark(),
                              gain=inputMovies.getGain())

            self.computePSDs(movie, aveMicFn, outMicFn,
                             outputFnCorrected=self._getPsdJpeg(movie))

        self._saveAlignmentPlots(movie)

    def _moveFiles(self, movie):
        # It really annoying that Relion default names changes if you use DW or not
        # if use DW, the default name are DW and the others noDW
        if self.doDW:
            pwutils.moveFile(self._getMovieOutFn(movie, '.mrc'),
                             self._getExtraPath(self._getOutputMicWtName(movie)))
            if self.saveNonDW:
                pwutils.moveFile(self._getMovieOutFn(movie, '_noDW.mrc'),
                                 self._getExtraPath(self._getOutputMicName(movie)))
        else:
            pwutils.moveFile(self._getMovieOutFn(movie, '.mrc'),
                             self._getExtraPath(self._getOutputMicName(movie)))

    def _getRange(self, movie):
        n = self._getNumberOfFrames(movie)
        iniFrame, _, indxFrame = movie.getFramesRange()
        first, last = self._getFrameRange(n, 'sum')

        if iniFrame != indxFrame:
            first -= iniFrame
            last -= iniFrame

        return first, last

    def _getNumberOfFrames(self, movie):
        _, lstFrame, _ = movie.getFramesRange()

        if movie.hasAlignment():
            _, lastFrmAligned = movie.getAlignment().getRange()
            if lastFrmAligned != lstFrame:
                return lastFrmAligned
        return movie.getNumberOfFrames()

    def _saveAlignmentPlots(self, movie):
        # Create plots and save as an image
        shiftsX, shiftsY = self._getMovieShifts(movie)
        first, _ = self._getFrameRange(movie.getNumberOfFrames(), 'align')
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first)
        plotter.savefig(self._getPlotGlobal(movie))
        plotter.close()


def createGlobalAlignmentPlot(meanX, meanY, first):
    """ Create a plotter with the shift per frame. """
    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax = figure.add_subplot(111)
    ax.grid()
    ax.set_title('Global shift')
    ax.set_xlabel('Shift x (pixels)')
    ax.set_ylabel('Shift y (pixels)')

    i = first
    skipLabels = ceil(len(meanX)/10.0)
    labelTick = 1

    for x, y in izip(meanX, meanY):
        if labelTick == 1:
            ax.text(x - 0.02, y + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    ax.plot(meanX, meanY, color='b')
    ax.plot(meanX, meanY, 'yo')

    plotter.tightLayout()

    return plotter
