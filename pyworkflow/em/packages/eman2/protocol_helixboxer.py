# **************************************************************************
# *
# * Authors:        Olivia Pfeil-Gardiner (zolivia@zedat.fu-berlin.de) [1]
#                   J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# *  [1] Oxford Particle Imaging Centre, Division of Structural Biology, University of Oxford
# *  [2] SciLifeLab, Stockholm University
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

import os

from pyworkflow.object import String
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.em.protocol import ProtParticlePicking
from pyworkflow.protocol.params import IntParam

import eman2
from pyworkflow.em.packages.eman2.convert import readFilaments
from convert import readSetOfCoordinates


class EmanProtHelixBoxer(ProtParticlePicking):
    """ Picks particles in a set of micrographs using eman2 boxer. """
    _label = 'helix boxer'
        
    def __init__(self, **kwargs):
        ProtParticlePicking.__init__(self, **kwargs)
        self.program = eman2.getEmanProgram("e2helixboxer.py")

    def _defineParams(self, form):
        ProtParticlePicking._defineParams(self, form)
        form.addParam('boxSize', IntParam, default=100,
                   label='Box size (px)')

    #--------------------------- INSERT steps functions ------------------------
    def _insertAllSteps(self):
        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUIStep', interactive=True)

    #--------------------------- STEPS functions -------------------------------
    def launchBoxingGUIStep(self):
        inputMics = self.getInputMicrographs()

        # Print the eman version, useful to report bugs
        self.runJob(eman2.getEmanProgram('e2version.py'), '')

        # Prepare the list of micrographs for the command line
        micList = [self.getRelPath(mic.getFileName()) for mic in inputMics]
        arguments = "--helix-width %d" % self.boxSize
        arguments += " --gui %s" % ' '.join(micList)

        self.runHelixBoxer(arguments)

        # Open dialog to request confirmation to create output
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
            self._createOutput(inputMics)

    def _createOutput(self, inputMics):
        self.info("Parsing .box files and creating the output SetOfFilaments.")

        filamentSet = self._createSetOfFilaments(inputMics)
        filamentSet.setBoxSize(self.boxSize.get())

        for mic in self.getInputMicrographs():
            micFn = self.getRelPath(mic.getFileName())
            baseName = pwutils.removeBaseExt(micFn)
            micCoordsFn = self._getExtraPath(baseName + '_coords.box')
            micCoordsFnRel = self.getRelPath(micCoordsFn)
            arguments = '--helix-coords=%s %s' % (micCoordsFnRel, micFn)
            result = self.runHelixBoxer(arguments)
            # If result=False it could means that there was not filaments
            # picked for this micrograph, so we just ignore it.
            # In the future we might implement a cleaner way to do this
            if result:
                readFilaments(mic, micCoordsFn, filamentSet)

        # Define that the resulting SetOfFilaments as output and that is was
        # produced from the input SetOfMicrographs
        self._defineOutputs(outputFilaments=filamentSet)
        self._defineSourceRelation(self.getInputMicrographsPointer(),
                                   filamentSet)


    #--------------------------- INFO functions --------------------------------
    def _validate(self):
        errors = []
        eman2.validateVersion(self, errors)
        return errors

    def _warnings(self):
        warnings = []
        firstMic = self.inputMicrographs.get().getFirstItem()
        fnLower = firstMic.getFileName().lower()
        
        if '.tif' in fnLower:
            warnings.append('Micrographs in .tif format are displayed upside down in EMAN2!!!')
            warnings.append('The generated coordinates will not be valid in Scipion.')
            warnings.append('TIP: a workaround could be importing the coordinates inverted')
            warnings.append('     in Xmipp particle picking GUI.')
        return warnings
    
    #--------------------------- UTILS functions -------------------------------

    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths

    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMics, coordSet)

    def getRelPath(self, filename):
        """ Get the relative micrograph filename from the protocol working dir.
        """
        return os.path.relpath(filename, self.getWorkingDir())

    def runHelixBoxer(self, arguments):
        """ Run e2helixboxer.py EMAN2 program, return True if succeed,
        False otherwise.
        """
        try:
            # Run the command with formatted parameters
            self.runJob(self.program, arguments, cwd=self.getWorkingDir())
            result = True
        except Exception, ex:
            result = False

        return result

