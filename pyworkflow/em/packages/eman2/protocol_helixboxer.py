# **************************************************************************
# *
# * Authors:        Olivia [1]....
#                   J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# *  [1] OPIC....
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
from pyworkflow.utils.path import join
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.em.protocol import ProtParticlePicking

import eman2
from pyworkflow.em.packages.eman2.convert import loadJson
from convert import readSetOfCoordinates


class EmanProtHelixBoxer(ProtParticlePicking):
    """ Picks particles in a set of micrographs using eman2 boxer. """
    _label = 'helix boxer'
        
    def __init__(self, **args):     
        ProtParticlePicking.__init__(self, **args)
        # The following attribute is only for testing
        self.importFolder = String(args.get('importFolder', None))
    
    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.inputMics = self.inputMicrographs.get()
        micList = [os.path.relpath(mic.getFileName(), self.workingDir.get()) for mic in self.inputMics]

        self._params = {'inputMics': ' '.join(micList)}
        # Launch Boxing GUI
        self._insertFunctionStep('launchBoxingGUIStep', interactive=True)

    #--------------------------- STEPS functions ---------------------------------------------------
    def launchBoxingGUIStep(self):
        # Print the eman version, useful to report bugs
        self.runJob(eman2.getEmanProgram('e2version.py'), '')
        # Program to execute and it arguments
        program = eman2.getEmanProgram("e2helixboxer.py")
        arguments = " --gui %(inputMics)s"
        # Run the command with formatted parameters
        self._log.info('Launching: ' + program + ' ' + arguments % self._params)
        try:
            self.runJob(program, arguments % self._params)
        except Exception, ex:
            print "Warning, error: ", ex

        # Open dialog to request confirmation to create output
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
            #self._leaveDir()# going back to project dir
            # TODO: Parse helix boxer outputs and generate the SetOfFilaments for Scipion
            for mic in self.getInputMicrographs():
                micFn = os.path.relpath(mic.getFileName(), self.workingDir.get())
                micCoordsFn = pwutils.replaceBaseExt(micFn, '.coords.txt')
                arguments = '--helix-coords=extra/%s %s' % (micCoordsFn, micFn)
                self.runJob(program, arguments)

        def _createOutput(self):
            # TODO: Parse helix boxer outputs and generate the SetOfFilaments for Scipion
            # e2helixboxer.py --helix-coords=microgname.txt ../000002_ProtImportMicrographs/extra/Feb22_20.11.33_aligned_mic_DW.mrc
            pass

    #--------------------------- INFO functions ---------------------------------------------------
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
    
    #--------------------------- UTILS functions ---------------------------------------------------
    def _runSteps(self, startIndex):
        # Redefine run to change to workingDir path
        # Change to protocol working directory
        self._enterWorkingDir()
        ProtParticlePicking._runSteps(self, startIndex)
    
    def getFiles(self):
        filePaths = self.inputMicrographs.get().getFiles() | ProtParticlePicking.getFiles(self)
        return filePaths

    def readSetOfCoordinates(self, workingDir, coordSet):
        readSetOfCoordinates(workingDir, self.inputMics, coordSet)
        
        
