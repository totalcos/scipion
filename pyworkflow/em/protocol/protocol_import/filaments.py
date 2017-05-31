# **************************************************************************
# *
# * Authors:        Olivia Pfeil-Gardiner (zolivia@zedat.fu-berlin.de) [1]
#                   J.M. De la Rosa Trevin (delarosatrevin@scilifelab.se) [2]
# *
# *  [1] OPIC, Division of Structural Biology, University of Oxford
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

from base import ProtImportFiles
from pyworkflow.protocol.params import IntParam, PointerParam, FloatParam, BooleanParam
from pyworkflow.utils.path import removeBaseExt

class ProtImportFilaments(ProtImportFiles):
    """ Protocol to import a set of filaments """
    _label = 'import filaments'

    IMPORT_FROM_EMAN = 0

    def _getImportChoices(self):
        """ Return a list of possible choices
        from which the import can be done.
        (usually packages formats such as: xmipp3, eman2, relion...etc.
        """
        return ['eman']

    #--------------------------- DEFINE param functions --------------------------------------------

    def _defineParams(self, form):

        ProtImportFiles._defineParams(self, form)

    def _defineImportParams(self, form):
        form.addParam('inputMicrographs', PointerParam, pointerClass='SetOfMicrographs',
                          label='Input micrographs',
                          help='Select the micrographs for which you want to import filaments.')

        form.addParam('boxSize', IntParam,
                      label='Box size')

        #form.addParam('scale', FloatParam,
        #              label='Scale', default=1,
        #              help='Factor to scale coordinates')

     #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):

        importFrom = self.importFrom.get()
        self._insertFunctionStep('createOutputStep', importFrom,
                                     self.filesPath.get())

    #--------------------------- STEPS functions ---------------------------------------------------

    def createOutputStep(self, importFrom, *args):

        inputMics = self.inputMicrographs.get()
        filamentsSet = self._createSetOfFilaments(inputMics)
        filamentsSet.setBoxSize(self.boxSize.get())
        ci = self.getImportClass()

        for filamentFile, fileId in self.iterFiles():
            mic = self.getMatchingMic(filamentFile,fileId)
            if mic is not None:
                def addFilament(fil):
                    fil.setMicrograph(mic)
                    filamentsSet.append(fil)
                # Parse the filaments in the given format for this micrograph
                ci.importFilaments(mic, filamentFile, filamentsSet)


        self._defineOutputs(outputFilaments=filamentsSet)
        self._defineSourceRelation(self.inputMicrographs, filamentsSet)

  #--------------------------- UTILS functions ---------------------------------------------------
    def getImportClass(self):
        """ Return the class in charge of importing the files. """

        #importFrom = self.ImportFrom.get()

        #if importFrom == self.IMPORT_FROM_EMAN:
        from pyworkflow.em.packages.eman2.dataimport import EmanImport
        return EmanImport(self)

        #elif importFrom == ....

        #else:
        #    self.importFilePath = ''
        #    return None

    def getMatchingMic(self, filamentFile, fileId):
        """ Given a filaments file check if there is a micrograph
        that matches this file. 
        """
        micSet = self.inputMicrographs.get()

        if fileId is None:
            coordBase = removeBaseExt(filamentFile)
            for mic in micSet:
                micBase = removeBaseExt(mic.getFileName())
                if coordBase in micBase or micBase in coordBase: #temporal use of in
                    return mic
            return None
        else:
            return micSet[fileId]



















