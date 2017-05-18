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

from pyworkflow.object import String, Pointer
from pyworkflow.utils.properties import Message
import pyworkflow.utils as pwutils
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.em.protocol import ProtExtractParticles
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import Coordinate, Filament, SetOfCoordinates, SetOfFilaments
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking

import eman2
from pyworkflow.em.packages.eman2.convert import readFilaments
from convert import readSetOfCoordinates
import math


class ProtSegmentHelices(ProtParticlePicking):
    """ Segments set of filaments into set of coordinates """
    _label = 'segment filaments'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        #what else should go here?

 #--------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputFilaments', params.PointerParam,
                      pointerClass='SetOfFilaments',
                      #important=True,
                      label="Input filaments",
                      help='Select the SetOfFilaments from which you want to extract coordinates ')

        form.addParam('overlap', params.IntParam,
                      label='Overlap (px)',                                  ##############change px to A!
                      validators=[params.Positive],
                      help='Overlap of two adjacent particles within a filament.'
                           'This defines how finely you want to segment the filament.')

        #inputFilaments = self.inputFilaments.get()

        form.addParam('boxSize', params.IntParam,
                      #default=self.inputFilaments.get().getBoxSize(),      how to add the boxsize used for picking as default?
                      label='Particle box size (px)',
                      validators=[params.Positive],
                      help='Size of the boxed particles (in pixels). '
                     )

    #--------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        #self._insertFunctionStep('convertInputStep')

        #self._insertFunctionStep('segmentFilamentsStep')

        self._insertFunctionStep('createOutputStep')


    #--------------- STEPS functions -----------------------

    #def convertInputStep(self):
    #    pass

    #def segmentFilamentsStep(self, params):
        #inputFilaments = self.inputFilaments.get()

    def createOutputStep(self):
        inputFilaments = self.inputFilaments.get()
        overlap = self.overlap.get()
        boxsize = self.boxSize.get()
        mics = inputFilaments.getMicrographs()
        outputCoords = self._createSetOfCoordinates(mics)

        micDic = {}
        for mic in mics.iterItems():
            micDic[mic.getObjId()] = mic.clone()

        for filament in inputFilaments.iterFilaments():
            startX = float(filament.getEndpoints()[0])
            startY = float(filament.getEndpoints()[1])
            endX = float(filament.getEndpoints()[2])
            endY = float(filament.getEndpoints()[3])
            length = filament.getLength()
            moveby = boxsize - overlap
            angle = filament.getAngle()
            moveX = math.cos(angle)*moveby
            moveY = math.sin(angle)*moveby
            amountSegments = int(length/moveby)
            mic = micDic[filament.getMicId()]
            print filament.getMicId(), micdic[filament.getMicId()]
            initialCoord=self.makeFilCoord(startX, startY, mic)
            initialCoord.filamentId = filament.getObjId() ##mark
            outputCoords.append(initialCoord)

            for counter in range(amountSegments):
                startX += moveX
                startY += moveY #newCoord.shiftY(moveY) doesn't work - maybe because of float vs int....
                newCoord = self.makeFilCoord(int(startX), int(startY), mic)
                newCoord.filamentId = filament.getObjId() ##mark
                outputCoords.append(newCoord)

            lastCoord = (self.makeFilCoord(endX, endY, mic))
            lastCoord.filamentId = filament.getObjId()##mark
            outputCoords.append(lastCoord)

        outputCoords.filamentsPointer = Pointer()
        outputCoords.filamentsPointer.set(inputFilaments) ##is this correct?
        outputCoords.setBoxSize(boxsize)
        self._defineOutputs(outputCoordinates=outputCoords)
        self._defineSourceRelation(self.inputFilaments, outputCoords)

    #--------------- INFO functions -------------------------

    def _validate(self):
        return []

    def _citations(self):
        return []

    def _summary(self):
        return []

    def _methods(self):
        return []

    #--------------- UTILS functions -------------------------

    def makeFilCoord(self,x,y,mic):
        newCoord = Coordinate()
        newCoord.setPosition(x,y)
        newCoord.setMicrograph(mic)
        return newCoord