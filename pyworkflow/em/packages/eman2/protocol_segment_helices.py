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

import math

from pyworkflow.object import Pointer, Integer
import pyworkflow.protocol.params as params
from pyworkflow.em.data import Coordinate
from pyworkflow.em.protocol.protocol_particles import ProtParticlePicking


class ProtSegmentHelices(ProtParticlePicking):
    """ Segments set of filaments into set of coordinates """
    _label = 'segment helices'

    def __init__(self, **kwargs):
        ProtParticlePicking.__init__(self, **kwargs)

 #--------------- DEFINE param functions ---------------

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputFilaments', params.PointerParam,
                      pointerClass='SetOfFilaments',
                      important=True,
                      label="Input filaments",
                      help="Select the SetOfFilaments from which you would "
                           "like to extract coordinates. ")

        form.addParam('boxSize', params.IntParam,
                      #default=self.inputFilaments.get().getBoxSize(),      how to add the boxsize used for picking as default?
                      label='Particle box size (px)',
                      validators=[params.Positive],
                      help='Size of the boxed particles (in pixels).')

        form.addParam('overlap', params.FloatParam,
                      label='Overlap (A)',
                      validators=[params.Positive],        ###a nice default would be 90% of default boxsize (cf. relion wiki)
                      help='Overlap of two adjacent particles within a '
                           'filament (in Angstrom).'
                           '(Measured along the filament line.)'
                           'This defines how finely you want to segment the '
                           'filament.')

    #--------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    #--------------- STEPS functions -----------------------

    def createOutputStep(self):
        inputFilaments = self.inputFilaments.get()
        mics = inputFilaments.getMicrographs()
        pixSize = mics.getScannedPixelSize()
        # transforms the overlap from A to pixels
        overlap = self.overlap.get() / pixSize
        boxsize = self.boxSize.get()
        outputCoords = self._createSetOfCoordinates(mics)

        micDict = {}
        for mic in mics.iterItems():
            micDict[mic.getObjId()] = mic.clone()

        for filament in inputFilaments.iterFilaments():
            startX, startY, endX, endY = filament.getEndpoints()
            length = filament.getLength()
            moveby = boxsize - overlap
            angle = filament.getAngle()*math.pi/180
            moveX = math.cos(angle)*moveby
            moveY = math.sin(angle)*moveby
            amountSegments = int(length/moveby)
            mic = micDict[filament.getMicId()]
            coord = self.makeFilCoord(startX, startY, mic, filament.getObjId())
            outputCoords.append(coord)

            for counter in range(amountSegments):
                # We can re-use the same Coordinate object but updating
                # the x and y values and clean
                startX += moveX
                startY += moveY
                coord.setX(round(startX))
                coord.setY(round(startY))
                coord.setObjId(None) # Reset id to insert as a new coordinate
                outputCoords.append(coord)

            lastCoord = self.makeFilCoord(endX, endY, mic, filament.getObjId())
            outputCoords.append(lastCoord)

        outputCoords.filamentsPointer = Pointer()
        outputCoords.filamentsPointer.set(inputFilaments)
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

    def makeFilCoord(self, x, y, mic, filamentId):
        newCoord = Coordinate()
        newCoord.setPosition(x,y)
        newCoord.setMicrograph(mic)
        newCoord.filamentId = Integer(filamentId)
        return newCoord
