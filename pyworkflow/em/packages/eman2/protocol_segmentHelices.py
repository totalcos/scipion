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
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.em.protocol import ProtExtractParticles
import pyworkflow.protocol.params as params
from pyworkflow.em.protocol import EMProtocol
from pyworkflow.em.data import Coordinate, Filament, SetOfCoordinates, SetOfFilaments

import eman2
from pyworkflow.em.packages.eman2.convert import readFilaments
from convert import readSetOfCoordinates
import math


class ProtSegmentHelices(EMProtocol):
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
                      label='overlap',
                      validators=[params.Positive],
                      help='Overlap of two adjacent particles within a filament.'
                           'This defines how finely you want to segment the filament.')

        form.addParam('boxSize', params.IntParam, default=inputFilaments.getBoxSize()
                     label='Particle box size (px)',
                      validators=[params.Positive],
                      help='Size of the boxed particles (in pixels). '
                     )

    #--------------- INSERT steps functions ----------------

    def _insertAllSteps(self):
        #self._insertFunctionStep('convertInputStep')

        #self.insertFunctionStep('segmentFilamentsStep')

        self.insertFunctionStep('createOutputStep')


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

        for filament in inputFilaments.iterFilaments():
            startX = float(filament.getEndpoints()[0])  #define some params to calculate coordinates
            startY = float(filament.getEndpoints()[1])
            endX = float(filament.getEndpoints()[2])
            endY = float(filament.getEndpoints()[3])
            filLength = math.sqrt((startX-endX)**2 + (startY-endY)**2) ######make a function 'filament.getLength'
            moveby = boxsize - overlap
            angle = math.tan(abs(startY-endY)/abs(startX-endX))
            moveX = math.cos(angle)*moveby
            moveY = math.sin(angle)*moveby
            amountSegments = filLength/
            newCoord=Coordinate(x=startX, y=startY)####HoToDefineACoordinate? --> Maybe with a function setFilCoord(x,y,filament)
            outputCoords.append(newCoord) ######HowToAddToSetOfCoordinates?
            for counter in range(amountSegments):
                x+=
                y=+
                outputCoords.append()
            outputCoords.append(endX, endY)



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
