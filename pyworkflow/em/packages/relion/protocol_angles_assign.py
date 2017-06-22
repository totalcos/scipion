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
from pyworkflow.em.protocol.protocol_particles import ProtProcessParticles
from pyworkflow.protocol.params import BooleanParam
import pyworkflow.em as em
from pyworkflow.em import metadata as md
from convert import rowToAlignment, alignmentToRow
from random import randint
import numpy

class ProtAnglesAssign(ProtProcessParticles):
    """Allows the user to assign a tilt angle of 90 deg and randomize the rotational angle of particles"""
    _label = 'assign angles'

    def __init__(self, **kwargs):
        ProtProcessParticles.__init__(self, **kwargs)

    #--------------- DEFINE param functions ---------------

    def _defineParams(self, form):

        ProtProcessParticles._defineParams(self, form)

    def _defineProcessParams(self,form):
        form.addParam('assignTilt90', BooleanParam, label='Assign out of plane angle (Tilt) to 90 deg?',
                      help='Assigns a value of 90 degrees to the tilt angle of each particle.')

        form.addParam('randomizeRot', BooleanParam,
                      label='Randomize the rotational angle around the particle axis (Rot)?',
                      help='Assigns random values to the rotational angles around the particle axis.')

    #--------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    #--------------------------- STEPS functions ---------------------------------------------------

    def createOutputStep(self):
        particles = self.inputParticles.get()
        partSet = self._createSetOfParticles()
        #partSet.copyInfo(particles.getMicrographs())
        if self.assignTilt90:
            partSet.setAlignmentProj()
            alignRow = md.Row()
            for part in particles.iterItems():

                #retrieving old alignment#
                if particles.hasAlignment2D():
                    row2D = md.Row() #maybe unnecessary?
                    transform = part.getTransform()
                    alignmentToRow(transform, row2D, em.ALIGN_2D)
                    alignRow.copyFromRow(row2D)
                if particles.hasAlignmentProj():
                    rowProj = md.Row()
                    transform = part.getTransform()
                    alignmentToRow(transform, rowProj, em.ALIGN_PROJ)
                    alignRow.copyFromRow(rowProj)

                alignRow.setValue(md.RLN_ORIENT_TILT, 90)
                transform = rowToAlignment(alignRow, em.ALIGN_PROJ)
                part.setTransform(transform)
                partSet.append(part)

        #This part is not ready yet. I am getting an error when trzing to iterate through new set
        #Strange, since I don't when I do it further down...
        #if self.randomizeRot:
         #   partSet2 = self._createSetOfParticles()
          #  partSet2.setAlignmentProj()
           # for part in partSet.iterItems():
            #    alignRow = md.Row()
             #   alignRow.setValue(md.RLN_ORIENT_ROT, randint(-180,180))
              #  transform = rowToAlignment(alignRow, em.ALIGN_PROJ)
               # part.setTransform(transform)
                #partSet2.append(part)

        #Checking if the old SetOfParticles (particles) has the new alignment:
        for part in particles.iterItems():
            row = md.Row()
            trans = part.getTransform()
            alignmentToRow(trans, row, em.ALIGN_PROJ)
            print row
        #Checking if the new SetOfParticles (particles) has the new alignment:
        for part in partSet.iterItems():
            row = md.Row()
            trans = part.getTransform()
            alignmentToRow(trans, row, em.ALIGN_PROJ)
            print row

        #ERROR: outputParticles has no sampling rate!!!
        
        self._defineOutputs(outputParticles=partSet)
        self._defineSourceRelation(self.inputParticles, partSet)
