# **************************************************************************
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
# **************************************************************************

from pyworkflow.em import *
from pyworkflow.em.packages.eman2 import ProtSegmentHelices
from pyworkflow.em.packages.relion import ProtRelionExtractParticles

from pyworkflow.tests import *
from test_workflow import TestWorkflow



class TestWorkflowHelical(TestWorkflow):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet('tmv_helix')

    def test_tmv(self):
        protImport = self.newProtocol(ProtImportMicrographs,
                                      objLabel='import TMV 4 mics',
                                      filesPath=self.ds.getFile('micrographs'),
                                      filesPattern='TMV*.mrc',
                                      samplingRateMode=0,
                                      samplingRate=1.062,
                                      voltage=300,
                                      sphericalAberration=2.7)
        self.launchProtocol(protImport)
        self.assertIsNotNone(protImport.outputMicrographs,
                             "There was a problem with the import")

        protImportFil = self.newProtocol(ProtImportFilaments,
                                         filesPath=self.ds.getFile('coords'),
                                         filesPattern='TMV*.box',
                                         boxSize=250)
        protImportFil.inputMicrographs.set(protImport.outputMicrographs)
        self.launchProtocol(protImportFil)

        protSegment = self.newProtocol(ProtSegmentHelices,
                                       boxSize=250,
                                       overlap=600)
        protSegment.inputFilaments.set(protImportFil.outputFilaments)
        self.launchProtocol(protSegment)

        protExtract = self.newProtocol(ProtRelionExtractParticles,
                                       boxSize=250,
                                       doInvert=True,
                                       doNormalize=True)
        protExtract.inputCoordinates.set(protSegment.outputCoordinates)
        self.launchProtocol(protExtract)



