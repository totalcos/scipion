# **************************************************************************
# *
# * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'jmdelarosa@cnb.csic.es'
# *
# **************************************************************************

import pyworkflow.tests as pwtests
import pyworkflow.em as em

from pyworkflow.em.packages.xmipp3 import (XmippProtHelicalSymmetrize,
                                           XmippProtHelixRefine)
from test_protocols_xmipp_3d import TestXmippBase



class TestXmippProtHelicalSymmetrize(TestXmippBase):
    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.ds = pwtests.DataSet.getDataSet('general')
        cls.vol = cls.ds.getFile('vol_helix')
        cls.protImport = cls.runImportVolumes(cls.vol, 1.0)

    def testHelicalParameters(self):
        print "Run symmetrize helical"
        protHelical = self.newProtocol(XmippProtHelicalSymmetrize,
                                       cylinderOuterRadius=20,
                                       dihedral=True,
                                       rot0=50,rotF=70,rotStep=5,
                                       z0=5,zF=10,zStep=0.5)

        protHelical.inputVolume.set(self.protImport.outputVolume)
        self.proj.launchProtocol(protHelical, wait=True)

        self.assertIsNotNone(protHelical.outputVolume,
                             "There was a problem with Helical output volume")
        self.assertIsNotNone(protHelical.deltaRot.get(),
                             "Output delta rot is None")
        self.assertIsNotNone(protHelical.deltaZ.get(),
                             "Output delta Z is None")
        self.assertAlmostEqual(protHelical.deltaRot.get(), 59.59, delta=1,
                               msg="Output delta rot is wrong")
        self.assertAlmostEqual(protHelical.deltaZ.get(), 6.628, delta=0.2,
                               msg="Output delta Z is wrong")


class TestXmippProtHelicalProjMatching(TestXmippBase):

    @classmethod
    def setUpClass(cls):
        pwtests.setupTestProject(cls)
        cls.dataset = pwtests.DataSet.getDataSet('relion_tutorial')
        cls.vol = cls.dataset.getFile('volume')

    def test_case1(self):
        print "Import Particles"
        protImportParts = self.newProtocol(em.ProtImportParticles,
                                 objLabel='Particles from scipion',
                                 importFrom=em.ProtImportParticles.IMPORT_FROM_SCIPION,
                                 sqliteFile=self.dataset.getFile('import/case2/particles.sqlite'),
                                 magnification=50000,
                                 samplingRate=7.08,
                                 haveDataBeenPhaseFlipped=True
                                 )
        self.launchProtocol(protImportParts)
        self.assertIsNotNone(protImportParts.getFiles(), "There was a problem with the import")
        
        print "Get a Subset of particles"
        protSubset = self.newProtocol(em.ProtSubSet,
                                         objLabel='100 Particles',
                                         chooseAtRandom=True,
                                         nElements=100)
        protSubset.inputFullSet.set(protImportParts.outputParticles)
        self.launchProtocol(protSubset)
        
        print "Import Volume"
        protImportVol = self.newProtocol(em.ProtImportVolumes,
                                         objLabel='Volume',
                                         filesPath=self.vol,
                                         samplingRate=7.08)
        self.launchProtocol(protImportVol)
        self.assertIsNotNone(protImportVol.getFiles(), "There was a problem with the import")
        
        print "Run Projection Matching"
        protProjMatch = self.newProtocol(XmippProtHelixRefine,
                                         ctfGroupMaxDiff=0.00001,
                                         mpiJobSize=10,
                                         numberOfIterations=2,
                                         numberOfThreads=2,
                                         numberOfMpi=3)
        protProjMatch.inputParticles.set(protSubset.outputParticles)
        protProjMatch.input3DReferences.set(protImportVol.outputVolume)
        self.launchProtocol(protProjMatch)
        self.assertIsNotNone(protProjMatch.outputVolume, "There was a problem with Projection Matching")


