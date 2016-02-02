# **************************************************************************
# *
# * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
# *  e-mail address 'xmipp@cnb.csic.es'
# *
# **************************************************************************

import pyworkflow.object as pwobj
from pyworkflow.em import *  
from xmipp import MetaData, MDL_ANGLE_ROT, MDL_SHIFT_Z
from pyworkflow.em.packages.xmipp3.convert import getImageLocation
from pyworkflow.protocol.constants import LEVEL_ADVANCED

from helix import HelicalFinder



class XmippProtHelicalSymmetrize(ProtPreprocessVolumes):
    """ Estimate helical parameters and symmetrize.
    
        Helical symmetry is defined as:
        V(r,rot,z)=V(r,rot+k*DeltaRot,z+k*Deltaz).
    """
    _label = 'helical symmetrize'

    #--------------------------- DEFINE param functions --------------------------------------------
    def _defineParams(self, form):
        form.addSection(label='General parameters')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input volume')
        form.addParam('cylinderInnerRadius', IntParam, default=-1,
                      label='Cylinder inner radius',
                      help="The helix is supposed to occupy this radius in "
                           "voxels around the Z axis. Leave it as -1 for "
                           "symmetrizing the whole volume")
        form.addParam('cylinderOuterRadius',IntParam,
                      label='Cylinder outer radius', default=-1,
                      help="The helix is supposed to occupy this radius in voxels "
                           "around the Z axis. Leave it as -1 for symmetrizing "
                           "the whole volume")
        form.addParam('dihedral',BooleanParam, default=False,
                      label='Apply dihedral symmetry')
        form.addParam('forceDihedralX',BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label='Force the dihedral axis to be in X',
                      help="If this option is chosen, then the dihedral axis is "
                           "not searched and it is assumed that it is around X.")

        form.addSection(label='Search limits')
        form.addParam('heightFraction',FloatParam,default=0.9,label='Height fraction',
                      help="The helical parameters are only sought using the fraction indicated by this number. " \
                           "In this way, you can avoid including planes that are poorly resolved at the extremes of the volume. " \
                           "However, note that the algorithm can perfectly work with a fraction of 1.")
        form.addParam('rot0',FloatParam,default=0,label='Minimum rotational angle',help="In degrees")
        form.addParam('rotF',FloatParam,default=360,label='Maximum rotational angle',help="In degrees")
        form.addParam('rotStep',FloatParam,default=5,label='Angular step',help="In degrees")
        form.addParam('z0',FloatParam,default=0,label='Minimum shift Z',help="In Angstroms")
        form.addParam('zF',FloatParam,default=10,label='Maximum shift Z',help="In Angstroms")
        form.addParam('zStep',FloatParam,default=0.5,label='Shift step',help="In Angstroms")
        self.deltaZ=Float()
        self.deltaRot=Float()
        form.addParallelSection(threads=4, mpi=0)

    #--------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('coarseSearchStep')
        self._insertFunctionStep('fineSearchStep')
        self._insertFunctionStep('symmetrizeStep')
        self._insertFunctionStep('createOutputStep')

        # Store some general attributes
        inputVol = self.inputVolume.get()
        self.fnVol = getImageLocation(inputVol)
        self.fnVolSym = self._getPath('volume_symmetrized.vol')
        self.height = inputVol.getXDim()
        mask = (self.cylinderInnerRadius.get(),
                self.cylinderOuterRadius.get(),
                self.height)
        self.helical = HelicalFinder(self, inputVol.getSampling(),
                                     self.heightFactor.get(),
                                     self.dihedral.get(), mask)

    #--------------------------- STEPS functions --------------------------------------------
    def convertInputStep(self):
        if self.dihedral:
            if not self.forceDihedralX:
                self.runJob("xmipp_transform_symmetrize",
                            "-i %s -o %s --sym dihedral --dont_wrap" % (self.fnVol, self.fnVolSym))
            else:
                self.runJob("xmipp_transform_geometry",
                            "-i %s -o %s --rotate_volume axis 180 1 0 0" % (self.fnVol, self.fnVolSym))
                self.runJob("xmipp_image_operate",
                            "-i %s --plus %s -o %s" % (self.fnVol, self.fnVolSym, self.fnVolSym))
                self.runJob("xmipp_image_operate",
                            "-i %s --mult 0.5" % self.fnVolSym)
        else:
            ImageHandler().convert(self.inputVolume.get(), self.fnVolSym)

    def coarseSearchStep(self):
        z = (self.z0.get(), self.zF.get(), self.zStep.get())
        rot = (self.rot0.get(), self.rotF.get(), self.rotStep.get())
        self.helical.searchCoarse(self.fnVolSym,
                                  self._getExtraPath('coarseParams.xmd'),
                                  z, rot, self.numberOfThreads.get())

    def fineSearchStep(self):
        z = (self.z0.get(), self.zF.get())
        rot = (self.rot0.get(), self.rotF.get())
        self.helical.searchFine(self.fnVolSym,
                                self._getExtraPath('coarseParams.xmd'),
                                self._getExtraPath('fineParams.xmd'),
                                z, rot)

    def symmetrizeStep(self):
        self.helical.symmetrize(self.fnVolSym,
                           self._getExtraPath('fineParams.xmd'), self.fnVolSym)\

    def createOutputStep(self):
        volume = Volume()
        volume.setFileName(self._getPath('volume_symmetrized.vol'))
        #volume.setSamplingRate(self.inputVolume.get().getSamplingRate())
        volume.copyInfo(self.inputVolume.get())
        self._defineOutputs(outputVolume=volume)
        self._defineTransformRelation(self.inputVolume, self.outputVolume)

        md = MetaData(self._getExtraPath('fineParams.xmd'))
        objId = md.firstObject()
        self._defineOutputs(deltaRot=pwobj.Float(md.getValue(MDL_ANGLE_ROT, objId)),
                            deltaZ=pwobj.Float(md.getValue(MDL_SHIFT_Z, objId)))

        #--------------------------- INFO functions --------------------------------------------
    def _summary(self):
        messages = []
        if self.deltaZ.hasValue():
            messages.append('DeltaZ=%f (voxels) %f (Angstroms)'%(self.deltaZ.get()/self.inputVolume.get().getSamplingRate(),self.deltaZ.get()))
            messages.append('DeltaRot=%f (degrees)'%self.deltaRot.get())
        return messages

    def _citations(self):
        papers=[]
        return papers

    def _methods(self):
        messages = []
        messages.append('We looked for the helical symmetry parameters of the volume %s using Xmipp [delaRosaTrevin2013].' % self.getObjectTag('inputVolume'))
        if self.deltaZ.hasValue():
            messages.append('We found them to be %f Angstroms and %f degrees.'%(self.deltaZ.get(),self.deltaRot.get()))
            messages.append('We symmetrized %s with these parameters and produced the volume %s.'%(self.getObjectTag('inputVolume'),
                                                                                                   self.getObjectTag('outputVolume')))
            if self.dihedral.get():
                messages.append('We applied dihedral symmetry.')
        return messages

