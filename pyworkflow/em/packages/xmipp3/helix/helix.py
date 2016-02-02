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

import pyworkflow.em.metadata as md


class HelicalFinder():
    """ Helper class to find helical symmetry and apply it. """

    def __init__(self, protocol, sampling, heightFraction, dihedral,
                 mask=(-1, -1, -1)):
        """
        Params:
            protocol: protocol object that will use the HelicalFinder.
            sampling: sampling rate (or pixel size)
            heightFraction: ???
            dihedral: boolean, True if dihedral symmetry.
            mask: 3-tuple: outerRadius, innerRadius and height
                (use -1 for no mask cilinder/tube)
        """
        self.protocol = protocol
        self.runJob = protocol.runJob # borrow runJob method from protocol
        self.sampling = sampling
        self.heightFraction = heightFraction
        self.dihedral = dihedral
        self.mask = mask

    def getSymmetry(self):
        return 'helical%s' % 'Dihedral' if self.dihedral else ''

    def hasMask(self):
        """ We assume to have mask if the outerRadius is > 0. """
        return self.mask[0] > 0

    def searchCoarse(self, fnVol, fnOut,
                     z0, zF, zStep,
                     rot0, rotF, rotStep,
                     numberOfThreads):
        args = self._commonArgs(fnVol, fnOut)
        args += " -z %f %f %f " % (z0, zF, zStep)
        args += " --rotHelical %f %f %f" % (rot0, rotF, rotStep)
        args += " --thr %s" % numberOfThreads
        self.runJob('xmipp_volume_find_symmetry', args)

    def searchFine(self, fnVol, fnCoarse, fnFine, z0, zF, rot0, rotF):
        rotInit, zInit = self._paramsFromMd(fnCoarse)
        args = self._commonArgs(fnVol, fnFine)
        args += " --localHelical %f %f " % (zInit, rotInit)
        args += " -z %f %f 1 --rotHelical %f %f 1" % (z0, zF, rot0, rotF)
        self.runJob('xmipp_volume_find_symmetry',args)

    def symmetrize(self, fnVol, fnParams, fnOut):
        rot0, z0 = self._paramsFromMd(fnParams)
        args = self._commonArgs(fnVol, fnOut, maskArgs=False)
        args += " --helixParams %f %f -o %s" % (z0, rot0)
        self.runJob('xmipp_transform_symmetrize', args)

        if self.hasMask():
            args = "-i %s %s" % (fnVol, self._maskArgs())
            self.runJob('xmipp_transform_mask',args)

    #--------- Internal functions -----------------
    def _maskArgs(self):
        arg = ''
        if self.hasMask():
            i, o, h = self.mask
            if self.mask[1] < 0:
                arg = ' --mask cylinder -%d -%d' % (o, h)
            else:
                arg = ' --mask tube -%d -%d -%d' % (i, o, h)
        return arg

    def _commonArgs(self, inFn, outFn, maskArgs=True):
        args =  ' -i %s -o %s --sym %s' % (inFn, outFn, self.getSymmetry())
        args += ' --heightFraction %f' % self.heightFraction
        args += ' --sampling %f' % self.sampling
        if maskArgs:
            args += self._maskArgs()

    def _paramsFromMd(self, fnMd):
        row = md.getFirstRow(fnMd)
        return row.getValue(md.MDL_ANGLE_ROT), row.getValue(md.MDL_SHIFT_Z)
