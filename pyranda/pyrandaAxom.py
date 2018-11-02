################################################################################
# Copyright (c) 2018, Lawrence Livemore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# LLNL-CODE-749864
# This file is part of pyranda
# For details about use and distribution, please read: pyranda/LICENSE
#
# Written by: Britton J. Olson, olson45@llnl.gov
################################################################################
import numpy 
import re
import sys
import time
from .pyrandaPackage import pyrandaPackage

from . import parcop

class pyrandaAxom(pyrandaPackage):

    def __init__(self,pysim):

        PackageName = 'Axom'
        pyrandaPackage.__init__(self,PackageName,pysim)

        

    def get_sMap(self):
        sMap = {}
        sMap['axom.stlfile('] =  "self.packages['Axom'].questInit("
        sMap['axom.levelset(']  =  "self.packages['Axom'].quest("
        self.sMap = sMap





    def quest_init(self,filename):

        parcop.parcop.init_quest( filename )


    def quest_get(self,x,y,z):

        
        distance = parcop.parcop.get_quest(x,y,z)


        return distance
