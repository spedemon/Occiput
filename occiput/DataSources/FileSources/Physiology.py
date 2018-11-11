# -*- coding: utf-8 -*-
# occiput
# Harvard University, Martinos Center for Biomedical Imaging
# Aalto University, Department of Computer Science

try:
    from Occiput_Interface_Biograph_mMR import Biograph_mMR_Physiology
except:
    Biograph_mMR_Physiology = None

try:
    from Occiput_Interface_Brain_PET import Brain_PET_Physiology
except:
    Brain_PET_Physiology = None
