#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 17:29:56 2022

@author: dianlyu
"""

import cv2
import numpy as np
import glob
import os


os.chdir('/Users/dianlyu/Dropbox/Stanford_Matters/data/SELF/CCEP/results/prelim_plots/explore5_3/render5')
files0 = np.array(glob.glob('*.jpg'))
time_seq = [];
for t,f in enumerate(files0):
    wordlist = f.split('_')
    time_seq.append(float(wordlist[2][0:-6]))
    
correct_order = np.argsort(time_seq)
files = files0[correct_order].tolist()

img_array = []
for filename in files:
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)


out = cv2.VideoWriter('OutflowCCEP.mp4', cv2.VideoWriter_fourcc(*'MP4V'), 12, size)
 
for i in range(len(img_array)):
    out.write(img_array[i])
out.release()