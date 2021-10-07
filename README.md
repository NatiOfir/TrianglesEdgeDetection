TrianglesEdgeDetection
======================

Faint edge detection with shape score based on triangle binary partition of the images.
Traingle implementation of the CVPR 2016 paper "Fast Detection of Curved Edges at Low SNR":
https://www.cv-foundation.org/openaccess/content_cvpr_2016/html/Ofir_Fast_Detection_of_CVPR_2016_paper.html

Usage in Matlab:

>> I  = im2double(imread('img.png'));
>> I = runIm(I);

See demo.m for example.

Important Files:
getPrm.m - Params of the algorithm

Good luck!
