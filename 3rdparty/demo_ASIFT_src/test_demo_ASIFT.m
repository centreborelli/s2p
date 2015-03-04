% *************************************************************************
% NOTE: The ASIFT SOFTWARE ./demo_ASIFT IS STANDALONE AND CAN BE EXECUTED
%       WITHOUT MATLAB. 
% *************************************************************************
%
% demo_ASIFT.m is a MATLAB interface for the ASIFT software. This
% script provides an example of using demo_ASIFT.m. The input/output of
% demo_ASIFT.m follow those of the ASIFT software and a description can
% found in README.txt. (Note that the ASIFT software support only the PNG format, 
% the Matlab interface reads most standard image formats.)
%
% The ASIFT C/C++ source code must be compiled before running the ASIFT software
% for the first time. See README.txt for more details. 
% 
% Copyright, Jean-Michel Morel, Guoshen Yu, 2008. 
%
% Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
%
% Reference: J.M. Morel and G.Yu, ASIFT: A New Framework for Fully Affine Invariant Image 
%           Comparison, SIAM Journal on Imaging Sciences, vol. 2, issue 2, pp. 438-469, 2009. 
% Reference: ASIFT online demo (You can try ASIFT with your own images online.) 
% http://www.ipol.im/pub/algo/my_affine_sift/
%
% 2010.08.17

file_img1 = 'adam1.png';
file_img2 = 'adam2.png';
imgOutVert = 'imgOutVert.png';
imgOutHori = 'imgOutHori.png';
matchings = 'matchings.txt';
keys1 = 'keys1.txt';
keys2 = 'keys2.txt';
flag_resize = 0;

demo_ASIFT(file_img1, file_img2, imgOutVert, imgOutHori, matchings, keys1, keys2, flag_resize);


