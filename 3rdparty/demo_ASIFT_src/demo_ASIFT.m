%*-------------------demo_ASIFT MATLAB interface   -------------------------*/
% 
% *************************************************************************
% NOTE: The ASIFT SOFTWARE ./demo_ASIFT IS STANDALONE AND CAN BE EXECUTED
%       WITHOUT MATLAB. 
% *************************************************************************
% 
% Detect corresponding points in two images with the ASIFT method. 
% Copyright, Jean-Michel Morel, Guoshen Yu, 2008. 
% Please report bugs and/or send comments to Guoshen Yu yu@cmap.polytechnique.fr
%
% Reference: J.M. Morel and G.Yu, ASIFT: A New Framework for Fully Affine Invariant Image 
%           Comparison, SIAM Journal on Imaging Sciences, vol. 2, issue 2, pp. 438-469, 2009. 
% Reference: ASIFT online demo (You can try ASIFT with your own images online.) 
% http://www.ipol.im/pub/algo/my_affine_sift/
%
% 2010.08.17
% ---------------------------------------------------------------------------*/

function demo_ASIFT(file_img1, file_img2, imgOutVert, imgOutHori, matchings, keys1, keys2, flag_resize)

if (nargin == 8) 
    if (isnumeric(flag_resize) && (flag_resize == 0))
        flag_resize = 0;
    else
        flag_resize = 1;
    end
elseif (nargin == 7)
    flag_resize = 1;
else
    disp('*******************************************************************************');
    disp('***************************  ASIFT image matching  ****************************');
    disp('*******************************************************************************');
    disp('Usage: ./demo_ASIFT imgIn1 imgIn2 imgOutVert.png imgOutHori.png');
    disp('          matchings.txt keys1.txt keys2.txt [Resize option: 0/1]'); 
    disp('- imgIn1, imgIn2: input images (in most standard formats).');
    disp('- imgOutVert.png, imgOutHori.png: output images (vertical/horizontal concatenated,');
    disp('  in PNG format.) The detected matchings are connected by write lines.');
    disp('- matchings.txt: coordinates of matched points (col1, row1, col2, row2).'); 
    disp('- keys1.txt keys2.txt: ASIFT keypoints of the two images.')
    disp('- [optional 0/1]. 1: input images resize to 800x600 (default). 0: no resize.'); 
    disp('*******************************************************************************'); 
    disp('*********************  Jean-Michel Morel, Guoshen Yu, 2010 ********************'); 
    disp('*******************************************************************************'); 
    assert(false);
end

imgIn1 = imread(file_img1);
imgIn2 = imread(file_img2);

% convert the image to png format 
file_img1_png = 'tmpASIFTinput1.png';
file_img2_png = 'tmpASIFTinput2.png';
imwrite(imgIn1, file_img1_png, 'png');
imwrite(imgIn2, file_img2_png, 'png');

% ASIFT command
command_ASIFT = './demo_ASIFT'; 
command_ASIFT = [command_ASIFT ' ' file_img1_png ' ' file_img2_png ' ' ...
  imgOutVert ' ' imgOutHori ' ' matchings ' ' keys1 ' ' keys2];
if (flag_resize == 0)
    command_ASIFT = [command_ASIFT ' 0'];
end
	
% get the number of processors 
% Mac
if (ismac == 1)
    [s, w] = unix('sysctl -n hw.ncpu');
    num_CPUs = str2num(w);
    % set the maximum OpenMP threads to the number of processors 
    set_threads = sprintf('export OMP_NUM_THREADS=%d;', num_CPUs);
     
    command = [set_threads ' ' command_ASIFT];

% Unix    
elseif (isunix == 1)
    [s, w] = unix('grep processor /proc/cpuinfo | wc -l');
     num_CPUs = str2num(w);
    % set the maximum OpenMP threads to the number of processors 
    set_threads = sprintf('export OMP_NUM_THREADS=%d;', num_CPUs);
      
    command = [set_threads ' ' command_ASIFT];
% Windows    
elseif (ispc == 1)
    [s, w] = dos('set NUMBER_OF_PROCESSORS');
    num_CPUs = sscanf(w, '%*21c%d', [1, Inf]);
    % set the maximum OpenMP threads to the number of processors 
    setenv('OMP_NUM_THREADS', num2str(num_CPUs));
        
    command = command_ASIFT;
else
    error('Unrecognized operating system. The operating system should be Windows, Linux/Unix, or Mac OS.');
end


status = system(command);
