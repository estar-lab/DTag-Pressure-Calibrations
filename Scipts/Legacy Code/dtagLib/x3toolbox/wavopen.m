function    wh = wavopen(fname,fs,nch)

%     wh = wavopen(fname,fs,nch)
%     Generate a new wav format file called fname.
%     fs is the sampling rate for the file.
%
%     Returns:
%     wh is a structure containing the file handle
%     and a cumulative sample count which will be used
%     to close the file.
%
%     This is a quick hack relying on the Matlab wavwrite
%     to create the base file. Only 16-bit wav files are
%     supported.
%
%     GNU General Public License, see README.txt
%     mark johnson, SOI / Univ. St. Andrews
%     mj26@st-andrews.ac.uk
%     November 2012

% create a wav file by writing some dummy data
wavwrite(zeros(10,nch),fs,16,fname) ;
wh.fid = fopen(fname,'r+','l') ;

% move file cursor back to the start of the data chunk
fseek(wh.fid,-10*2*nch,'eof') ;

wh.nch = nch ;
wh.ns = 0 ;
