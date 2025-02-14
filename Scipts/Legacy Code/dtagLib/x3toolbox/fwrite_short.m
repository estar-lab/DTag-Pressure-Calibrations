function    n = fwrite_short(fid,p)

%    n = fwrite_short(fid,p)
%    Trick to avoid endian issues when writing 16-bit numbers
%    to a binary file in Matlab
%

p = [floor(p'/256);rem(p',256)] ;
fwrite(fid,p(:),'uchar') ;
return
