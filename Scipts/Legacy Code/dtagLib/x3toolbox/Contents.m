% X3 Toolbox for Matlab, Version 1.0, 14 March 2013
%
% Copyright and license information
%  README.txt
%  GPL_COPYING.txt
%
% Archive reading and writing functions:
%  wav_to_x3   Generates an X3 archive from a wav file.
%  x3_to_wav   Unpacks an X3 archive to generate a wav and metadata xml file.
%
% File i/o functions:
%  x3new       Creates a new X3 archive and initializes the metadata
%  x3open      Opens an existing X3 archive and extracts the configuration metadata
%  fwrite_short   Endian-independent 16-bit binary file write
%  getnextframe   Extract the next frame of data from an X3 archive
%  wavopen     Create a new wav file
%  wavappend   Add audio data to a wav file
%  wavclose    Complete and close a wav file
%  readx3xml   Read an X3 metadata file or string and convert into a structure
%
% Functions for handling frames of X3 data
%  x3makeframe    Compress a frame of audio data and add it to an X3 archive 
%  x3makemetaframe   Add a frame of metadata to an X3 archive
%  x3decodeframe  Decode a frame of audio or metadata from an X3 archive
%  parsemetadata  Extract and act on configuration information in a metadata frame 
%
% The actual compression and uncompression functions:
%  x3compress  Compresses a matrix of audio data into packed integers.
%  x3uncompress  Uncompresses packed integers into a matrix of audio data.
%
% Utilities called by the above functions:
%  bpackv      Vectorized bit packer
%  bunpackv    Vectorized bit unpacker
%  bunpackrice Bit unpacker for variable length codes
%  makericecode   Generate a structure of information about a Rice code
%  crc16       Compute the 16-bit CCITT cyclic redundancy code for a vector
%  makeccitt16tab Compute a lookup table for crc16
%  gettime     Example function to make a time field
%
% XML4MATv2    XML parser (see separate copyright notice in this directory).
%
% To use, add the x3 directory and sub-directories to the Matlab path.
