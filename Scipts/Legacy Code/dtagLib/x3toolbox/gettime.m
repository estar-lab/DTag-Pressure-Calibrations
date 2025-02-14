function    T = gettime(dvec)

%     T = gettime(dvec)
%     Example function for generating an X3 time field. The X3 time field
%     comprises 4 x 16-bit words but is not currently defined. It could
%     contain a cumulative sample count or a time-of-day with microsecond
%     resolution. This function produces the latter type of time field by
%     converting a 6-element datevector [year month day hour minute second]
%     into (i) a 32 bit UNIX datenumber, i.e., the number of seconds 
%     since midnight Jan 1 1970, and (ii) the fractional seconds in 
%     microseconds. If no argument is given, the current time is returned.
%

if nargin==0,
   dvec = clock ;
end

dd = datenum([1970 1 1 0 0 0]) ;
ut = round((datenum(round(dvec))-dd)*3600*24) ;
ms = round(1e6*rem(dvec(6),1)) ;

% convert ut and ms to 4 x 16-bit words
T = zeros(4,1) ;
T(1) = floor(ut*2^(-16)) ;
T(2) = rem(ut,2^16) ;
T(3) = floor(ms*2^(-16)) ;
T(4) = rem(ms,2^16) ;
return
