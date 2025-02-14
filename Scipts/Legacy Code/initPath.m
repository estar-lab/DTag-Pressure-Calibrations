function initPath(recdir)

if ispc
    settagpath('cal',[recdir, '\cal']);
    %settagpath('raw','c:/tag/data/raw');
    settagpath('prh',[recdir, '\prh']);
    settagpath('audit',[recdir, '\audit']) ;
else
    settagpath('cal',[recdir, '/cal']);
    %settagpath('raw','c:/tag/data/raw');
    settagpath('prh',[recdir, '/prh']);
    settagpath('audit',[recdir, '/audit']) ;
end

end