function evalScript(fileNameScript)

    disp(['Executing ' fileNameScript])
    warning off
    eval(['run ' fileNameScript])
end