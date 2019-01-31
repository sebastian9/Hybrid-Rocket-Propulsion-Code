function [result] = NASACEAinp1 (OF, Pc, Pip)
    fid = fopen ('test.inp', 'w');
    fprintf (fid, ' prob ro equilibrium \n\n');
    fprintf (fid, '  ! iac problem \n');
    fprintf (fid, ' o/f %f \n', OF);
    fprintf (fid, ' p,atm = %f \n', Pc);
    fprintf (fid, ' pip %f \n', Pip);
    fprintf (fid, ' reac \n');
    fprintf (fid, '   fuel  C32H66(a) wt%%=100 t,k=298.15 \n');
    fprintf (fid, '   oxid  N2O wt%%=100.  t,k=298.15 \n');
    fprintf (fid, ' output trace=1e-51 \n');
    fprintf (fid, ' end \n');
    fclose(fid);
    
    f = fopen('temp', 'w');
    fprintf (f, 'test\n');
    fclose (f);
    
    dos('FCEA2.exe < temp');
    
    text = fileread('test.out');
    text = extractAfter(text,strfind(text,'Ae/At')+35);
    text = extractBefore(text,strfind(text,'CSTAR')-3);
    AEAT = str2double(text);
    text = fileread('test.out');
    text = extractAfter(text,strfind(text,'CSTAR')+35);
    text = extractBefore(text,strfind(text,'CF')-3);
    CSTAR = str2double(text);
    text = fileread('test.out');
    text = extractAfter(text,strfind(text,'CF')+35);
    text = extractBefore(text,strfind(text,'Ivac')-3);
    CF = str2double(text);
    
    [result] = [AEAT, CSTAR, CF];
end