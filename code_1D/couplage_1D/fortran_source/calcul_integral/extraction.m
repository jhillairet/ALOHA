    % extraction des données du fichier S_plasma.dat


    clear all
    format long
    
    long = 198*4;

    fichier_data = fopen('data_W.dat','r');

    fread(fichier_data,1,'int32');
    
    
    for ind = 1:long
     
            W(ind) = fread(fichier_data,1,'double');

    end

 
    fclose(fichier_data);
    
    ai=W(1:long/4);
    bi=W(1+long/4:2*long/4);
    ei=W(1+2*long/4:3*long/4);
    ri=W(1+3*long/4:4*long/4);
    
 
