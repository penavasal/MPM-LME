%funcion para imprimir ficheros con formato
function print_file_tens( file1 , file2, obj )

    col=size(obj,2);
    
    myfile1 = fopen( file1 , 'a' );
    myfile2 = fopen( file2 , 'a' );
    
    for i = 1:3
        for j = 1:col
            for k=1:3
	            %fprintf( myfile1 , '%e\t' , obj(j).ZPiola1new(i,k) );
                fprintf( myfile1 , '%e\t' , obj(j).ZCauchynew(i,k) );
                %fprintf( myfile2 , '%e\t' , obj(j).ZDeformgradnew(i,k) );
                fprintf( myfile2 , '%e\t' , obj(j).ZEpsilonNew(i,k) );
            end
        end
    
        fprintf( myfile1 , '\n' );
        fprintf( myfile2 , '\n' );
    end
    
    fclose( myfile1 );
    fclose( myfile2 );

end