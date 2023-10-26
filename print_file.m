%funcion para imprimir ficheros con formato
function print_file( file , mat )

[row col] = size (mat);

myfile = fopen( file , 'a' );

for i = 1:row
	for j = 1:col
	fprintf( myfile , '%e\t' , mat(i,j) );
	end
fprintf( myfile , '\n' );
end

fclose( myfile );
