ip = 'localhost';
	port = 20320;
	% Create and open a UDP object that connects to BCI2000.
	u = udp( ip, 20319, 'LocalPort', port, 'Terminator', 'CR/LF', 'Timeout', 10 );
	fopen( u );
	% Read data until timeout occurs.
	s = fgetl( u );
	while( s~=-1 )
	  fprintf( '%s', s );
	  s = fgetl( u );
	end
	% Close and delete the UDP object.
	fclose( u );
	delete( u );