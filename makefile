all: bin/cOMet

bin/cOMet:
	mkdir bin
	g++ src/*.cpp -o bin/cOMet

clean: 
	rm -rf bin
	
parallel:
	mkdir bin
	g++ src_parallel/*.cpp -o bin/cOMet
