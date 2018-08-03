all: bin/cOMet

bin/cOMet:
	mkdir bin
	g++ src/*.cpp -o bin/cOMet

clean: 
	rm -rf bin
	
parallel:
	mkdir bin
	mkdir bin/results
	g++ src_parallel/*.cpp -o bin/cOMet
