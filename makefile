all: bin/cOMet

clean: 
	rm -rf bin

bin/cOMet:
	mkdir bin
	g++ src/*.cpp -o bin/cOMet

