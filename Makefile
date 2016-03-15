all: bin/fractla



dist/main.o: dist/ src/main.cpp
	g++ -std=c++11 -I./src/ -c src/main.cpp -o dist/main.o

bin/fractla: bin/ dist/main.o
	g++ dist/main.o -lm -o bin/fractla

bin/: 
	mkdir bin

dist/:
	mkdir dist

clean: 
	rm -rf dist/*.o bin/pngm
	rm -R dist bin


