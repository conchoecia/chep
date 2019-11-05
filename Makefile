CPPF=-O3

all: bin bin/chep_pileup_to_array bin/chep_windowed_het bin/chep_plot

bin:
	mkdir bin

bin/chep_pileup_to_array: src/pileup_to_array.cpp src/pileup_to_array.hpp
	g++ ${CPPF} -o bin/pileup_to_array src/pileup_to_array.cpp

bin/chep_windowed_het: src/windowed_het.cpp src/windowed_het.hpp
	g++ -o bin/chep_windowed_het src/windowed_het.cpp

bin/chep_plot: scripts/heterozygosity_matrix.py
	chmod ugoo+x scripts/heterozygosity_matrix.py; \
	cd bin; \
	ln -s ../scripts/heterozygosity_matrix.py chep_plot;\
	cd ..
     
#clean:
#   rm hello.o hello.exe
