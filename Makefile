CPPF=-O3 -std=c++11

all: bin bin/chep_pileup_to_array bin/chep_windowed_het bin/chep_plot bin/chep_gff2all

bin:
	mkdir bin

bin/chep_pileup_to_array: src/pileup_to_array.cpp src/pileup_to_array.hpp
	g++ ${CPPF} -o bin/chep_pileup_to_array src/pileup_to_array.cpp

bin/chep_windowed_het: src/windowed_het.cpp src/windowed_het.hpp
	g++ -o bin/chep_windowed_het src/windowed_het.cpp

bin/chep_plot: scripts/heterozygosity_matrix.py
	chmod ugo+x scripts/heterozygosity_matrix.py; \
	cd bin; \
	ln -s ../scripts/heterozygosity_matrix.py chep_plot;\
	cd ..

bin/chep_gff2all: scripts/gff_to_intron_bed.py
	chmod ugo+x scripts/gff_to_intron_bed.py; \
	cd bin; \
	ln -s ../scripts/gff_to_intron_bed.py chep_gff2all;\
	cd ..

#clean:
#   rm hello.o hello.exe
