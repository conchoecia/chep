CPPF=-O3 -std=c++11

all: bin bin/chep_pileup_to_array bin/chep_windowed_het bin/chep_plot bin/chep_gff2all bin/chep_genome_analysis bin/merge_pileups

bin:
	mkdir bin

bin/chep_pileup_to_array: src/pileup_to_array.cpp src/pileup_to_array.hpp
	g++ ${CPPF} -o bin/chep_pileup_to_array src/pileup_to_array.cpp

bin/chep_windowed_het: src/windowed_het.cpp src/windowed_het.hpp
	g++ ${CPPF} -o bin/chep_windowed_het src/windowed_het.cpp

bin/merge_pileups: src/merge_pileups.cpp
	g++ ${CPPF} -o bin/merge_pileups src/merge_pileups.cpp -lz

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

bin/chep_genome_analysis: scripts/snakemake_scripts/Snakefile_chep_genome_analysis scripts/bash_scripts/chep_genome_analysis.sh
	chmod ugo+x scripts/bash_scripts/chep_genome_analysis.sh; \
	cd bin; \
	ln -s ../scripts/bash_scripts/chep_genome_analysis.sh chep_genome_analysis;\
	cd ..

#clean:
#   rm hello.o hello.exe
