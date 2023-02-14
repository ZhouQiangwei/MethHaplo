CC = g++
FLAGS = -O2 
LDFLAGS = -O2
LIBS = -m64 -I./submodules/libbm/ -L./submodules/libbm/ -I./src/samtools-0.1.18/ -L./src/samtools-0.1.18/ -lbam -lz -lpthread -lgsl -lgslcblas -lBinaMeth
srcdir = ./src
builddir = ./build
subdir = ./submodules/libbm
PWD := $(shell pwd)

TARGET = $(srcdir)/methyhap
TOCOMPILE = $(builddir)/methyhap.o $(builddir)/paired.o $(builddir)/processPairedBlock.o -Wl,-rpath $(PWD)/$(subdir)

TARGET_asm = $(srcdir)/ASM
TOCOMPILE_asm = $(builddir)/ASM.o

TARGET_bsmerge = $(srcdir)/bsmerge
TOCOMPILE_bsmerge = $(builddir)/bsmerge.o

TARGET_bsmergehic = $(srcdir)/bsmergehic
TOCOMPILE_bsmergehic = $(builddir)/bsmergehic.o

TARGET_asman = $(srcdir)/ASManno
TOCOMPILE_asman = $(builddir)/ASManno.o

TARGET_asmsitesan = $(srcdir)/ASMannoSites
TOCOMPILE_asmsitesan = $(builddir)/ASMannoSites.o

TARGET_others = $(srcdir)/mergehap
TOCOMPILE_others = $(builddir)/mergehap.o

TARGET_homo = $(srcdir)/homometh
TOCOMPILE_homo = $(builddir)/homometh.o

TARGET_bam2md = $(srcdir)/bam2md
TOCOMPILE_bam2md = $(builddir)/bam2md.o -Wl,-rpath $(PWD)/$(subdir)

TARGET_filgenome = $(srcdir)/filgenome
TOCOMPILE_filgenome = $(builddir)/filgenome.o

TARGET_splitmr = $(srcdir)/splitmr
TOCOMPILE_splitmr = $(builddir)/splitmr.o -Wl,-rpath $(PWD)/$(subdir)

all: ${TOCOMPILE} ${TOCOMPILE_asm} ${TOCOMPILE_asman} ${TOCOMPILE_others} ${TOCOMPILE_homo} ${TOCOMPILE_asmsitesan} ${TOCOMPILE_bsmerge} ${TOCOMPILE_bsmergehic} ${TOCOMPILE_filgenome} ${TOCOMPILE_splitmr} ${TOCOMPILE_bam2md}
	make -C submodules/libbm
	${CC} $(LDFLAGS) -o $(TARGET_asm) ${TOCOMPILE_asm}
	${CC} $(LDFLAGS) -o $(TARGET_asman) ${TOCOMPILE_asman}
	${CC} $(LDFLAGS) -o $(TARGET_asmsitesan) ${TOCOMPILE_asmsitesan}
	${CC} $(LDFLAGS) -o $(TARGET_homo) ${TOCOMPILE_homo}
	${CC} $(LDFLAGS) -o $(TARGET_bsmerge) ${TOCOMPILE_bsmerge}
	${CC} $(LDFLAGS) -o $(TARGET_bsmergehic) ${TOCOMPILE_bsmergehic}
	${CC} $(LDFLAGS) -o $(TARGET_others) ${TOCOMPILE_others}
	${CC} $(LDFLAGS) -o $(TARGET_filgenome) ${TOCOMPILE_filgenome}
	${CC} $(LDFLAGS) -o $(TARGET_bam2md) ${TOCOMPILE_bam2md} ${LIBS}
	${CC} $(LDFLAGS) -o $(TARGET) ${TOCOMPILE} ${LIBS}
	${CC} $(LDFLAGS) -o $(TARGET_splitmr) ${TOCOMPILE_splitmr} ${LIBS}
	cd submodules/HapCUT2 && make

$(TOCOMPILE): $(srcdir)/methyhap.cpp $(srcdir)/methyhap.hpp $(srcdir)/paired.cpp $(srcdir)/paired.hpp $(srcdir)/processPairedBlock.cpp $(srcdir)/processPairedBlock.hpp $(srcdir)/common.hpp $(srcdir)/alignRead.h -Wl,-rpath $(PWD)/$(subdir) | $(builddir) $(subdir)
	$(CC) -c $(srcdir)/methyhap.cpp -o $(builddir)/methyhap.o
	$(CC) -c $(srcdir)/paired.cpp -o $(builddir)/paired.o
	$(CC) -c $(srcdir)/processPairedBlock.cpp -o $(builddir)/processPairedBlock.o

$(TOCOMPILE_asm): $(srcdir)/ASM.cpp
	$(CC) -c $(srcdir)/ASM.cpp -o $(builddir)/ASM.o

$(TOCOMPILE_filgenome): $(srcdir)/filgenome.cpp
	$(CC) -c $(srcdir)/filgenome.cpp -o $(builddir)/filgenome.o

$(TOCOMPILE_bam2md): $(srcdir)/bam2md.cpp
	$(CC) -c $(srcdir)/bam2md.cpp -o $(builddir)/bam2md.o -Wl,-rpath $(PWD)/$(subdir)

$(TOCOMPILE_splitmr): $(srcdir)/splitmr.cpp
	$(CC) -c $(srcdir)/splitmr.cpp -o $(builddir)/splitmr.o -Wl,-rpath $(PWD)/$(subdir)

$(TOCOMPILE_asman): $(srcdir)/ASManno.cpp
	$(CC) -c $(srcdir)/ASManno.cpp -o $(builddir)/ASManno.o

$(TOCOMPILE_asmsitesan): $(srcdir)/ASMannoSites.cpp
	$(CC) -c $(srcdir)/ASMannoSites.cpp -o $(builddir)/ASMannoSites.o

$(TOCOMPILE_homo): $(srcdir)/homometh.cpp
	$(CC) -c $(srcdir)/homometh.cpp -o $(builddir)/homometh.o

$(TOCOMPILE_others): $(srcdir)/mergehap.cpp
	$(CC) -c $(srcdir)/mergehap.cpp -o $(builddir)/mergehap.o

$(TOCOMPILE_bsmerge): $(srcdir)/bsmerge.cpp
	$(CC) -c $(srcdir)/bsmerge.cpp -o $(builddir)/bsmerge.o

$(TOCOMPILE_bsmergehic): $(srcdir)/bsmergehic.cpp
	$(CC) -c $(srcdir)/bsmergehic.cpp -o $(builddir)/bsmergehic.o

.c.o:
	$(CC) ${FLAGS} ${LIBS} -c $(srcdir)/$*.c

clean:
	rm -f $(builddir)/*.o $(TARGET) $(TARGET_asm) $(TARGET_asman) $(TARGET_homo) $(TARGET_others) $(TARGET_asmsitesan) $(TARGET_bsmerge) $(TARGET_bsmergehic) $(TARGET_bam2md) $(TARGET_filgenome) $(TARGET_splitmr)
	rm -rf submodules/HapCUT2/build/
	make clean -C submodules/libbm/
	if [ -d "bin" ]; then rm -r bin; fi

install:
	if [ -d "bin" ]; then echo bin exists; else mkdir bin; fi
	cp src/methyhap ./bin/
	cp src/ASM ./bin/
	cp src/ASManno ./bin/
	cp src/homometh ./bin/
	cp src/mergehap ./bin/
	cp src/ASMannoSites ./bin/
	cp src/bsmerge ./bin/
	cp src/bsmergehic ./bin/
	cp submodules/HapCUT2/build/HAPCUT2 ./bin/hapcut2
	cp submodules/HapCUT2/build/extractHAIRS ./bin/extracthairs
	cp scripts/methHaplo ./bin/
	cp scripts/bamStrand ./bin/
	cp src/bam2md ./bin/
	cp src/filgenome ./bin/
	cp src/splitmr ./bin/
