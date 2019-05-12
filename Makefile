CC = g++
FLAGS = -O2 
LDFLAGS = -O2
LIBS = -m64 -I./src/samtools-0.1.18/ -L./src/samtools-0.1.18/ -lbam -lz -lpthread -lgsl -lgslcblas
srcdir = ./src
builddir = ./build

TARGET = $(srcdir)/methyhaplo
TOCOMPILE = $(builddir)/methyhaplo.o $(builddir)/paired.o $(builddir)/processPairedBlock.o

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

all: ${TOCOMPILE} ${TOCOMPILE_asm} ${TOCOMPILE_asman} ${TOCOMPILE_others} ${TOCOMPILE_homo} ${TOCOMPILE_asmsitesan} ${TOCOMPILE_bsmerge} ${TOCOMPILE_bsmergehic}
	${CC} $(LDFLAGS) -o $(TARGET_asm) ${TOCOMPILE_asm}
	${CC} $(LDFLAGS) -o $(TARGET_asman) ${TOCOMPILE_asman}
	${CC} $(LDFLAGS) -o $(TARGET_asmsitesan) ${TOCOMPILE_asmsitesan}
	${CC} $(LDFLAGS) -o $(TARGET_homo) ${TOCOMPILE_homo}
	${CC} $(LDFLAGS) -o $(TARGET_bsmerge) ${TOCOMPILE_bsmerge}
	${CC} $(LDFLAGS) -o $(TARGET_bsmergehic) ${TOCOMPILE_bsmergehic}
	${CC} $(LDFLAGS) -o $(TARGET_others) ${TOCOMPILE_others}
	${CC} $(LDFLAGS) -o $(TARGET) ${TOCOMPILE} ${LIBS}
	cd submodules/HapCUT2 && make

$(TOCOMPILE): $(srcdir)/methyhaplo.cpp $(srcdir)/methyhaplo.hpp $(srcdir)/paired.cpp $(srcdir)/paired.hpp $(srcdir)/processPairedBlock.cpp $(srcdir)/processPairedBlock.hpp $(srcdir)/common.hpp $(srcdir)/alignRead.h | $(builddir)
	$(CC) -c $(srcdir)/methyhaplo.cpp -o $(builddir)/methyhaplo.o
	$(CC) -c $(srcdir)/paired.cpp -o $(builddir)/paired.o
	$(CC) -c $(srcdir)/processPairedBlock.cpp -o $(builddir)/processPairedBlock.o

$(TOCOMPILE_asm): $(srcdir)/ASM.cpp
	$(CC) -c $(srcdir)/ASM.cpp -o $(builddir)/ASM.o

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
	rm -f $(builddir)/*.o $(TARGET) $(TARGET_asm) $(TARGET_asman) $(TARGET_homo) $(TARGET_others) $(TARGET_asmsitesan) $(TARGET_bsmerge) $(TARGET_bsmergehic)
	rm -rf submodules/HapCUT2/build/

install:
	cp src/methyhaplo ./
	cp src/ASM ./
	cp src/ASManno ./
	cp src/homometh ./
	cp src/mergehap ./
	cp src/ASMannoSites ./
	cp src/bsmerge ./
	cp src/bsmergehic ./
	cp submodules/HapCUT2/build/HAPCUT2 ./hapcut2
	cp submodules/HapCUT2/build/extractHAIRS ./extracthairs
