#ifndef CPP
#$(error CPP variable undefined)
#endif

CPP = $(shell pwd)
SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS =

CC = gcc
AR = ar
RANLIB = ranlib
CFLAGS = -g -w -O3 -Wsign-compare
#-fopenmp
#-Wall
## changed: = instaed of ?= above 4 lines
## -fopenmp for multi-threads
LIBS = -lm -lz
EXTRA_CFLAGS_PIC = -fpic
LDFLAGS_SUB = -lz
LDLIBS =
INCLUDES = 

# Create a simple test-program to check if gcc can compile with curl
tmpfile:=$(shell mktemp --suffix=.c)
$(file >$(tmpfile),#include <curl/curl.h>)
$(file >>$(tmpfile),int main() { return 0; })
HAVE_CURL:=$(shell $(CC) $(CFLAGS) $(EXTRA_CFLAGS_PIC) $(LIBS) -lcurl $(tmpfile) -o /dev/null >/dev/null 2>&1 && echo "YES")
$(shell rm $(tmpfile))

ifeq ($(HAVE_CURL),YES)
	# If yes, add the library
	LIBS += -lcurl
else
	# and if not, disable CURL specific code compilation
	CFLAGS += -DNOCURL
endif


prefix = /usr/local
includedir = $(prefix)/include
libdir = $(exec_prefix)/lib

.PHONY: all clean lib test doc

.SUFFIXES: .c .o .pico

all: lib $(PROGS)

lib: lib-static lib-shared

lib-static: libBinaMeth.a

lib-shared: libBinaMeth.so

doc:
	doxygen

OBJS = io.o bmValues.o bmRead.o bmStats.o bmWrite.o

.c.o:
	$(CC) -I. $(CFLAGS) $(INCLUDES) -c -o $@ $<

.c.pico:
	$(CC) -I. $(CFLAGS) $(INCLUDES) $(EXTRA_CFLAGS_PIC) -c -o $@ $<

libBinaMeth.a: $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)
	$(RANLIB) $@

libBinaMeth.so: $(OBJS:.o=.pico)
	$(CC) -shared -o $@ $(OBJS:.o=.pico) $(LDLIBS) $(LIBS)
#$(CC) -shared $(LDFLAGS) -o $@ $(OBJS:.o=.pico) $(LDLIBS) $(LIBS)

test/testLocal: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testLocal.c libBinaMeth.a $(LIBS)

test/testRemoteManyContigs: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testRemoteManyContigs.c libBinaMeth.a $(LIBS)

test/testRemote: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testRemote.c libBinaMeth.a $(LIBS)

test/testWrite: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testWrite.c libBinaMeth.a $(LIBS)

bmtools: libBinaMeth.so
	$(CC) -o $@ -I. -L. $(CFLAGS) bmtools.c -lBinaMeth $(LIBS) -Wl,-rpath $(CPP) -lpthread

bmDMR:
	g++ $(CFLAGS) -c -o regression.o regression.cpp
	g++ $(CFLAGS) -o bmDMR bmDMR.cpp regression.o -I. -L. -lBinaMeth -Wl,-rpath $(CPP) -lgsl -lgslcblas -lm -lz

bam2bm:
	g++ $(CFLAGS) bam2bm.cpp -o bam2bm -m64 -I. -L. -lz -lBinaMeth -Wl,-rpath $(CPP) $(LDFLAGS_SUB)

test/exampleWrite: libBinaMeth.so
	$(CC) -o $@ -I. -L. $(CFLAGS) test/exampleWrite.c -lBinaMeth $(LIBS) -Wl,-rpath .

test/testBigBed: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testBigBed.c libBinaMeth.a $(LIBS)

test/testIterator: libBinaMeth.a
	$(CC) -o $@ -I. $(CFLAGS) test/testIterator.c libBinaMeth.a $(LIBS)

install: bam2bm bmtools bmDMR

test: test/testLocal test/testRemote test/testWrite test/testLocal bmtools test/exampleWrite test/testRemoteManyContigs test/testBigBed test/testIterator
	./test/test.py

clean:
	rm -f *.o libBinaMeth.a libBinaMeth.so *.pico test/testLocal test/testRemote test/testWrite bmtools bmDMR bam2bm test/exampleWrite test/testRemoteManyContigs test/testBigBed test/testIterator example_output.bw

install-old: libBinaMeth.a libBinaMeth.so
	install -d $(prefix)/lib $(prefix)/include
	install libBinaMeth.a $(prefix)/lib
	install libBinaMeth.so $(prefix)/lib
	install *.h $(prefix)/include

