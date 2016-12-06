################################################################
## scPipe Makefile
## Version 0.1
## Author: Luyi Tian (tian.l@wehi.edu.au)
################################################################

DIR = $(shell pwd)
HAS_HTS = 1
HTS_DIR = $(DIR)/hts


#CXX = g++-5
#CC = gcc-5
CPPFLAGS = -std=c++11 -static-libstdc++ -O2
ifeq ($(HAS_HTS),1)
	LIBS = -lhts -lz
else
	LIBS = -I$(HTS_DIR)/include -L$(HTS_DIR)/lib -Wl,-rpath,$(HTS_DIR)/lib -lhts -lz
endif
OBJS = utils.o cellbarcode.o parsebam.o parsecount.o trimbarcode.o transcriptmapping.o
RM = rm -rf

%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(LIBS) -o $@ $<

ifeq ($(HAS_HTS),1)
all: buildPipe 
else
all: buildPipe buildhts
endif


buildhts:
	mkdir -p $(HTS_DIR) &&\
	cd htslib-1.3 && \
	make clean && \
	./configure --disable-libcurl --prefix=$(HTS_DIR) CC=$(CC) CXX=$(CXX) && \
	make && \
	make install && \
	cd ..

buildtest: $(OBJS)
	$(CXX) -o test/basic_test test/basic_test.cpp $(LIBS) $(OBJS) $(CPPFLAGS)

runtest:
	./test/basic_test

test: buildtest runtest

buildPipe: $(OBJS)
	mkdir -p bin
	$(CXX) -o bin/sc_demultiplex sc_demultiplex.cpp $(LIBS) $(OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_gene_counting sc_gene_counting.cpp $(LIBS) $(OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_trim_barcode sc_trim_barcode.cpp $(LIBS) $(OBJS) $(CPPFLAGS)
	$(CXX) -o bin/sc_exon_mapping sc_exon_mapping.cpp $(LIBS) $(OBJS) $(CPPFLAGS)


utils.o: utils.cpp
cellbarcode.o: cellbarcode.cpp utils.o
parsebam.o: parsebam.cpp cellbarcode.o utils.o
parsecount.o: parsecount.cpp cellbarcode.o utils.o
trimbarcode.o: trimbarcode.cpp utils.o
transcriptmapping.o: transcriptmapping.cpp utils.o

clean:
	$(RM) $(OBJS)
	$(RM) test/basic_test
	$(RM) ./bin