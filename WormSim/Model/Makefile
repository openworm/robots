# WormSim makefile

all: worm wormpy

CC=g++
CXX=$(CC) -O2 -fPIC -I ../sundials-2.3.0/src/ida -I ../sundials-2.3.0/include -I ../sundials-2.3.0/src -Wno-deprecated
LIBS=-L /usr/local/lib -l sundials_ida -l sundials_cvode -l sundials_nvecserial -lm
SRCS=$(shell ls *.cpp)
OBJS=$(subst .cpp,.o,$(SRCS))

worm: $(OBJS)
	$(CXX) $(OBJS) $(LIBS) -o worm

wormpy: $(OBJS) wormpy.i
	swig -c++ -python wormpy.i
	$(CC) -O2 -fPIC -c wormpy_wrap.cxx -I/home/tom/anaconda2/include/python2.7
	$(CXX) -shared $(OBJS) wormpy_wrap.o $(LIBS) -o _wormpy.so

check:
	@echo SRCS = $(SRCS)
	@echo OBJS = $(OBJS)

clean:
	rm -rf $(TARGET) $(OBJS)

depend:
	$(CXX) -MM $(SRCS) > .depends.mk

-include .depends.mk
