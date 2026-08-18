include ../config.mk

CXXFLAGS+=-g

.PHONY: all clean check

all: check

check: tests
	./tests

tests: tests.o test_finite_difference.o

clean: 
	rm -f tests *.o 
