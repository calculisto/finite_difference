DEPENDENCIES_HEADERS=

PROJECT=finite_difference
LINK.o=${LINK.cc}
CXXFLAGS+=-Wall -Wextra -Wconversion -Wshadow $(foreach dir, ${DEPENDENCIES_HEADERS}, -I../${dir})

