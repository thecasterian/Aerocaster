CC = gcc
AR = ar r

CFLAGS = -Wall -Wextra -O2 -g -rdynamic
INCLUDE = -I../include
LIBS = -lcgns -lhdf5

SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c, ../build/%.o, $(SRCS))

TARGET = ../build/libaerocaster.a

.PHONY: clean mkbuild

all: mkbuild $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $@ $(OBJS)

$(OBJS): ../build/%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE) $(LIBS)

clean:
	rm -rf ../build

mkbuild: ../build

../build:
	mkdir -p ../build