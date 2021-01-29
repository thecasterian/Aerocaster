CC = gcc
AR = ar r

CFLAGS = -Wall -Wextra -O2
INCLUDE = -I../include -I/usr/include/glib-2.0

SRCS = $(wildcard *.c)
OBJS = $(patsubst %.c, ../build/%.o, $(SRCS))

TARGET = ../build/libaerocaster.a

.PHONY: clean mkbuild

all: mkbuild $(TARGET)

$(TARGET): $(OBJS)
	$(AR) $@ $(OBJS)

$(OBJS): ../build/%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE)

clean:
	rm -rf ../build

mkbuild: ../build

../build:
	mkdir -p ../build