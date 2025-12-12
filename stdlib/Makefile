include Makefile.config

C = gcc
ASM = nasm
CPP = g++
LINK = ld

BIN = bin/
STATIC = bin/static/
INT = int/
SRC = src/
TEST = test/
INC = inc/

C_FLAGS = -Wall -std=c99 $(C_DEFINES) -Iinc
CPP_FLAGS = -Wall -std=c++23 $(CPP_DEFINES) -Iinc
ASM_FLAGS = -Wall -f elf64 $(ASM_DEFINES) -Iinc


default: prep math cpu string vector
	#$(LINK) --shared -o $(BIN)stdlib.so \
	#	$(STATIC)math.o \
	#	$(STATIC)cpu.o \
	#	$(STATIC)string.o
	$(LINK) --relocatable -o $(BIN)stdlib.o \
		$(STATIC)math.o \
		$(STATIC)cpu.o \
		$(STATIC)string.o \
		$(STATIC)vector.o

prep:
	mkdir -p int
	mkdir -p bin
	mkdir -p bin/static

clean:
	rm -rf bin
	rm -rf int

math:
	$(ASM) $(ASM_FLAGS) -o $(INT)arithmetic.o $(SRC)arithmetic.asm
	$(ASM) $(ASM_FLAGS) -o $(INT)trig.o $(SRC)trig.asm
	$(ASM) $(ASM_FLAGS) -o $(INT)exp.o $(SRC)exp.asm
	$(LINK) --relocatable -o $(STATIC)math.o \
		$(INT)arithmetic.o $(INT)trig.o $(INT)exp.o

cpu:
	$(ASM) $(ASM_FLAGS) -o $(INT)cpuid.o $(SRC)cpuid.asm
	$(LINK) --relocatable -o $(STATIC)cpu.o \
    		$(INT)cpuid.o

string:
	$(ASM) $(ASM_FLAGS) -o $(INT)string.o $(SRC)string.asm
	$(LINK) --relocatable -o $(STATIC)string.o \
    		$(INT)string.o

vector:
	$(ASM) $(ASM_FLAGS) -o $(INT)vector.o $(SRC)vector.asm
	$(LINK) --relocatable -o $(STATIC)vector.o \
    		$(INT)vector.o

all: default
	$(C) $(TEST)main.c $(BIN)stdlib.o $(C_FLAGS) -o $(BIN)main