# Copyright 2014 Ives Rey-Otero <ivesreyotero@gmail.com>

# compilers configuration
CC = gcc
OFLAGS = -g -O3
LIBS = -lpng -L/usr/local/lib -lm

CFLAGS = -Wall -Wno-write-strings  -pedantic -std=c99 -D_POSIX_C_SOURCE=200809L

# Source files with executables.
SRC_ALGO = sift_cli

SRC_MATCH = match_cli

# TEMP normalized_patch
SRC_DEMO = demo_extract_patch

SRCa = lib_sift.c \
	   lib_sift_anatomy.c \
	   lib_scalespace.c \
	   lib_description.c \
       lib_discrete.c \
	   lib_keypoint.c \
	   lib_util.c

SRCb = lib_io_scalespace.c

SRCc = lib_matching.c

SRCDIR = src
OBJDIR = src
BINDIR = bin

OBJa = $(addprefix $(OBJDIR)/,$(SRCa:.c=.o))
OBJb = $(addprefix $(OBJDIR)/,$(SRCb:.c=.o))
OBJc = $(addprefix $(OBJDIR)/,$(SRCc:.c=.o))

OBJ = $(OBJa) $(OBJb) $(OBJc)

BIN = $(addprefix $(BINDIR)/,$(SRC_ALGO))
BINMATCH = $(addprefix $(BINDIR)/,$(SRC_MATCH))
BINDEMO = $(addprefix $(BINDIR)/,$(SRC_DEMO))

sift= $(BIN)
match= $(BINMATCH)
demo= $(BINDEMO)
default: $(OBJDIR) $(BINDIR) $(sift) $(match) $(demo)

#---------------------------------------------------------------
#  SIFT CLI
#

$(BIN) : $(BINDIR)/% : $(SRCDIR)/%.c $(OBJDIR)/lib_sift_anatomy.o  $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_scalespace.o $(OBJDIR)/lib_description.o   $(OBJDIR)/lib_discrete.o  $(OBJDIR)/lib_io_scalespace.o  $(OBJDIR)/lib_util.o   $(OBJDIR)/io_png.o $(OBJDIR)/lib_sift.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(OBJDIR):
	    -mkdir -p $(OBJDIR)

$(BINDIR):
	    -mkdir -p $(BINDIR)

#---------------------------------------------------------------
#  LIB_SIFT
#
$(OBJDIR)/lib_sift.o : $(SRCDIR)/lib_sift.c $(OBJDIR)/lib_sift_anatomy.o $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_scalespace.o : $(SRCDIR)/lib_scalespace.c $(OBJDIR)/lib_discrete.o  $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_discrete.o : $(SRCDIR)/lib_discrete.c $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_description.o : $(SRCDIR)/lib_description.c $(OBJDIR)/lib_discrete.o $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_keypoint.o : $(SRCDIR)/lib_keypoint.c $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_sift_anatomy.o : $(SRCDIR)/lib_sift_anatomy.c $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_discrete.o $(OBJDIR)/lib_scalespace.o $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_util.o : $(SRCDIR)/lib_util.c
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

#--------------------------------------------------------------
#   IN (image) and OUT (scalespace)
#
$(OBJDIR)/io_png.o : $(SRCDIR)/io_png.c
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/lib_io_scalespace.o : $(SRCDIR)/lib_io_scalespace.c $(OBJDIR)/io_png.o  $(OBJDIR)/lib_scalespace.o $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<


#-------------------------------------------------------------
#   Matching algorithm
#
$(OBJDIR)/lib_matching.o : $(SRCDIR)/lib_matching.c $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(BINMATCH) : $(SRCDIR)/match_cli.c $(OBJDIR)/lib_keypoint.o $(OBJDIR)/lib_matching.o   $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) -lm

#-------------------------------------------------------------
#  Tools used in the demo 
#
$(BINDEMO) : $(BINDIR)/% :	 $(SRCDIR)/demo_extract_patch.c  $(OBJDIR)/lib_discrete.o $(OBJDIR)/io_png.o $(OBJDIR)/lib_util.o 
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)


#-------------------------------------------------------------------------------------
#  clean
#
cleanobj:
	-rm -f $(OBJ) $(OBJDIR)/io_png.o

clean: cleanobj
	-rm -f $(BIN) $(BINMATCH) $(BINDEMO)
