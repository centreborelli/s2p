# compilers configuration
CC = gcc
OFLAGS = -O3
LIBS = -lpng -lm
#-lefence -ltiff -ljpeg

#CFLAGS = -g3 -rdynamic -Wall -Werror -Wno-write-strings  -Wno-deprecated -pedantic -std=c99 -D_POSIX_C_SOURCE=200809L
CFLAGS = -Wall -Werror -Wno-write-strings  -pedantic -std=c99 -D_POSIX_C_SOURCE=200809L

SRC_ALGO = sift matching extract_thumbnail_bis

#source files, for object files.
SRC2 = sift_scalespace.c sift_description.c \
       discrete_representation.c sift_matching.c std_image_lib.c sift_main.c sift_keypoint.c

SRCDIR = src
OBJDIR = src
BINDIR = bin

BIN = $(addprefix $(BINDIR)/,$(SRC_ALGO)) 
#BINDEMO = $(addprefix $(BINDIR)/,$(SRC_DEMO))
OBJ = $(addprefix $(OBJDIR)/,$(SRC2:.c=.o)) $(OBJDIR)/io_png.o 

default: $(BIN)
#demo : $(BINDEMO)

# binary files
$(BIN) : $(BINDIR)/% : $(SRCDIR)/%.c $(OBJ)
	-mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

$(BINDEMO) : $(BINDIR)/% : $(SRCDIR)/%.c $(OBJ)
	-mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

# object files with specific prerequisites
$(OBJDIR)/sift_scalespace.o : $(SRCDIR)/sift_scalespace.c $(OBJDIR)/discrete_representation.o  #$(OBJDIR)/sift_param.o #$(OBJDIR)/interpolation_DCT.o
	-mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/sift_matching.o : $(SRCDIR)/sift_matching.c $(OBJDIR)/sift_keypoint.o
	-mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

$(OBJDIR)/std_image_lib.o : $(SRCDIR)/std_image_lib.c $(OBJDIR)/io_png.o 
	-mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<


$(OBJDIR)/io_png.o : $(SRCDIR)/io_png.c
	-mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

	
$(OBJDIR)/%.o : $(SRCDIR)/%.c $(OBJDIR)/io_png.o
	-mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(OFLAGS) -c -o $@ $<

	
	
cleanobj:
	-rm -f $(OBJ)
	
clean: cleanobj
	-rm -f $(BIN)
