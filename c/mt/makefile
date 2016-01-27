# Copyright 2010 Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
#
# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the copyright notice and this notice are preserved.  This file is
# offered as-is, without any warranty.

# source code, C language
CSRC	= mt.c example/rand.c
# source code, all languages (only C here)
SRC	= $(CSRC)
# object files (partial compilation)
OBJ	= $(CSRC:.c=.o)
# binary executable program
BIN	= example/rand

# standard C compiler optimization options
COPT	= -O2 -funroll-loops -fomit-frame-pointer
# complete C compiler options
CFLAGS	= -ansi -pedantic -Wall -Wextra -Werror -pipe $(COPT)

# default target: the binary executable program
default: $(BIN)

# partial C compilation xxx.c -> xxx.o
%.o	: %.c
	$(CC) $< -c $(CFLAGS) -I. -o $@

# final link of the partially compiled files
$(BIN)	: $(OBJ)
	$(CC) $(OBJ) -o $@

# cleanup
.PHONY	: clean distclean scrub
clean	:
	$(RM) $(OBJ)
	$(RM) *.flag
distclean	: clean
	$(RM) $(BIN)
	$(RM) -r srcdoc

################################################
# extra tasks

PROJECT	= mt
DATE	= $(shell date -u +%Y%m%d)
RELEASE_TAG   = 0.$(DATE)

.PHONY	: srcdoc lint beautify test release
# source documentation
srcdoc	: $(SRC)
	doxygen doc/doxygen.conf
# code cleanup
beautify	: $(CSRC)
	for FILE in $^; do \
		expand $$FILE | sed 's/[ \t]*$$//' > $$FILE.$$$$ \
		&& indent -kr -i4 -l78 -nut -nce -sob -sc \
			$$FILE.$$$$ -o $$FILE \
		&& rm $$FILE.$$$$; \
	done
# static code analysis
lint	: $(CSRC)
	for FILE in $^; do splint -ansi-lib -weak -I. $$FILE; done;
	for FILE in $^; do splint -posix-lib -weak -I. $$FILE; done;
	for FILE in $^; do clang --analyze -ansi -I. $$FILE; done; $(RM) *.plist
# code tests
test	: $(CSRC)
	sh -e test/run.sh && echo SUCCESS || ( echo ERROR; return 1)
# release tarball
release	: beautify lint test distclean
	git archive --format=tar --prefix=$(PROJECT)-$(RELEASE_TAG)/ HEAD \
	        | gzip > ../$(PROJECT)-$(RELEASE_TAG).tar.gz
