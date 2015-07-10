# Copyright 2010 Nicolas Limare <nicolas.limare@cmla.ens-cachan.fr>
#
# Copying and distribution of this file, with or without
# modification, are permitted in any medium without royalty provided
# the copyright notice and this notice are preserved.  This file is
# offered as-is, without any warranty.
#
# This makefile imports and compiles libraries libraries
# for local use in static form

DEST_DIR = build
LIB_DIR	= $(DEST_DIR)/lib
INC_DIR	= $(DEST_DIR)/include

#
# COMMON
#

default	: all
all	: zlib libpng libtiff libjpeg

$(LIB_DIR) $(INC_DIR)	:
	mkdir -p $@
$(LIB_DIR)/% $(INC_DIR)/%	:
	cp $(filter %.h %.a, $^) $@


-include makefile.zlib
-include makefile.libpng
-include makefile.libjpeg
-include makefile.libtiff

#
# HOUSEKEEPING
#

.PHONY	: clean distclean scrub
clean	:
	$(RM) -r $(ZLIB_DIR)	
	$(RM) -r $(LIBPNG_DIR)
	$(RM) -r $(LIBJPEG_DIR)
	$(RM) -r $(LIBTIFF_DIR)
distclean	: clean
	$(RM) -r $(DEST_DIR)
scrub	: distclean
	$(RM) $(ZLIB_ARC)
	$(RM) $(LIBPNG_ARC)
	$(RM) $(LIBJPEG_ARC)
	$(RM) $(LIBTIFF_ARC)
