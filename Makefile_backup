#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#   Note: this Makefile was modified from the version distributed with
#   TAP-plugins, an excellent collection by Tom Szilagyi.
#   http://tap-plugins.sf.net


#####################################################################
# PLEASE CHANGE THIS to your preferred installation location!
#
# Change this if you want to install somewhere else. In particular

INSTALL_PLUGINS_DIR	=	/usr/lib64/ladspa/

# NO EDITING below this line is required
# if all you want to do is install and use the plugins.
#####################################################################



# GENERAL

CC		=	gcc
CFLAGS		=	-I. -O3 -Wall -fomit-frame-pointer -fstrength-reduce -funroll-loops -ffast-math -c -fPIC -DPIC
LDFLAGS		=	-nostartfiles -shared -Wl,-Bsymbolic -lc -lm -lrt

PLUGINS		=	autotalent.so

all: $(PLUGINS)

# RULES TO BUILD PLUGINS FROM C CODE

autotalent.so: autotalent.c ladspa.h
	$(CC) $(CFLAGS) autotalent.c mayer_fft.c
	$(CC) $(LDFLAGS) -o autotalent.so autotalent.o mayer_fft.o

# OTHER TARGETS

install: targets
	-mkdir -p		$(INSTALL_PLUGINS_DIR)
	cp *.so 		$(INSTALL_PLUGINS_DIR)

targets:	$(PLUGINS)

always:	

clean:
	-rm -f `find . -name "*.so"`
	-rm -f `find . -name "*.o"`
	-rm -f `find . -name "*~"`
