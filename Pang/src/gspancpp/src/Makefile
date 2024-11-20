#################################################################
###   Copyright (C) 2004 Taku Kudo, All rights reserved.
###     This is free software with ABSOLUTELY NO WARRANTY.
###  
###   This program is free software; you can redistribute it and/or modify
###     it under the terms of the GNU General Public License as published by
###     the Free Software Foundation; either version 2 of the License, or
###     (at your option) any later version.
###    
###   This program is distributed in the hope that it will be useful,
###     but WITHOUT ANY WARRANTY; without even the implied warranty of
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###     GNU General Public License for more details.
###    
###   You should have received a copy of the GNU General Public License
###     along with this program; if not, write to the Free Software
###     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
###     02111-1307, USA
###
###     modified by Marisa Thoma (extracted the original version from matlab; 
###                                                    extended to CORK pruning)
###    release: 12/21/2009
###    version: 1.0
#################################################################

CPP = g++
BIN = ../bin/

###
### The values below should not be modified
###

# DEBUG and RELEASE settings
CPPFLAGS = -DNDEBUG -O3 -g -ggdb -Wall -Wno-deprecated
IGRAPHLIB = -L/usr/local/lib -ligraph 
IGRAPHINCLUDE = -I/usr/local/include/igraph 

OBJ = ${BIN}graph.o ${BIN}misc.o ${BIN}dfs.o ${BIN}gspan.o ${BIN}gSpanCORK.o

### Targets

all:  ${BIN}gSpanCORK

${BIN}gSpanCORK: ${OBJ} gSpanCORK.cpp gspan.h
	${CPP} ${CPPFLAGS} ${IGRAPHLIB} ${IGRAPHINCLUDE} gSpanCORK.cpp ${OBJ} -o ${BIN}gSpanCORK

${BIN}gSpanCORK.o: gspan.cpp gspan.h
	${CPP} ${CPPFLAGS} ${IGRAPHINCLUDE} -c gSpanCORK.cpp -o ${BIN}gSpanCORK.o	
${BIN}gspan.o: gspan.cpp gspan.h
	${CPP} ${CPPFLAGS} ${IGRAPHINCLUDE} -c gspan.cpp -o ${BIN}gspan.o
${BIN}graph.o: graph.cpp gspan.h
	${CPP} ${CPPFLAGS} ${IGRAPHINCLUDE} -c graph.cpp -o ${BIN}graph.o
${BIN}misc.o: misc.cpp gspan.h
	${CPP} ${CPPFLAGS} ${IGRAPHINCLUDE} -c misc.cpp -o ${BIN}misc.o
${BIN}dfs.o: dfs.cpp gspan.h
	${CPP} ${CPPFLAGS} ${IGRAPHINCLUDE} -c dfs.cpp -o ${BIN}dfs.o

clean:
	rm -f ${BIN}*.o

