/*
Copyright (C) 2003-2006 Tommi Junttila

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License version 2
 as published by the Free Software Foundation.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef HEAP_HH
#define HEAP_HH

namespace igraph {

class Heap
{
#if defined(CONSISTENCY_CHECKS)
  unsigned int N;
#endif
  unsigned int n;
  unsigned int *array;
  void upheap(unsigned int k);
  void downheap(unsigned int k);
public:
  Heap() {array = 0; n = 0; }
  ~Heap();
  void init(unsigned int size);

  bool is_empty() const {return(n==0); }
  void clear() {n = 0;}
  void insert(unsigned int v);
  unsigned int remove();
};

}

#endif
