/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <igraph.h>
#include <stdlib.h>

void print_vector(igraph_vector_t *v, FILE *f) {
  long int i;
  for (i=0; i<igraph_vector_size(v); i++) {
    fprintf(f, " %li", (long int) VECTOR(*v)[i]);
  }
  fprintf(f, "\n");
}

int main() {
  
  igraph_vector_t v1, v2, v3;
  igraph_real_t min, max;
  long int imin, imax;
  int i;
  
  igraph_vector_init_seq(&v1, 1, 10);
  igraph_vector_init_seq(&v2, 0, 9);
  
  igraph_vector_swap(&v1, &v2);
  print_vector(&v1, stdout);
  print_vector(&v2, stdout);
  
  igraph_vector_swap_elements(&v1, 0, 9);
  igraph_vector_swap_elements(&v1, 3, 6);
  print_vector(&v1, stdout);
  
  igraph_vector_reverse(&v2);
  print_vector(&v2, stdout);
  igraph_vector_reverse(&v2);
  print_vector(&v2, stdout);
  
  igraph_vector_destroy(&v1);
  igraph_vector_destroy(&v2);
  
  igraph_vector_init(&v1, 10);
  igraph_vector_init(&v2, 10);
  igraph_vector_fill(&v1, 4);
  igraph_vector_fill(&v2, 2);
  
  igraph_vector_add(&v1, &v2);
  print_vector(&v1, stdout);
  igraph_vector_sub(&v1, &v2);
  print_vector(&v1, stdout);
  igraph_vector_div(&v1, &v2);
  print_vector(&v1, stdout);
  igraph_vector_mul(&v1, &v2);
  print_vector(&v1, stdout);
  
  igraph_vector_minmax(&v1, &min, &max);
  igraph_vector_which_minmax(&v1, &imin, &imax);
  printf("%g %g %i %i\n", min, max, imin, imax);
  
  igraph_vector_destroy(&v1);
  igraph_vector_destroy(&v2);
  
  igraph_vector_init_seq(&v1, 1, 10);
  igraph_vector_init(&v2, 10);
  for (i=0; i<10; i++) { VECTOR(v2)[i] = 10-i; }

  igraph_vector_minmax(&v1, &min, &max);
  igraph_vector_which_minmax(&v1, &imin, &imax);
  printf("%g %g %i %i\n", min, max, imin, imax);
  igraph_vector_minmax(&v2, &min, &max);
  igraph_vector_which_minmax(&v2, &imin, &imax);
  printf("%g %g %i %i\n", min, max, imin, imax);
  
  if (igraph_vector_isnull(&v1)) {
    return 1;
  }
  igraph_vector_null(&v1);
  if (!igraph_vector_isnull(&v1)) {
    return 2;
  }

  igraph_vector_destroy(&v1);
  igraph_vector_destroy(&v2);

  igraph_vector_init_int(&v1, 10, 3, 5, 6, 6, 6, 7, 8, 8, 9, 10);
  igraph_vector_init_int(&v2, 10, 1, 3, 3, 6, 6, 9, 12, 15, 17, 20);
  igraph_vector_init(&v3, 0);

  igraph_vector_intersect_sorted(&v1, &v2, &v3, 0);
  print_vector(&v3, stdout);
  igraph_vector_intersect_sorted(&v1, &v2, &v3, 1);
  print_vector(&v3, stdout);

  igraph_vector_destroy(&v1);
  igraph_vector_destroy(&v2);
  igraph_vector_destroy(&v3);

  if (IGRAPH_FINALLY_STACK_SIZE() != 0) return 3; 

  return 0;
}
