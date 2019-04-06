#pragma once
/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 *
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/

/* stack.c
 *
 * Implementation of generic stack structures.
 */

#include <stdlib.h>

#include "squid.hpp"


/************************************************************
 * intstack_s implementation.
 *
 * Functions: InitIntStack() - returns ptr to new stack
 *            PushIntStack() - (void)
 *            PopIntStack()  - returns 1 on success, 0 if stack empty
 *            FreeIntStack() - returns number of elements free'd, or 0 if
 *                             stack was empty.
 *           
 * Implementation of the pushdown stack for storing single
 * integers.
 *************************************************************/ 
struct intstack_s *
InitIntStack();


void
PushIntStack(
  struct intstack_s *stack,
  int data);


int
PopIntStack(
  struct intstack_s  *stack,
  int *ret_data);


void
ReverseIntStack(
  struct intstack_s *stack);


int
FreeIntStack(
  struct intstack_s *stack );