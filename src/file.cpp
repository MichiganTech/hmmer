/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"
#include "file.h"


char*
FileDirname(
  char *file
){
  char *dirname;
  char *lastslash;
  int   len;
  
  lastslash = strrchr(file, '/');
  len =  (lastslash == NULL) ? 0 : (int) (lastslash - file);
  dirname = (char *) MallocOrDie (sizeof(char) * (len+2));
  if (len > 0)                strncpy(dirname, file, len);
  else if (*file != '/') { *dirname = '.';      len = 1; }
  else                        { *dirname = '/'; len = 1; }
  dirname[len] = '\0';
  return dirname;
}

     
char* 
FileTail(
  char *file, 
  bool noextension
){
  char *tail;
  char *lastslash;
  char *lastdot;
				/* remove directory prefix */
  lastslash = strrchr(file, '/');
  tail = (char *) MallocOrDie (sizeof(char) * (strlen(file)+1));
  if (lastslash == NULL) strcpy(tail, file);
  else                   strcpy(tail, lastslash+1);
				/* remove trailing suffix */
  if (noextension) {
    if ((lastdot = strrchr(tail, '.')) != NULL)
      *lastdot = '\0';
  }

  return tail;
}


char*
FileSameDirectory(
  char *file1, 
  char *file2
){
  char *path;
  char *tail;
  char *result;
  int   seems_ok = 1;

  path  = FileDirname(file1);
  tail  = FileTail(file2, false);
  if (strcmp(file2, tail) != 0) seems_ok = 0; /* ut-oh, file2 *had* a path */
  result = FileConcat(path, tail);
  if (! seems_ok && strcmp(result, file2) != 0) {
    free(result); result = NULL; 
  }
  free(path);
  free(tail);
  return result;
}


char*
FileConcat(
  char *dir, 
  char *file
){
  char *full;

  full = (char *) MallocOrDie (sizeof(char) * (strlen(dir)+strlen(file)+2));
  if (*file == '/') strcpy(full, file); /* file = "/foo", ignore directory. */
  else                   sprintf(full, "%s%c%s", dir, '/', file);
  return full;
}


char*
FileAddSuffix(
  char *filename, 
  char *sfx
){
  char *new;
  new = MallocOrDie(strlen(filename) + strlen(sfx) + 2);
  sprintf(new, "%s.%s", filename, sfx);
  return new;
}


FILE*
EnvFileOpen(
  char *fname, 
  char *env, 
  char **ret_dir
){
  FILE *fp;
  char *path;
  char *s;                      /* ptr to indiv element in env list */
  char  full[1024];             /* constructed file name */

  if (env == NULL) return NULL;
  if ((path = strdup(getenv(env))) == NULL) return NULL;
  
  fp = NULL;
  s  = strtok(path, ":");
  while (s != NULL)
    {
      if (((int) strlen(fname) + (int) strlen(s) + 2) > 1024) 
	{ free(path); return NULL; }
      sprintf(full, "%s%c%s", s, '/', fname);
      if ((fp = fopen(full, "r")) != NULL) break;
      s = strtok(NULL, ":");
    }

  /* Return the path we used, if caller wants it
   */
  if (ret_dir != NULL) *ret_dir = strdup(s);
  free(path);
  
  return fp;
}


bool
FileExists(
  char *filename
){
  FILE *fp;
  if ((fp = fopen(filename, "r"))) { fclose(fp); return true; }
  return false;
}
