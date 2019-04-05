/*****************************************************************
 * SQUID - a library of functions for biological sequence analysis
 * Copyright (C) 1992-2002 Washington University School of Medicine
 * 
 *     This source code is freely distributed under the terms of the
 *     GNU General Public License. See the files COPYRIGHT and LICENSE
 *     for details.
 *****************************************************************/


/* file.c
 * 
 * File operation utilities, dealing with pathnames, directories,
 * and environment variables.
 * 
 * The goal is to have these be platform-independent but they
 * currently are UNIX-specific: i.e. this file is currently POSIX compliant
 * but it is NOT ANSI C compliant. (The sole offender is getenv().)
 *
 */

#include "squid.h"


/* Function: FileDirname()
 * 
 * Purpose:  Returns the path from a filename:
 *             "/foo/bar/baz"  -> "/foo/bar"
 *             "foo/bar"       -> "foo" 
 *             "foo"           -> "."
 *             "/"             -> "/" 
 *           i.e. the string will be non-NULL; it will
 *           contain the string up to but not including the
 *           last '/' character; returns "." if
 *           there are no '/' characters, and returns "/"
 *           if the last slash is the first character.
 *           Modeled on Tcl's "file dirname" command.
 *
 * Args:     file - name of file "/foo/bar/baz".
 *           
 * Return:   ptr to malloc'ed string "/foo/bar".          
 */
char *
FileDirname(
  char *file);


/* Function: FileTail()
 * 
 * Purpose:  Return everything after the DIRSLASH:
 *             "/foo/bar/baz.1"  -> "baz.1"
 *             "foo/bar"         -> "bar" 
 *             "foo"             -> "foo"
 *             "/"               -> "" 
 *           If noextension is TRUE, removes a trailing ".foo" extension
 *           too.
 *           
 * Args:     file        - name of file "/foo/bar/baz.1"         
 *           noextension - TRUE to also remove extensions
 *           
 * Return:   ptr to malloc'ed string "baz.1"          
 */      
char* 
FileTail(
  char *file, 
  bool noextension);


/* Function: FileSameDirectory()
 *
 * Purpose:  Given a path to one file, and the 
 *           name of another file in the same directory,
 *           concat the path from file1 onto file2, and
 *           return the result. Caller must free the ptr
 *           that's returned. 
 *           
 *           Written for SSI - SSI indices contain filenames
 *           without paths, and we will need to convert that
 *           to a full path.
 *
 * Args:     file1 - a path to a file, e.g. "/foo/bar/baz.1"
 *           file2 - a simple filename, e.g. "quux.2"
 *           
 * Returns:  path to file2: e.g. "/foo/bar/quux.2"
 *           Returns NULL if file2 already has a path, and the result
 *           would be a different place.
 */
char*
FileSameDirectory(
  char *file1, 
  char *file2);


/* Function: FileConcat()
 * 
 * Purpose:  Concatenate a directory path and a file name,
 *           returning a pointer to a malloc'ed string with the
 *           full filename. This isn't just a string concat,
 *           because we're careful about the dir slash.
 */
char *
FileConcat(
  char *dir, 
  char *file);


/* Function: FileAddSuffix()
 *
 * Purpose:  Add a suffix to a filename, return a malloc'ed
 *           string containing the new filename.sfx name.
 *           Example:
 *             FileAddSuffix("genbank", "ssi")
 *           returns "genbank.ssi".  
 */
char *
FileAddSuffix(
  char *filename, 
  char *sfx);


/* Function: EnvFileOpen()
 * 
 * Purpose:  Open a file, given a file name and an environment
 *           variable that contains a directory path. Files
 *           are opened read-only. Does not look at current directory
 *           unless "." is explicitly in the path specified by env.
 *           
 *           For instance: 
 *             fp = EnvFileOpen("BLOSUM45", "BLASTMAT", NULL);
 *           or:
 *             fp = EnvFileOpen("swiss", "BLASTDB", NULL);  
 *             
 *           Environment variables may contain a colon-delimited
 *           list of more than one path; e.g.
 *             setenv BLASTDB /nfs/databases/foo:/nfs/databases/bar
 *             
 *           Sometimes a group of files may be found in
 *           one directory; for instance, an index file with a
 *           database. The caller can EnvFileOpen() the main
 *           file, and ask to get the name of the
 *           directory back in ret_dir, so it can construct
 *           the other auxiliary file names and fopen() them. (If it called
 *           EnvFileOpen(), it might get confused by 
 *           file name clashes and open files in different
 *           directories.
 *             
 * Args:     fname   - name of file to open
 *           env     - name of environment variable containing path
 *           ret_dir - if non-NULL, RETURN: name of dir that was used.
 *  
 * Return:   FILE * to open file, or NULL on failure -- same as fopen()
 *           Caller must free ret_dir if it passed a non-NULL address.
 */
FILE *
EnvFileOpen(
  char *fname, 
  char *env, 
  char **ret_dir);


/* Function: FileExists()
 * 
 * Purpose:  Return TRUE if filename exists.
 *           Testing fopen() is the only possible platform-independent test
 *           I'm aware of.  
 */
bool
FileExists(
  char *filename);
