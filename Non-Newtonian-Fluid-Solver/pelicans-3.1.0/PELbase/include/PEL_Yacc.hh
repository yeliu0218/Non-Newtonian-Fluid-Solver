/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     PEL__IDENTIF = 258,
     PEL__STRING = 259,
     PEL__REAL = 260,
     PEL__INTEGER = 261,
     PEL__EOF = 262,
     PEL__ZERO = 263,
     PEL__MODULE = 264,
     PEL__END = 265,
     PEL__TRUE = 266,
     PEL__FALSE = 267,
     PEL_INCLUDE = 268,
     PEL__CONCAT = 269,
     PEL__OR = 270,
     PEL__AND = 271,
     PEL__IF = 272,
     PEL__LE = 273,
     PEL__GE = 274,
     PEL__NEQ = 275,
     PEL__LAST = 276,
     UMINUS = 277,
     UNOT = 278
   };
#endif
/* Tokens.  */
#define PEL__IDENTIF 258
#define PEL__STRING 259
#define PEL__REAL 260
#define PEL__INTEGER 261
#define PEL__EOF 262
#define PEL__ZERO 263
#define PEL__MODULE 264
#define PEL__END 265
#define PEL__TRUE 266
#define PEL__FALSE 267
#define PEL_INCLUDE 268
#define PEL__CONCAT 269
#define PEL__OR 270
#define PEL__AND 271
#define PEL__IF 272
#define PEL__LE 273
#define PEL__GE 274
#define PEL__NEQ 275
#define PEL__LAST 276
#define UMINUS 277
#define UNOT 278




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE PELlval;

