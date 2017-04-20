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
     YYDOC_COMMENT = 258,
     YYDOC_INCLUDE = 259,
     YYDOC_CLASS = 260,
     YYDOC_IDENTIF = 261,
     YYDOC_CATEGORY = 262,
     YYDOC_PUBLIC = 263,
     YYDOC_PROTECTED = 264,
     YYDOC_PRIVATE = 265,
     YYDOC_CONST = 266,
     YYDOC_VOID = 267,
     YYDOC_VIRTUAL = 268,
     YYDOC_STATIC = 269,
     YYDOC_OPERATOR = 270,
     YYDOC_TYPE_BASE = 271,
     YYDOC_PEL_ASSERT = 272,
     YYDOC_PEL_CHECK_PRE = 273,
     YYDOC_PEL_CHECK_POST = 274,
     YYDOC_PEL_LABEL = 275,
     YYDOC_USING = 276,
     YYDOC_FRIEND = 277,
     YYDOC_MUTABLE = 278,
     YYDOC_CHIFFRE = 279,
     YYDOC_ENUM = 280,
     YYDOC_STRUCT = 281,
     YYDOC_CHAINE = 282,
     YYDOC_FLECHE = 283,
     YYDOC_FORALLEXISTS = 284,
     YYDOC_TYPEDEF = 285,
     YYDOC_CAST = 286,
     YYDOC_POINTPOINT = 287,
     YYDOC_NEW = 288,
     YYDOC_EXTERN = 289,
     YYDOC_EGALEGAL = 290,
     YYDOC_DIFFEGAL = 291,
     YYDOC_INFEGAL = 292,
     YYDOC_SUPEGAL = 293,
     YYDOC_OROR = 294,
     YYDOC_ANDAND = 295,
     UMINUS = 296,
     NOT = 297,
     COMP = 298,
     TYPE = 299
   };
#endif
/* Tokens.  */
#define YYDOC_COMMENT 258
#define YYDOC_INCLUDE 259
#define YYDOC_CLASS 260
#define YYDOC_IDENTIF 261
#define YYDOC_CATEGORY 262
#define YYDOC_PUBLIC 263
#define YYDOC_PROTECTED 264
#define YYDOC_PRIVATE 265
#define YYDOC_CONST 266
#define YYDOC_VOID 267
#define YYDOC_VIRTUAL 268
#define YYDOC_STATIC 269
#define YYDOC_OPERATOR 270
#define YYDOC_TYPE_BASE 271
#define YYDOC_PEL_ASSERT 272
#define YYDOC_PEL_CHECK_PRE 273
#define YYDOC_PEL_CHECK_POST 274
#define YYDOC_PEL_LABEL 275
#define YYDOC_USING 276
#define YYDOC_FRIEND 277
#define YYDOC_MUTABLE 278
#define YYDOC_CHIFFRE 279
#define YYDOC_ENUM 280
#define YYDOC_STRUCT 281
#define YYDOC_CHAINE 282
#define YYDOC_FLECHE 283
#define YYDOC_FORALLEXISTS 284
#define YYDOC_TYPEDEF 285
#define YYDOC_CAST 286
#define YYDOC_POINTPOINT 287
#define YYDOC_NEW 288
#define YYDOC_EXTERN 289
#define YYDOC_EGALEGAL 290
#define YYDOC_DIFFEGAL 291
#define YYDOC_INFEGAL 292
#define YYDOC_SUPEGAL 293
#define YYDOC_OROR 294
#define YYDOC_ANDAND 295
#define UMINUS 296
#define NOT 297
#define COMP 298
#define TYPE 299




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;

