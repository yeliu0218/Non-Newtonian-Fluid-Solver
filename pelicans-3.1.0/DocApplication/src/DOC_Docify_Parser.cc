/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

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

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



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




/* Copy the first part of user declarations.  */
#line 1 "Docify.bi"

#include <list>
#include <string>
#include <stdio.h>
#include <DOC_Argument.hh>
#include <DOC_Text.hh>
#include <DOC_Category.hh>
#include <DOC_Class.hh>
#include <DOC_Enum.hh>
#include <DOC_FriendDeclaration.hh>
#include <DOC_Function.hh>
#include <DOC_Sequence.hh>
#include <DOC_Method.hh>
#include <DOC_Struct.hh>
#include <DOC_Tools.hh>
#include <DOC_Type.hh>
#include <DOC_Typedef.hh>
#include <PEL_assertions.hh>
#include <PEL_Root.hh>

extern int yylex( void ) ;
extern void yyerror( const char * s )  ;
extern void switch_mode( void ) ;

DOC_Attribute::Protection  protection_courante = DOC_Attribute::Private ;
DOC_Category const* category_courante = 0 ;
DOC_Text * comment_classe = 0 ;
DOC_Text * comment_classe_old = 0 ;
DOC_Text * comment_courant = 0 ;
DOC_Method * methode_courante = 0 ;


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 226 "Docify.bi.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  3
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   703

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  65
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  46
/* YYNRULES -- Number of rules.  */
#define YYNRULES  144
/* YYNRULES -- Number of states.  */
#define YYNSTATES  285

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   299

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    64,     2,     2,     2,    45,    47,     2,
      61,    62,    43,    41,    59,    42,    54,    44,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    58,    55,
      48,    60,    49,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    56,    46,    57,    63,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    50,    51,    52,    53
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     6,     9,    11,    13,    15,    17,
      19,    21,    23,    25,    27,    31,    38,    48,    50,    51,
      54,    56,    58,    61,    64,    66,    68,    70,    73,    75,
      77,    79,    85,    86,    89,    93,    95,    97,   102,   109,
     115,   118,   123,   124,   127,   129,   131,   133,   139,   146,
     149,   152,   157,   160,   163,   166,   167,   170,   171,   174,
     178,   180,   182,   186,   188,   191,   194,   197,   202,   206,
     208,   212,   214,   216,   218,   220,   223,   228,   229,   232,
     235,   237,   240,   242,   247,   250,   253,   255,   258,   261,
     264,   267,   272,   276,   278,   282,   286,   290,   294,   298,
     302,   306,   310,   314,   318,   322,   326,   330,   334,   338,
     342,   346,   349,   351,   353,   359,   367,   369,   373,   377,
     382,   384,   389,   394,   399,   404,   406,   411,   412,   416,
     419,   420,   422,   426,   433,   441,   442,   444,   448,   450,
     454,   456,   458,   459,   461
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      66,     0,    -1,    67,    -1,    -1,    67,    68,    -1,    78,
      -1,    69,    -1,    70,    -1,    71,    -1,    95,    -1,   101,
      -1,   102,    -1,     1,    -1,     4,    -1,    72,     6,    55,
      -1,    72,     6,    56,    73,    57,    55,    -1,    72,     6,
      58,     8,     6,    56,    73,    57,    55,    -1,     5,    -1,
      -1,    73,    74,    -1,    78,    -1,     7,    -1,    83,    58,
      -1,    84,    55,    -1,    85,    -1,    77,    -1,    80,    -1,
      75,    55,    -1,    71,    -1,    79,    -1,     1,    -1,    26,
      88,    56,    76,    57,    -1,    -1,    76,    85,    -1,    76,
      84,    55,    -1,    22,    -1,     3,    -1,    30,    89,    88,
      55,    -1,    25,    88,    56,    81,    57,    55,    -1,    25,
      56,    81,    57,    55,    -1,    88,    82,    -1,    81,    59,
      88,    82,    -1,    -1,    60,    97,    -1,     8,    -1,    10,
      -1,     9,    -1,    88,    61,    92,    62,    87,    -1,    89,
      88,    61,    92,    62,    87,    -1,    13,    84,    -1,    14,
      84,    -1,    89,     6,    86,    55,    -1,    14,    85,    -1,
      23,    85,    -1,    26,    85,    -1,    -1,    60,    97,    -1,
      -1,    87,    11,    -1,    87,    60,    24,    -1,     6,    -1,
      15,    -1,     6,    32,    88,    -1,    91,    -1,    89,    47,
      -1,    89,    43,    -1,    89,    11,    -1,    88,    48,    90,
      49,    -1,    89,    32,    88,    -1,    89,    -1,    90,    59,
      89,    -1,    12,    -1,    16,    -1,    88,    -1,    12,    -1,
      94,    93,    -1,    92,    59,    94,    93,    -1,    -1,    60,
      97,    -1,    89,    88,    -1,    84,    -1,   100,    55,    -1,
      96,    -1,    89,    88,    93,    55,    -1,    14,    96,    -1,
      11,    96,    -1,    24,    -1,    63,    97,    -1,    43,    97,
      -1,    47,    97,    -1,    42,    97,    -1,    61,    91,    62,
      97,    -1,    61,    97,    62,    -1,    27,    -1,    41,    41,
      97,    -1,    97,    41,    41,    -1,    97,    41,    97,    -1,
      97,    43,    97,    -1,    97,    44,    97,    -1,    97,    45,
      97,    -1,    97,    42,    97,    -1,    97,    35,    97,    -1,
      97,    36,    97,    -1,    97,    39,    97,    -1,    97,    40,
      97,    -1,    97,    47,    97,    -1,    97,    46,    97,    -1,
      97,    48,    97,    -1,    97,    49,    97,    -1,    97,    37,
      97,    -1,    97,    38,    97,    -1,    64,    97,    -1,   105,
      -1,    99,    -1,    33,    89,    61,   104,    62,    -1,    31,
      48,    89,    49,    61,    97,    62,    -1,    88,    -1,    99,
      28,    99,    -1,    99,    54,    99,    -1,    99,    61,   104,
      62,    -1,    98,    -1,    18,    61,    97,    62,    -1,    19,
      61,    97,    62,    -1,    17,    61,    97,    62,    -1,    20,
      61,    27,    62,    -1,    21,    -1,    34,    56,   103,    57,
      -1,    -1,   103,    84,    55,    -1,   103,    69,    -1,    -1,
      97,    -1,   104,    59,    97,    -1,    29,    61,   106,    59,
      97,    62,    -1,    61,   107,    55,    97,    55,   110,    62,
      -1,    -1,   108,    -1,   107,    59,   108,    -1,   109,    -1,
     109,    60,    97,    -1,    99,    -1,    94,    -1,    -1,    97,
      -1,   110,    59,    97,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,    88,    88,    90,    91,    93,    94,    95,    96,    97,
      98,    99,   100,   105,   111,   113,   122,   143,   154,   155,
     164,   165,   168,   169,   170,   171,   172,   173,   174,   175,
     176,   178,   189,   190,   194,   199,   208,   219,   229,   237,
     248,   252,   257,   258,   260,   263,   266,   270,   282,   294,
     298,   303,   313,   317,   321,   324,   325,   327,   328,   332,
     337,   338,   339,   348,   349,   350,   351,   352,   355,   359,
     360,   363,   364,   366,   369,   370,   375,   382,   383,   385,
     389,   393,   394,   396,   397,   398,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   411,   413,   415,   417,   419,
     421,   423,   425,   427,   429,   431,   433,   435,   437,   439,
     441,   443,   444,   445,   446,   448,   452,   453,   455,   457,
     459,   461,   464,   467,   472,   479,   480,   482,   483,   484,
     486,   487,   490,   494,   497,   500,   501,   502,   507,   508,
     514,   514,   516,   517,   520
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "YYDOC_COMMENT", "YYDOC_INCLUDE",
  "YYDOC_CLASS", "YYDOC_IDENTIF", "YYDOC_CATEGORY", "YYDOC_PUBLIC",
  "YYDOC_PROTECTED", "YYDOC_PRIVATE", "YYDOC_CONST", "YYDOC_VOID",
  "YYDOC_VIRTUAL", "YYDOC_STATIC", "YYDOC_OPERATOR", "YYDOC_TYPE_BASE",
  "YYDOC_PEL_ASSERT", "YYDOC_PEL_CHECK_PRE", "YYDOC_PEL_CHECK_POST",
  "YYDOC_PEL_LABEL", "YYDOC_USING", "YYDOC_FRIEND", "YYDOC_MUTABLE",
  "YYDOC_CHIFFRE", "YYDOC_ENUM", "YYDOC_STRUCT", "YYDOC_CHAINE",
  "YYDOC_FLECHE", "YYDOC_FORALLEXISTS", "YYDOC_TYPEDEF", "YYDOC_CAST",
  "YYDOC_POINTPOINT", "YYDOC_NEW", "YYDOC_EXTERN", "YYDOC_EGALEGAL",
  "YYDOC_DIFFEGAL", "YYDOC_INFEGAL", "YYDOC_SUPEGAL", "YYDOC_OROR",
  "YYDOC_ANDAND", "'+'", "'-'", "'*'", "'/'", "'%'", "'|'", "'&'", "'<'",
  "'>'", "UMINUS", "NOT", "COMP", "TYPE", "'.'", "';'", "'{'", "'}'",
  "':'", "','", "'='", "'('", "')'", "'~'", "'!'", "$accept", "file",
  "liste_elements", "element", "directive_include", "class_predeclaration",
  "class_declaration", "a__classe", "liste_directives_classe",
  "directives_classe", "struct_dec", "liste_attribute", "class_friend",
  "comment", "type_def", "enum", "liste_identificateurs",
  "initialiseur_enum", "niveau_protection", "prototype", "attribute",
  "defaut_attribute", "liste_modificateur", "identificateur", "type",
  "liste_type", "type_simple", "liste_arguments", "val_defaut", "argument",
  "method_implementation", "var_decl", "expression", "cast", "identifiant",
  "instruction", "using", "extern_def", "liste_prototype", "liste_args",
  "forall_exp", "for_desc", "for_init", "for_init_item",
  "for_init_identif", "expression_for", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,    43,    45,    42,    47,    37,   124,    38,    60,    62,
     296,   297,   298,   299,    46,    59,   123,   125,    58,    44,
      61,    40,    41,   126,    33
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    65,    66,    67,    67,    68,    68,    68,    68,    68,
      68,    68,    68,    69,    70,    71,    71,    72,    73,    73,
      74,    74,    74,    74,    74,    74,    74,    74,    74,    74,
      74,    75,    76,    76,    76,    77,    78,    79,    80,    80,
      81,    81,    82,    82,    83,    83,    83,    84,    84,    84,
      84,    85,    85,    85,    85,    86,    86,    87,    87,    87,
      88,    88,    88,    89,    89,    89,    89,    89,    89,    90,
      90,    91,    91,    91,    92,    92,    92,    93,    93,    94,
      95,    95,    95,    96,    96,    96,    97,    97,    97,    97,
      97,    97,    97,    97,    97,    97,    97,    97,    97,    97,
      97,    97,    97,    97,    97,    97,    97,    97,    97,    97,
      97,    97,    97,    97,    97,    98,    99,    99,    99,    99,
      99,   100,   100,   100,   100,   101,   102,   103,   103,   103,
     104,   104,   104,   105,   106,   107,   107,   107,   108,   108,
     109,   109,   110,   110,   110
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     0,     2,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     3,     6,     9,     1,     0,     2,
       1,     1,     2,     2,     1,     1,     1,     2,     1,     1,
       1,     5,     0,     2,     3,     1,     1,     4,     6,     5,
       2,     4,     0,     2,     1,     1,     1,     5,     6,     2,
       2,     4,     2,     2,     2,     0,     2,     0,     2,     3,
       1,     1,     3,     1,     2,     2,     2,     4,     3,     1,
       3,     1,     1,     1,     1,     2,     4,     0,     2,     2,
       1,     2,     1,     4,     2,     2,     1,     2,     2,     2,
       2,     4,     3,     1,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     2,     1,     1,     5,     7,     1,     3,     3,     4,
       1,     4,     4,     4,     4,     1,     4,     0,     3,     2,
       0,     1,     3,     6,     7,     0,     1,     3,     1,     3,
       1,     1,     0,     1,     3
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       3,     0,     0,     1,    12,    36,    13,    17,    60,     0,
      71,     0,     0,    61,    72,     0,     0,     0,     0,   125,
       0,     4,     6,     7,     8,     0,     5,    80,    73,     0,
      63,     9,    82,     0,    10,    11,     0,     0,    73,     0,
      85,     0,    49,     0,    50,    84,     0,     0,     0,     0,
     127,     0,     0,     0,    66,     0,    65,    64,    77,    81,
      62,    77,     0,    86,    93,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   116,     0,   120,   113,   112,
       0,     0,     0,     0,    14,    18,     0,    69,     0,    71,
       0,     0,    77,    68,     0,     0,     0,     0,     0,     0,
       0,    90,    88,    89,   116,     0,     0,    87,   111,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   123,     0,     0,   130,   121,   122,
     124,   126,   129,     0,     0,     0,    67,     0,    79,     0,
      57,    75,    78,     0,    83,   135,     0,     0,   130,    94,
       0,    92,   101,   102,   109,   110,   103,   104,    95,    96,
     100,    97,    98,    99,   106,   105,   107,   108,   117,   118,
     131,     0,   128,    30,    21,    44,    46,    45,     0,    35,
       0,     0,     0,     0,     0,    28,     0,    19,     0,    25,
      20,    29,    26,     0,     0,    24,     0,     0,    70,    77,
      47,    57,    73,   141,   140,     0,   136,   138,     0,     0,
       0,    91,     0,   119,     0,    52,     0,    53,     0,     0,
       0,    54,    73,     0,    15,     0,    27,    22,    23,    55,
      18,    76,    58,     0,    48,     0,     0,     0,     0,     0,
     114,   132,    55,     0,    42,     0,    32,     0,     0,     0,
       0,    59,     0,   137,   139,   133,     0,     0,     0,     0,
      40,     0,     0,    37,    56,    51,     0,   142,   115,    39,
      42,    43,     0,    31,     0,    33,    16,   143,     0,    41,
      38,    34,     0,   134,   144
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,     2,    21,    22,    23,   185,   186,   134,   187,
     188,   262,   189,   190,   191,   192,   243,   260,   193,    44,
     195,   249,   200,    75,    90,    88,    30,    91,    96,    92,
      31,    45,   170,    77,    78,    33,    34,    35,    83,   171,
      79,   146,   205,   206,   207,   278
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -143
static const yytype_int16 yypact[] =
{
    -143,    48,   409,  -143,  -143,  -143,  -143,  -143,    47,   290,
    -143,   324,   195,  -143,  -143,    25,    30,    49,    69,  -143,
      43,  -143,  -143,  -143,  -143,   105,  -143,  -143,   -40,    62,
    -143,  -143,  -143,    63,  -143,  -143,    61,   290,    94,    62,
    -143,   324,  -143,    62,  -143,  -143,    98,    98,    98,   119,
    -143,    68,    86,   282,  -143,    61,  -143,  -143,   -43,  -143,
    -143,    90,    97,  -143,  -143,   102,   121,    86,   126,    98,
      98,    98,   141,    98,    98,  -143,   356,  -143,   -12,  -143,
     448,   476,   109,   223,  -143,  -143,   167,    13,   -29,    28,
      62,    57,    90,  -143,    98,   282,   118,   116,    86,    14,
      98,  -143,  -143,  -143,   123,   129,   504,  -143,  -143,    98,
      98,    98,    98,    98,    98,   317,    98,    98,    98,    98,
      98,    98,    98,    98,  -143,    16,    16,    98,  -143,  -143,
    -143,  -143,  -143,   137,   360,   190,  -143,    86,  -143,    86,
    -143,  -143,   624,    93,  -143,   243,   139,   117,    98,   639,
      98,  -143,  -143,  -143,    -3,    -3,   216,   216,  -143,   639,
     639,   654,   654,   654,   216,   216,    -3,    -3,   142,   142,
     624,   127,  -143,  -143,  -143,  -143,  -143,  -143,   632,  -143,
     319,    -1,   319,    86,   161,  -143,   217,  -143,   173,  -143,
    -143,  -143,  -143,   164,   175,  -143,    74,   176,    13,    90,
       2,  -143,   165,  -143,   -12,    23,  -143,   180,    98,   181,
     128,  -143,    98,  -143,   319,  -143,   319,  -143,   101,    61,
     185,  -143,    15,    62,  -143,     8,  -143,  -143,  -143,   -17,
    -143,  -143,  -143,   220,     2,    98,   243,    98,   532,    98,
    -143,   624,   197,    38,   201,    61,  -143,   208,    98,   211,
     447,  -143,   588,  -143,   624,  -143,   560,   221,    61,    98,
    -143,    79,   269,  -143,   624,  -143,   222,    98,  -143,  -143,
     201,   624,   231,  -143,   232,  -143,  -143,   624,   156,  -143,
    -143,  -143,    98,  -143,   624
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -143,  -143,  -143,  -143,   206,  -143,   276,   288,    70,  -143,
    -143,  -143,  -143,   289,  -143,  -143,    58,    37,  -143,     0,
    -142,  -143,   146,    -2,    17,  -143,   252,   233,   -91,  -133,
    -143,    50,   199,  -143,  -122,  -143,  -143,  -143,  -143,   179,
    -143,  -143,  -143,   107,  -143,  -143
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -117
static const yytype_int16 yytable[] =
{
      28,   141,    27,   168,   169,     8,   199,    38,    52,    28,
      28,    42,   203,   232,    13,    36,   125,    94,    95,    29,
     136,    53,     8,   204,    54,    54,    39,    58,    43,    29,
     137,    13,   109,   110,    60,    38,   215,    61,   217,    28,
     221,    62,   126,   248,   -60,    55,    55,    66,     3,   127,
      38,    38,    32,    93,    39,   219,    56,    56,    43,    40,
      57,    57,   233,    52,    85,    38,    86,     8,     8,    87,
     104,   246,   221,    54,   215,   148,    13,    13,   235,    36,
     229,    28,   236,   133,    99,    54,    46,   -74,   138,    13,
     -74,    47,     8,    38,    55,   257,    38,   258,    10,    50,
      43,    13,    14,   203,     8,    56,    55,   242,   231,    57,
      48,    51,    54,    13,   204,   147,   139,    56,    59,   140,
     275,    57,    63,    84,    85,    64,    86,    65,    54,    66,
      49,    67,    28,    55,   194,    38,   272,    38,   258,    68,
      69,    70,    52,   202,    56,    71,    82,     8,    57,    55,
      94,   196,   139,    10,   198,   201,    13,    14,    95,    72,
      56,    73,    74,    97,    57,    63,   209,   100,    64,    98,
      65,   130,    66,   144,    67,   135,    28,   145,    38,   220,
     222,    38,    68,    69,    70,   -73,   212,   212,    71,   213,
     240,   150,   172,  -116,    62,   196,   197,   218,   208,   218,
     223,     8,    72,   127,    73,    74,     9,    10,    11,    12,
      13,    14,    38,    52,    38,   282,   224,   244,   283,  -116,
    -116,   247,   227,   225,  -116,  -116,  -116,     6,   226,     8,
     228,   218,   230,   218,   202,    10,    11,    41,    13,    14,
     237,   245,   239,   244,   251,    76,    80,    81,    28,     8,
     194,   109,   110,   111,   112,    10,   270,   248,    13,    14,
      28,   259,   274,   263,   122,   123,   265,   196,   101,   102,
     103,   106,   107,   108,    66,     8,   269,   276,    24,   196,
     131,    10,    11,   178,    13,    14,   280,   281,     8,   132,
      25,    26,   180,   142,    89,   214,     8,    13,    14,   149,
     250,     9,    10,   261,    37,    13,    14,   279,   152,   153,
     154,   155,   156,   157,   159,   160,   161,   162,   163,   164,
     165,   166,   167,     8,   105,     8,   273,   210,   143,     0,
       8,    10,    13,   216,    13,    14,    10,    11,    41,    13,
      14,    63,   180,   253,    64,   214,    65,   234,    66,   211,
      67,     0,     0,     0,     0,     0,     0,     0,   158,    69,
      70,   173,     0,     5,    71,     7,     8,   174,   175,   176,
     177,     0,    10,    11,   178,    13,    14,     0,    72,     0,
      73,    74,   179,   180,     0,   181,   182,     0,     0,     0,
     183,   109,   110,   111,   112,   113,   114,   115,   116,   117,
     118,   119,   120,   121,   122,   123,     0,   238,     0,    -2,
       4,   241,     5,     6,     7,     8,     0,   184,   124,     0,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,     0,     0,     0,   252,     0,   254,     0,   256,     0,
       0,     0,     0,    20,     0,     0,     0,   264,   173,     0,
       5,     0,     7,     8,   174,   175,   176,   177,   271,    10,
      11,   178,    13,    14,     0,     0,   277,     0,     0,   179,
     180,     0,   181,   182,     0,     0,     0,   183,     0,     0,
       0,   284,     0,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,   119,   120,   121,   122,   123,     0,     0,
       0,     0,     0,     0,   266,     0,     0,     0,     0,     0,
     128,   109,   110,   111,   112,   113,   114,   115,   116,   117,
     118,   119,   120,   121,   122,   123,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   129,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     120,   121,   122,   123,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   151,   109,   110,   111,
     112,   113,   114,   115,   116,   117,   118,   119,   120,   121,
     122,   123,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   255,   109,   110,   111,   112,   113,
     114,   115,   116,   117,   118,   119,   120,   121,   122,   123,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   268,   109,   110,   111,   112,   113,   114,   115,
     116,   117,   118,   119,   120,   121,   122,   123,     8,     0,
       0,     0,     0,   267,    10,    11,   178,    13,    14,     0,
       0,     0,     0,     0,     0,   180,     0,     0,   214,   109,
     110,   111,   112,   113,   114,   115,   116,   117,   118,   119,
     120,   121,   122,   123,   109,   110,   111,   112,   113,   114,
       0,     0,   117,   118,   119,   120,   121,   122,   123,   109,
     110,   111,   112,   113,   114,     0,     0,     0,     0,     0,
     120,   121,   122,   123
};

static const yytype_int16 yycheck[] =
{
       2,    92,     2,   125,   126,     6,   139,     9,    48,    11,
      12,    11,   145,    11,    15,    32,    28,    60,    61,     2,
      49,    61,     6,   145,    11,    11,     9,    29,    11,    12,
      59,    15,    35,    36,    36,    37,   178,    39,   180,    41,
     182,    43,    54,    60,    61,    32,    32,    31,     0,    61,
      52,    53,     2,    55,    37,    56,    43,    43,    41,     9,
      47,    47,    60,    48,    56,    67,    58,     6,     6,    52,
      72,    56,   214,    11,   216,    61,    15,    15,    55,    32,
       6,    83,    59,    83,    67,    11,    61,    59,    90,    15,
      62,    61,     6,    95,    32,    57,    98,    59,    12,    56,
      83,    15,    16,   236,     6,    43,    32,     6,   199,    47,
      61,     6,    11,    15,   236,    98,    59,    43,    55,    62,
     262,    47,    24,    55,    56,    27,    58,    29,    11,    31,
      61,    33,   134,    32,   134,   137,    57,   139,    59,    41,
      42,    43,    48,   145,    43,    47,    27,     6,    47,    32,
      60,   134,    59,    12,   137,    62,    15,    16,    61,    61,
      43,    63,    64,    61,    47,    24,    49,    41,    27,    48,
      29,    62,    31,    55,    33,     8,   178,    61,   180,   181,
     182,   183,    41,    42,    43,    62,    59,    59,    47,    62,
      62,    62,    55,    28,   196,   178,     6,   180,    59,   182,
     183,     6,    61,    61,    63,    64,    11,    12,    13,    14,
      15,    16,   214,    48,   216,    59,    55,   219,    62,    54,
      55,   223,    58,     6,    59,    60,    61,     4,    55,     6,
      55,   214,    56,   216,   236,    12,    13,    14,    15,    16,
      60,    56,    61,   245,    24,    46,    47,    48,   250,     6,
     250,    35,    36,    37,    38,    12,   258,    60,    15,    16,
     262,    60,   262,    55,    48,    49,    55,   250,    69,    70,
      71,    72,    73,    74,    31,     6,    55,    55,     2,   262,
      57,    12,    13,    14,    15,    16,    55,    55,     6,    83,
       2,     2,    23,    94,    12,    26,     6,    15,    16,   100,
     230,    11,    12,   245,    14,    15,    16,   270,   109,   110,
     111,   112,   113,   114,   115,   116,   117,   118,   119,   120,
     121,   122,   123,     6,    72,     6,    57,   148,    95,    -1,
       6,    12,    15,    14,    15,    16,    12,    13,    14,    15,
      16,    24,    23,   236,    27,    26,    29,   201,    31,   150,
      33,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    41,    42,
      43,     1,    -1,     3,    47,     5,     6,     7,     8,     9,
      10,    -1,    12,    13,    14,    15,    16,    -1,    61,    -1,
      63,    64,    22,    23,    -1,    25,    26,    -1,    -1,    -1,
      30,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    -1,   208,    -1,     0,
       1,   212,     3,     4,     5,     6,    -1,    57,    62,    -1,
      11,    12,    13,    14,    15,    16,    17,    18,    19,    20,
      21,    -1,    -1,    -1,   235,    -1,   237,    -1,   239,    -1,
      -1,    -1,    -1,    34,    -1,    -1,    -1,   248,     1,    -1,
       3,    -1,     5,     6,     7,     8,     9,    10,   259,    12,
      13,    14,    15,    16,    -1,    -1,   267,    -1,    -1,    22,
      23,    -1,    25,    26,    -1,    -1,    -1,    30,    -1,    -1,
      -1,   282,    -1,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    -1,    -1,
      -1,    -1,    -1,    -1,    57,    -1,    -1,    -1,    -1,    -1,
      62,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    62,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    62,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      48,    49,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    62,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,    48,    49,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    62,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,     6,    -1,
      -1,    -1,    -1,    55,    12,    13,    14,    15,    16,    -1,
      -1,    -1,    -1,    -1,    -1,    23,    -1,    -1,    26,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    35,    36,    37,    38,    39,    40,
      -1,    -1,    43,    44,    45,    46,    47,    48,    49,    35,
      36,    37,    38,    39,    40,    -1,    -1,    -1,    -1,    -1,
      46,    47,    48,    49
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    66,    67,     0,     1,     3,     4,     5,     6,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      34,    68,    69,    70,    71,    72,    78,    84,    88,    89,
      91,    95,    96,   100,   101,   102,    32,    14,    88,    89,
      96,    14,    84,    89,    84,    96,    61,    61,    61,    61,
      56,     6,    48,    61,    11,    32,    43,    47,    88,    55,
      88,    88,    88,    24,    27,    29,    31,    33,    41,    42,
      43,    47,    61,    63,    64,    88,    97,    98,    99,   105,
      97,    97,    27,   103,    55,    56,    58,    89,    90,    12,
      89,    92,    94,    88,    60,    61,    93,    61,    48,    89,
      41,    97,    97,    97,    88,    91,    97,    97,    97,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    45,
      46,    47,    48,    49,    62,    28,    54,    61,    62,    62,
      62,    57,    69,    84,    73,     8,    49,    59,    88,    59,
      62,    93,    97,    92,    55,    61,   106,    89,    61,    97,
      62,    62,    97,    97,    97,    97,    97,    97,    41,    97,
      97,    97,    97,    97,    97,    97,    97,    97,    99,    99,
      97,   104,    55,     1,     7,     8,     9,    10,    14,    22,
      23,    25,    26,    30,    57,    71,    72,    74,    75,    77,
      78,    79,    80,    83,    84,    85,    89,     6,    89,    94,
      87,    62,    88,    94,    99,   107,   108,   109,    59,    49,
     104,    97,    59,    62,    26,    85,    14,    85,    89,    56,
      88,    85,    88,    89,    55,     6,    55,    58,    55,     6,
      56,    93,    11,    60,    87,    55,    59,    60,    97,    61,
      62,    97,     6,    81,    88,    56,    56,    88,    60,    86,
      73,    24,    97,   108,    97,    62,    97,    57,    59,    60,
      82,    81,    76,    55,    97,    55,    57,    55,    62,    55,
      88,    97,    57,    57,    84,    85,    55,    97,   110,    82,
      55,    55,    59,    62,    97
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 88 "Docify.bi"
    { comment_classe = comment_courant = 0 ; ;}
    break;

  case 4:
#line 91 "Docify.bi"
    { protection_courante = DOC_Attribute::Private ; ;}
    break;

  case 8:
#line 96 "Docify.bi"
    { comment_classe_old =  0 ; ;}
    break;

  case 12:
#line 101 "Docify.bi"
    {
             yyerror( "unknown element" ) ;
          ;}
    break;

  case 13:
#line 105 "Docify.bi"
    { comment_classe = 0 ;
				  comment_courant = 0 ;
                                  DOC_Text const* txt =  (yyvsp[(1) - (1)])->to_text() ;
				  if( DOC_Tools::read( DOC_Tools::file_include( txt->text() ) ) ) 
					YYACCEPT ; ;}
    break;

  case 15:
#line 113 "Docify.bi"
    {
	   string const& name = (yyvsp[(2) - (6)])->to_text()->text() ;
	   DOC_Class::create(   name, 
			        comment_classe, 
			     	DOC_Tools::file(), 
				(yyvsp[(2) - (6)])->to_text()->line_number(),	      
				(yyvsp[(4) - (6)])->to_sequence() ) ;
	   comment_classe = 0 ;
	   (yyval) = 0 ; ;}
    break;

  case 16:
#line 122 "Docify.bi"
    {
	   string const& name = (yyvsp[(2) - (9)])->to_text()->text() ;
	   string const& name_mother = (yyvsp[(5) - (9)])->to_text()->text() ;
	   DOC_Class * mother = DOC_Class::search( name_mother ) ;
	   if( mother==0 ) mother = 
	     DOC_Class::create( name_mother, 
				0, 
			     	"", 
				0,	      
				0,
				0 ) ;
	   DOC_Class::create( name, 
				comment_classe, 
			     	DOC_Tools::file(), 
				(yyvsp[(2) - (9)])->to_text()->line_number(),	      
				(yyvsp[(7) - (9)])->to_sequence(),
				mother ) ;
	   
	   comment_classe = 0 ;
	   (yyval) = 0 ; ;}
    break;

  case 17:
#line 143 "Docify.bi"
    { 
	if( comment_classe==0 ) // Sinon, comment classes friend
	{
		comment_classe = comment_courant ; 
  		comment_courant=0 ; 
	}       
        else if( comment_classe_old==0 ) 
               comment_classe_old = comment_classe ;

	(yyval) = (yyvsp[(1) - (1)]) ; ;}
    break;

  case 18:
#line 154 "Docify.bi"
    { (yyval) = DOC_Sequence::create(0) ; ;}
    break;

  case 19:
#line 155 "Docify.bi"
    { 
	   if( (yyvsp[(2) - (2)]) !=0 ) 
	   {
	      (yyvsp[(1) - (2)])->to_sequence()->list()->append( (yyvsp[(2) - (2)])->to_class_item() ) ;
	      (yyvsp[(2) - (2)])->set_owner( (yyvsp[(1) - (2)]) ) ;
              (yyval) = (yyvsp[(1) - (2)]) ; 
           } ;}
    break;

  case 20:
#line 164 "Docify.bi"
    { (yyval) = 0 ; ;}
    break;

  case 21:
#line 165 "Docify.bi"
    { category_courante = DOC_Category::create( 
                                             (yyvsp[(1) - (1)])->to_text()->text() ) ; 
                                            (yyval) = 0 ; ;}
    break;

  case 22:
#line 168 "Docify.bi"
    { (yyval) = 0 ; ;}
    break;

  case 28:
#line 174 "Docify.bi"
    { (yyval) = 0 ; comment_classe = comment_classe_old ; ;}
    break;

  case 30:
#line 176 "Docify.bi"
    { (yyval) = 0 ; ;}
    break;

  case 31:
#line 178 "Docify.bi"
    {
                string const& id = (yyvsp[(2) - (5)])->to_text()->text() ;
		(yyval) = DOC_Struct::create( id, 
				  protection_courante,
				  category_courante,
				  comment_courant,
				  (yyvsp[(4) - (5)])->to_sequence()->list() ) ;

	        comment_courant = 0 ;}
    break;

  case 32:
#line 189 "Docify.bi"
    { (yyval) = DOC_Sequence::create(PEL_Root::object()) ; ;}
    break;

  case 33:
#line 191 "Docify.bi"
    { (yyvsp[(1) - (2)])->to_sequence()->list()->append( (yyvsp[(2) - (2)]) ) ; 
		    (yyvsp[(2) - (2)])->set_owner((yyvsp[(1) - (2)])) ;
                    (yyval)=(yyvsp[(1) - (2)]) ; ;}
    break;

  case 34:
#line 195 "Docify.bi"
    { (yyvsp[(1) - (3)])->to_sequence()->list()->append( (yyvsp[(2) - (3)]) ) ; 
                    (yyvsp[(2) - (3)])->set_owner((yyvsp[(1) - (3)])) ;
		    (yyval)=(yyvsp[(1) - (3)]) ; ;}
    break;

  case 35:
#line 199 "Docify.bi"
    { string const& str = (yyvsp[(1) - (1)])->to_text()->text() ; 
			       DOC_FriendDeclaration * fr =  DOC_FriendDeclaration::create( 
					str, 
					protection_courante,
					category_courante,
					comment_courant ) ;
				(yyval) = fr ;
				comment_courant = 0 ;}
    break;

  case 36:
#line 208 "Docify.bi"
    { if( comment_courant!=0 )
					   {
						comment_courant->append( 
						   DOC_Tools::text_comment( (yyvsp[(1) - (1)])->to_text()->text() ) ) ;
					   }
					   else
					   {
						comment_courant=DOC_Text::create( PEL_Root::object(),
						   DOC_Tools::text_comment( (yyvsp[(1) - (1)])->to_text()->text() ) ) ;
					   }
					   (yyval) = 0 ; ;}
    break;

  case 37:
#line 219 "Docify.bi"
    { 
		DOC_Type * typ = (yyvsp[(2) - (4)])->to_type() ; 
                string const& str = (yyvsp[(3) - (4)])->to_text()->text() ;
		(yyval) = DOC_Typedef::create( str, 
  			          typ,  
				  protection_courante,
				  category_courante,
				  comment_courant ) ;
	        comment_courant = 0 ; ;}
    break;

  case 38:
#line 230 "Docify.bi"
    {  string const& str = (yyvsp[(2) - (6)])->to_text()->text() ;
	 (yyval) = DOC_Enum::create( str, 
			protection_courante,
			category_courante,
			comment_courant,
                        (yyvsp[(4) - (6)])->to_sequence()->list() ) ;
         comment_courant = 0 ; ;}
    break;

  case 39:
#line 238 "Docify.bi"
    {  string const& str = "" ;
	 (yyval) = DOC_Enum::create( str, 
			protection_courante,
			category_courante,
			comment_courant,
                        (yyvsp[(3) - (5)])->to_sequence()->list() ) ;
         comment_courant = 0 ; ;}
    break;

  case 40:
#line 249 "Docify.bi"
    { DOC_Sequence * lst= DOC_Sequence::create( PEL_Root::object() ) ;
          lst->list()->append( (yyvsp[(1) - (2)]) ) ; 
	  (yyval) = lst ; ;}
    break;

  case 41:
#line 253 "Docify.bi"
    { DOC_Sequence * lst=(yyvsp[(1) - (4)])->to_sequence() ;
          lst->list()->append( (yyvsp[(3) - (4)]) ) ; 
	  (yyval) = lst ; ;}
    break;

  case 44:
#line 260 "Docify.bi"
    { category_courante = 0 ; 
				    protection_courante = DOC_Attribute::Public ; 
				    (yyval) = 0 ; ;}
    break;

  case 45:
#line 263 "Docify.bi"
    { category_courante = 0 ; 
				    protection_courante = DOC_Attribute::Private ; 
				    (yyval) = 0 ; ;}
    break;

  case 46:
#line 266 "Docify.bi"
    { category_courante = 0 ; 
				    protection_courante = DOC_Attribute::Protected ; 
				    (yyval) = 0 ; ;}
    break;

  case 47:
#line 270 "Docify.bi"
    { 
	   DOC_Method * meth =  DOC_Method::create( (yyvsp[(1) - (5)])->to_text()->text(), 
					    DOC_Type::create( PEL_Root::object(), "" ), 
					    protection_courante,
					    category_courante,
					    comment_courant,
					    (yyvsp[(3) - (5)])->to_sequence(),
					    (yyvsp[(5) - (5)])->to_sequence(),
					    (yyvsp[(1) - (5)])->to_text()->line_number(),
					    DOC_Tools::file() ) ; 
	   comment_courant = 0 ;
	   (yyval) = meth ; ;}
    break;

  case 48:
#line 282 "Docify.bi"
    { 
	   DOC_Method * meth = DOC_Method::create( (yyvsp[(2) - (6)])->to_text()->text(), 
					 (yyvsp[(1) - (6)])->to_type(), 
					 protection_courante,
					 category_courante,
 					 comment_courant,
					 (yyvsp[(4) - (6)])->to_sequence(),
					 (yyvsp[(6) - (6)])->to_sequence(),
					 (yyvsp[(2) - (6)])->to_text()->line_number(),
					 DOC_Tools::file() ) ; 
	   comment_courant = 0 ;
	   (yyval) = meth ; ;}
    break;

  case 49:
#line 294 "Docify.bi"
    { DOC_Method * meth = (yyvsp[(2) - (2)])->to_class_item()->method() ; 
					  PEL_ASSERT( meth!=0 ) ;
                                          meth->set_virtual() ; 
                                          (yyval)=meth ; ;}
    break;

  case 50:
#line 298 "Docify.bi"
    { DOC_Method * meth = (yyvsp[(2) - (2)])->to_class_item()->method() ; 
					  PEL_ASSERT( meth!=0 ) ;
                                          meth->set_static() ; 
                                          (yyval)=meth ; ;}
    break;

  case 51:
#line 303 "Docify.bi"
    { DOC_Type * typ = (yyvsp[(1) - (4)])->to_type() ; 
                                  string const& str = (yyvsp[(2) - (4)])->to_text()->text() ; 
				  DOC_Attribute * att =  DOC_Attribute::create( str, 
						     typ,  
						     protection_courante,
						     category_courante,
						     comment_courant ) ;
				  if( (yyvsp[(3) - (4)])!=0 ) att->initialize( (yyvsp[(3) - (4)]) ) ;
				  (yyval) = att ;
				  comment_courant = 0 ;}
    break;

  case 52:
#line 313 "Docify.bi"
    { DOC_Attribute * att = (yyvsp[(2) - (2)])->to_class_item()->attribute() ;
				 PEL_ASSERT( att!=0 ) ;
                                 att->set_static() ;
				 (yyval) = att ; ;}
    break;

  case 53:
#line 317 "Docify.bi"
    { DOC_Attribute * att = (yyvsp[(2) - (2)])->to_class_item()->attribute() ;
				 PEL_ASSERT( att!=0 ) ;
                                 // att->mutable() ;
				 (yyval) = att ; ;}
    break;

  case 54:
#line 321 "Docify.bi"
    { DOC_Attribute * att = (yyvsp[(2) - (2)])->to_class_item()->attribute() ;
				 PEL_ASSERT( att!=0 ) ;
				 (yyval) = att ; ;}
    break;

  case 55:
#line 324 "Docify.bi"
    { (yyval) = 0 ; ;}
    break;

  case 56:
#line 325 "Docify.bi"
    { (yyval) = (yyvsp[(2) - (2)]) ; ;}
    break;

  case 57:
#line 327 "Docify.bi"
    { (yyval) = DOC_Sequence::create( PEL_Root::object() ) ; ;}
    break;

  case 58:
#line 328 "Docify.bi"
    { 
   				 DOC_Sequence * lst = (yyvsp[(1) - (2)])->to_sequence() ;
				 lst->list()->append( DOC_Text::create( PEL_Root::object(), "const" ) ) ;
				 (yyval) = lst ; ;}
    break;

  case 59:
#line 332 "Docify.bi"
    { 
   				 DOC_Sequence * lst = (yyvsp[(1) - (3)])->to_sequence() ;
				 lst->list()->append( DOC_Text::create( PEL_Root::object(), "abstract" ) ) ;
				 (yyval) = lst ; ;}
    break;

  case 62:
#line 339 "Docify.bi"
    { string const& str1 = (yyvsp[(1) - (3)])->to_text()->text() ;
					     string const& str2 = (yyvsp[(3) - (3)])->to_text()->text() ; 
					     (yyval) = DOC_Text::create( PEL_Root::object(),
								  str1+"::"+str2,
							          (yyvsp[(1) - (3)])->to_text()->line_number() ) ; ;}
    break;

  case 64:
#line 349 "Docify.bi"
    { DOC_Type * typ = (yyvsp[(1) - (2)])->to_type() ; typ->set_reference() ; (yyval) = typ ; ;}
    break;

  case 65:
#line 350 "Docify.bi"
    { DOC_Type * typ = (yyvsp[(1) - (2)])->to_type() ; typ->set_pointer() ; (yyval) = typ ; ;}
    break;

  case 66:
#line 351 "Docify.bi"
    { DOC_Type * typ = (yyvsp[(1) - (2)])->to_type() ; typ->set_constant() ; (yyval) = typ ; ;}
    break;

  case 67:
#line 352 "Docify.bi"
    { (yyval) = DOC_Type::create( PEL_Root::object(),
						 (yyvsp[(1) - (4)])->to_text()->text()+"<"+
                                                 (yyvsp[(3) - (4)])->to_text()->text()+">" )  ;}
    break;

  case 68:
#line 355 "Docify.bi"
    { (yyval) = DOC_Type::create( PEL_Root::object(),
						 (yyvsp[(1) - (3)])->to_type()->full_type_name()+"::"+
                                                 (yyvsp[(3) - (3)])->to_text()->text() )  ;}
    break;

  case 69:
#line 359 "Docify.bi"
    { (yyval) = DOC_Text::create( PEL_Root::object(), (yyvsp[(1) - (1)])->to_type()->full_type_name() ) ; ;}
    break;

  case 70:
#line 360 "Docify.bi"
    { (yyval) = DOC_Text::create( PEL_Root::object(), 
				      (yyvsp[(1) - (3)])->to_text()->text()+","+(yyvsp[(3) - (3)])->to_type()->full_type_name() ) ; ;}
    break;

  case 71:
#line 363 "Docify.bi"
    { (yyval) = DOC_Type::create( PEL_Root::object(), "void" ) ; ;}
    break;

  case 72:
#line 364 "Docify.bi"
    { string const& str = (yyvsp[(1) - (1)])->to_text()->text() ; 
				(yyval) = DOC_Type::create( PEL_Root::object(), str ) ; ;}
    break;

  case 73:
#line 366 "Docify.bi"
    { string const& str = (yyvsp[(1) - (1)])->to_text()->text() ; 
				(yyval) = DOC_Type::create( PEL_Root::object(), str ) ; ;}
    break;

  case 74:
#line 369 "Docify.bi"
    { (yyval) = DOC_Sequence::create( PEL_Root::object() ) ; ;}
    break;

  case 75:
#line 370 "Docify.bi"
    { DOC_Sequence* lst = DOC_Sequence::create( PEL_Root::object() ) ; 
				 DOC_Argument * arg = (yyvsp[(1) - (2)])->to_argument() ;
				 if( (yyvsp[(2) - (2)])!=0 ) arg->initialize( (yyvsp[(2) - (2)])->text() ) ;
                                 lst->list()->append( arg ) ; 
                                 (yyval) = lst ; ;}
    break;

  case 76:
#line 375 "Docify.bi"
    { DOC_Sequence* lst = (yyvsp[(1) - (4)])->to_sequence() ; 
				                    DOC_Argument * arg = (yyvsp[(3) - (4)])->to_argument() ;
				                    if( (yyvsp[(4) - (4)])!=0 ) 
							arg->initialize( (yyvsp[(4) - (4)])->text() ) ;
                                                    lst->list()->append( arg ) ; 
                                         	    (yyval) = lst ; ;}
    break;

  case 77:
#line 382 "Docify.bi"
    { (yyval) = 0 ; ;}
    break;

  case 78:
#line 383 "Docify.bi"
    { (yyval) = (yyvsp[(2) - (2)]) ; ;}
    break;

  case 79:
#line 385 "Docify.bi"
    { string const& str = (yyvsp[(2) - (2)])->to_text()->text() ; 
                                 DOC_Type * typ = (yyvsp[(1) - (2)])->to_type() ; 
                                 (yyval) = DOC_Argument::create( PEL_Root::object(), typ, str ) ; ;}
    break;

  case 80:
#line 389 "Docify.bi"
    {  if( (yyvsp[(1) - (1)])!=0 && (yyvsp[(1) - (1)])->owner()!=0 )
					    { methode_courante = (yyvsp[(1) - (1)])->to_class_item()->method() ;}
					 else
					    { methode_courante=0 ; if( (yyvsp[(1) - (1)])) (yyvsp[(1) - (1)])->destroy() ; } ;}
    break;

  case 81:
#line 393 "Docify.bi"
    { switch_mode() ; ;}
    break;

  case 87:
#line 402 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "~", (yyvsp[(2) - (2)]) ) ; ;}
    break;

  case 88:
#line 403 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "*", (yyvsp[(2) - (2)]) ) ; ;}
    break;

  case 89:
#line 404 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "&", (yyvsp[(2) - (2)]) ) ; ;}
    break;

  case 90:
#line 405 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "-", (yyvsp[(2) - (2)]) ) ; ;}
    break;

  case 91:
#line 406 "Docify.bi"
    { (yyval) = (yyvsp[(4) - (4)]) ; ;}
    break;

  case 92:
#line 407 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "", (yyvsp[(2) - (3)]) ) ; ;}
    break;

  case 94:
#line 410 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "++", (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 95:
#line 412 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "++", (yyvsp[(1) - (3)]) ) ; ;}
    break;

  case 96:
#line 414 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "+", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 97:
#line 416 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "*", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 98:
#line 418 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "/", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 99:
#line 420 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "%", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 100:
#line 422 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "-", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 101:
#line 424 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "==", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 102:
#line 426 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "!=", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 103:
#line 428 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "||", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 104:
#line 430 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "&&", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 105:
#line 432 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "&", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 106:
#line 434 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "|", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 107:
#line 436 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "<", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 108:
#line 438 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), ">", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 109:
#line 440 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "<=", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 110:
#line 442 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), ">=", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 111:
#line 443 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "!", (yyvsp[(2) - (2)]) ) ; ;}
    break;

  case 114:
#line 446 "Docify.bi"
    { (yyval) = 0 ; ;}
    break;

  case 115:
#line 448 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), (yyvsp[(1) - (7)])->to_text()->text(),
								     (yyvsp[(3) - (7)]),
								     (yyvsp[(6) - (7)]), false ) ; ;}
    break;

  case 116:
#line 452 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), (yyvsp[(1) - (1)])->to_text()->text() ) ; ;}
    break;

  case 117:
#line 454 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "->", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 118:
#line 456 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), ".", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ; ;}
    break;

  case 119:
#line 458 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), (yyvsp[(1) - (4)])->to_function(), (yyvsp[(3) - (4)])->to_sequence()->list() ) ; ;}
    break;

  case 121:
#line 462 "Docify.bi"
    { if( methode_courante ) 
                                       methode_courante->add_condition( (yyvsp[(3) - (4)])->to_function(), DOC_Method::pre ) ; ;}
    break;

  case 122:
#line 465 "Docify.bi"
    { if( methode_courante ) 
                                       methode_courante->add_condition( (yyvsp[(3) - (4)])->to_function(), DOC_Method::post ) ; ;}
    break;

  case 123:
#line 468 "Docify.bi"
    { if( methode_courante )
				   {
                                       methode_courante->add_condition( (yyvsp[(3) - (4)])->to_function(), DOC_Method::assert ) ; 
				    };}
    break;

  case 124:
#line 473 "Docify.bi"
    { if( methode_courante )
				   {
                                       methode_courante->set_label( (yyvsp[(3) - (4)])->to_text()->text() ) ; 
				    };}
    break;

  case 126:
#line 480 "Docify.bi"
    { (yyval)=0 ; ;}
    break;

  case 127:
#line 482 "Docify.bi"
    { (yyval) =0 ;;}
    break;

  case 128:
#line 483 "Docify.bi"
    { if((yyvsp[(2) - (3)])) (yyvsp[(2) - (3)])->destroy() ; (yyval)=0 ; ;}
    break;

  case 129:
#line 484 "Docify.bi"
    { (yyval) =0;;}
    break;

  case 130:
#line 486 "Docify.bi"
    { (yyval) = DOC_Sequence::create( PEL_Root::object() ) ; ;}
    break;

  case 131:
#line 487 "Docify.bi"
    { DOC_Sequence * lst = DOC_Sequence::create( PEL_Root::object() ) ; 
					lst->list()->append( (yyvsp[(1) - (1)]) ) ; 
					(yyval) = lst ; ;}
    break;

  case 132:
#line 490 "Docify.bi"
    { DOC_Sequence * lst = (yyvsp[(1) - (3)])->to_sequence() ; 
					lst->list()->append( (yyvsp[(3) - (3)]) ) ; 
					(yyval) = lst ; ;}
    break;

  case 133:
#line 494 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), (yyvsp[(1) - (6)])->to_text()->text(),
							      (yyvsp[(3) - (6)]), (yyvsp[(5) - (6)]), false ) ; ;}
    break;

  case 134:
#line 497 "Docify.bi"
    { (yyval) = DOC_Function::create( PEL_Root::object(), "",
							      (yyvsp[(2) - (7)]), (yyvsp[(4) - (7)]), (yyvsp[(6) - (7)]) ) ; ;}
    break;

  case 135:
#line 500 "Docify.bi"
    { (yyval) = DOC_Sequence::create( PEL_Root::object() ) ; ;}
    break;

  case 137:
#line 502 "Docify.bi"
    {
			DOC_Sequence * lst = (yyvsp[(1) - (3)])->to_sequence() ;
			lst->list()->append( (yyvsp[(3) - (3)]) ) ;
			(yyval) = lst ; ;}
    break;

  case 139:
#line 508 "Docify.bi"
    {
			DOC_Sequence * lst = DOC_Sequence::create( PEL_Root::object() ) ;
			DOC_Function * fonc = DOC_Function::create( PEL_Root::object(), "=", (yyvsp[(1) - (3)]), (yyvsp[(3) - (3)]) ) ;
			lst->list()->append( fonc ) ;
			(yyval) = lst ; ;}
    break;

  case 142:
#line 516 "Docify.bi"
    { (yyval) = DOC_Sequence::create( PEL_Root::object() ) ; ;}
    break;

  case 143:
#line 517 "Docify.bi"
    {  DOC_Sequence * lst = DOC_Sequence::create( PEL_Root::object() ) ;
			lst->list()->append( (yyvsp[(1) - (1)]) ) ;
			(yyval) = lst ; ;}
    break;

  case 144:
#line 520 "Docify.bi"
    {
			DOC_Sequence * lst = (yyvsp[(1) - (3)])->to_sequence() ;
			lst->list()->append( (yyvsp[(3) - (3)]) ) ;
			(yyval) = lst ; ;}
    break;


/* Line 1267 of yacc.c.  */
#line 2534 "Docify.bi.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



