lcomment = #

kw1 = .DEFAULT, .DELETE_ON_ERROR, .EXPORT_ALL_VARIABLES, .IGNORE, .INTERMEDIATE
kw1 += .PHONY, .POSIX, .PRECIOUS, .SECONDARY, .SILENT, .SUFFIXES,

kw2 = define, endef, ifdef, ifndef, ifeq, ifneq, include, endif, else
kw2 += override, export, unexport, vpath
kw2 += subst, patsubst, strip, findstring, filter, filter-out
kw2 += sort, dir, notdir, suffix, basename, addsuffix, addprefix
kw2 += join, word, words, firstword, shell,  wildcard, origin, foreach
kw2 += $@, $%, $<, $?, $^, $+, $*
kw2 += $(@D), $(@F)
kw2 += $(*D), $(*F)
kw2 += $(%D), $(%F)
kw2 += $(<D), $(<F)
kw2 += $(^D), $(^F)
kw2 += $(+D), $(+F)
kw2 += $(?D), $(?F)

kw3 = MAKEFILES, VPATH, SHELL, MAKESHELL, MAKE, MAKELEVEL, MAKEFLAGS
kw3 += MAKECMDGOALS, CURDIR, SUFFIXES


