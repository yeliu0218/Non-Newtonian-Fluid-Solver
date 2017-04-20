;;
;;  Copyright 1995-2010 by IRSN
;;
;;  This software is an application framework, with a set of integrated  
;;  reusable components, whose purpose is to simplify the task of developing 
;;  softwares of numerical mathematics and scientific computing.
;; 
;;  This software is governed by the CeCILL-C license under French law and 
;;  abiding by the rules of distribution of free software. You can use, modify 
;;  and/or redistribute the software under the terms of the CeCILL-C license  
;;  as circulated by CEA, CNRS and INRIA at the following URL 
;;  "http://www.cecill.info". 
;;
;;  As a counterpart to the access to the source code and rights to copy,  
;;  modify and redistribute granted by the license, users are provided only 
;;  with a limited warranty and the software's author, the holder of the  
;;  economic rights, and the successive licensors have only limited liability. 
;;
;;  In this respect, the user's attention is drawn to the risks associated  
;;  with loading, using, modifying and/or developing or reproducing the  
;;  software by the user in light of its specific status of free software,
;;  that may mean that it is complicated to manipulate, and that also  
;;  therefore means that it is reserved for developers and experienced 
;;  professionals having in-depth computer knowledge. Users are therefore 
;;  encouraged to load and test the software's suitability as regards their 
;;  requirements in conditions enabling the security of their systems and/or 
;;  data to be ensured and, more generally, to use and operate it in the same 
;;  conditions as regards security. 
;;
;;  The fact that you are presently reading this means that you have had 
;;  knowledge of the CeCILL-C license and that you accept its terms.
;;

;;  pel.el --- PELICANS Hierarchical Data Structure mode 
;;             for GNU Emacs and GNU XEmacs.


;;  usage: put this file in (or define a symlink to) your personal directory 
;;  ~/lisp (for example), and be sure that you have the following in 
;;  your .emacs-file :
;;
;;     ;; -1- add path where this fiel is stored
;;     ;; -2- autoload the function pel-mode from pel.el
;;     ;; -3- associate pel suffix with pel-mode                               
;;     (add-to-list load-path  (expand-file-name "~/lisp"))
;;     (autoload 'pel-mode "pel")
;;     (add-to-list 'auto-mode-alist '("\\.pel\\'" . pel-mode))
;;
 
(setq pel-mode-hook
      '(lambda ()
	 (setq pel-indentation 3
	       pel-auto-keyword-case 'upcase-word
	       pel-startup-message t
               indent-tabs-mode nil
	       )
	 )
      )

;; User options
(defcustom pel-indentation 3
  "*Extra indentation applied to MODULE blocks."
  :type 'integer
  :group 'pel
)

(defcustom pel-comment-region    "//"
  "*String inserted by \\[pel-comment-region]\ at start of each line in region."
  :type 'string
  :group 'pel
)

(defconst pel-indented-comment  "//"
  "*Comments to be indented like code.")

(defcustom pel-directive-comment "///"
  "*String of comment-like directive not to be indented."
  :type 'string
  :group 'pel
)

(defcustom pel-smart-end 'fill
  "*From an END statement, check and fill the end using matching block start.
Allowed values are 'fill and nil, which determine whether to fill the matching beginning."
  :type '(choice (const :tag "Fill" 'fill) (const :tag "none" nil))
  :group 'pel
)

(defcustom pel-startup-message t
  "*Non-nil displays a startup message when PEL mode is first called."
  :type 'boolean
  :group 'pel
)

(defvar pel-intrinsic-constants  
  '(
    "true"
    "false"
    "concrete_name"
    "type"
    ) "List of PELICANS data file intrinsic constants"
)

(defcustom pel-split-length 16
  "*Maximum number of elements in a pel operator (sub)menu."
  :type 'integer
  :group 'pel
)

(defvar pel-operators
  '(
     (" !" nil "operators")
     (" && " nil "operators")
     (" || " nil "operators")
     ("(?:)" nil "operators")
     (" << " nil "operators")
     (" <" nil "operators")
     (" <= " nil "operators")
     (" = " nil "operators")
     (" >= " nil "operators")
     (" > " nil "operators")
     (" * " nil "operators")
     (" + " nil "operators")
     (" - " nil "operators")
     (" / " nil "operators")
     ("abs(DS|IS)" nil "a...b")
     ("acos(DS)" "" "a...b")
     ("acosh(DS)" "" "a...b")
     ("apply(<vector>,<expression>,<variable>[,<component>]" "" "a...b")
     ("array(<list of IV|DV>)" "" "a...b")
     ("asin(DS)" "" "a...b")
     ("asinh(DS)" "" "a...b")
     ("atan(DS)" "" "a...b")
     ("atan2(DS,DS)" "" "a...b")
     ("atanh(DS)" "" "a...b")
     ("basename(SS)" "" "a...b")
     ("binary(SS,SS,IS)" "" "a...b")
     ("ceil(DS)" "" "c...d")
     ("component(DV|IV|BV|SV,IS)" "" "c...d")
     ("conditional_vector(<list of pairs: BS,value with all values of the same type>)" "" "c...d")
     ("cos(DS)" "" "c...d")
     ("cosh(DS)" "" "c...d")
     ("d(<expression>,SS)" "" "c...d")
     ("data_with_context(<expression>[,SS,<value>])" "" "c...d")
     ("default_roundoff(DV,IS,DS)" "" "c...d")
     ("dirname(SS)" "" "c...d")
     ("dnum(<expression>,SS)" "" "c...d")
     ("double(IS)" "" "c...d")
     ("double_equality(DS,DS,DS,DS)" "" "c...d")
     ("e()" "" "e...h")
     ("empty(SS)" "" "e...h")
     ("erf(DS)" "" "e...h")
     ("erfc(DS)" "" "e...h")
     ("euler()" "" "e...h")
     ("exp(DS)" "" "e...h")     
     ("extracted_data(SS[,SS])" "" "e...h")
     ("extracted_module(SS[,SS])" "" "e...h")
     ("floor(DS)" "" "e...h")
     ("gamma(DS)" "" "e...h")
     ("geometric_sequence(DS,DS,IS)" "" "e...h")
     ("getcwd()" "" "e...h")
     ("getenv(SS)" "" "e...h")
     ("getpid()" "" "e...h")
     ("greater(DV|IV,DS|IS)" "" "e...h")
     ("has_data(SS)" "" "e...h")
     ("has_module(SS)" "" "e...h")
     ("host_name()" "" "e...h")
     ("in_box(DV,DV,DV)" "" "i...l")
     ("in_range(DS|IS,DV|IV)" "" "i...l")
     ("incomplete_gamma(DS,DS[,DS,DS,IS])" "" "i...l")
     ("increasing(DV|IV)" "" "i...l")
     ("int(DS)" "" "i...l")
     ("interpol(DV,DV,DV) or interpol(IS,DV) " "" "i...l")
     ("is_defined(SS)" "" "i...l")
     ("j0(DS)" "" "i...l")
     ("j1(DS)" "" "i...l")
     ("jn(IS,DS)" "" "i...l")
     ("join(<list of SS>)" "" "i...l")
     ("lgamma(DS)" "" "i...l")
     ("log(DS)" "" "i...l")
     ("log10(DS)" "" "i...l")
     ("matching_color(SS,SS)" "" "m...r")
     ("max(DS,DS) " "" "m...r")
     ("middle_point(DS,DV) " "" "m...r")
     ("middle_points([DV,]DV) " "" "m...r")
     ("min(DS,DS) " "" "m...r")
     ("modulo(IS,IS)" "" "m...r")
     ("nb_ranks()" "" "m...r")
     ("nvector(IS,<scalar_type>)" "" "m...r")
     ("path_name_separator()" "" "m...r")
     ("pi()" "" "m...r")
     ("pow(DS,DS)" "" "m...r")
     ("rand()" "" "m...r")
     ("random_double()" "" "m...r")
     ("rank()" "" "m...r")
     ("regular_vector(DS,IS,DS) or regular_vector(IS,IS,IS)" "" "m...r")
     ("reverse(DV|IV|SV|BV)" "" "m...r")
     ("segm2D_sort(<DV with 2 elems>,DV,IS,DV,IS)" "" "s...t")
     ("segm3D_sort(<DV with 3 elems>,DV,IS,DV,IS,DV,IS)" "" "s...t")
     ("segm_sort(DS,DV,IS) " "" "s...t")
     ("sin(DS)" "" "s...t")
     ("sinh(DS)" "" "s...t")
     ("size(DV|IV|BV|SV)" "" "s...t")
     ("sort(DV,\">\"|\"<\")" "" "s...t")
     ("sqr(DS)" "" "s...t")
     ("sqrt(DS)" "" "s...t")
     ("stretched_vector(DS,DS,DS,DS)" "" "s...t")
     ("sum(DV|IV)" "" "s...t")
     ("tan(DS)" "" "s...t")
     ("tanh(DS)" "" "s...t")
     ("this_file_dir()" "" "s...t")
     ("to_string(DS|IS)" "" "s...t")
     ("uname()" "" "u...z")
     ("unit_sort(DS,DS,DS,IS)" "" "u...z")
     ("valid_color(SS)" "" "u...z")
     ("value(SS,<default value>)" "" "u...z")
     ("vector(<list of values with the same type>)" "" "u...z")
     ("x_cut_points(DV,DS[,DS])" "" "u...z")
     ("y_cut_points(DV,DS[,DS])" "" "u...z")
     ("y0(DS)" "" "u...z")
     ("y1(DS)" "" "u...z")
     ("yn(IS,DS)" "" "u...z")
     ("z_cut_points(DV,DS,DS)" "" "u...z")
   )
  "PELICANS Operators definition : menu + help + submenu"
)

(defvar pel-special
  '(
   "AdvectionDiffusion1_value($DV_X,$DS_a,$DS_k,$DS_uin,$DS_out)"
   "AdvectionDiffusion2_value($DV_X,$DS_T,$DV_a,$DV_X0,$DV_K)"
   "Avula_pressure_3D($DV_X,$DS_T,$DS_R,$DS_L,$DS_P1,$DS_P0,$DS_PERIOD,$DS_CUT,$DS_RHO,$DS_MU,$DS_TOL)"
   "Avula_pressure_axi($DV_X,$DS_T,$DS_R,$DS_L,$DS_P1,$DS_P0,$DS_PERIOD,$DS_CUT,$DS_RHO,$DS_MU,$DS_TOL)"
   "Avula_velocity_3D($DV_X,$DS_T,$DS_R,$DS_L,$DS_P1,$DS_P0,$DS_PERIOD,$DS_CUT,$DS_RHO,$DS_MU,$DS_TOL)"
   "Avula_velocity_axi($DV_X,$DS_T,$DS_R,$DS_L,$DS_P1,$DS_P0,$DS_PERIOD,$DS_CUT,$DS_RHO,$DS_MU,$DS_TOL)"
   "Beltrami1_grad_velocity($DV_X,$DS_T,$DV_K,$DS_A,$DS_NU)"
   "Beltrami1_pressure($DV_X,$DS_T,$DV_K,$DS_A,$DS_NU)"
   "Beltrami1_velocity($DV_X,$DS_T,$DV_K,$DS_A,$DS_NU)"
   "GreenTaylor_grad_velocity($DV_X,$DS_T,$DS_NU)"
   "GreenTaylor_pressure($DV_X,$DS_T,$DS_NU)"
   "GreenTaylor_velocity($DV_X,$DS_T,$DS_NU)"
   "NavierStokes1_force($DV_X,$DS_ALPHA,$DS_MU)"
   "NavierStokes1_grad_velocity($DV_X,$DS_ALPHA,$DS_MU)"
   "NavierStokes1_pressure($DV_X,$DS_ALPHA,$DS_MU)"
   "NavierStokes1_velocity($DV_X,$DS_ALPHA,$DS_MU)"
   "TransientDiffusion1($DS_T)"
   "SolidBodyRotation_value($DV_X,$DS_T)"
   "SolidBodyRotation_velocity($DV_X,$DS_T)"
   "Stokes2_force($DV_X)"
   "Stokes2_grad_velocity($DV_X)"
   "Stokes2_pressure($DV_X)"
   "Stokes2_velocity($DV_X)"
   "VariableDensityFlow1_grad_rho($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"
   "VariableDensityFlow1_grad_rho_velocity($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"
   "VariableDensityFlow1_grad_velocity($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"
   "VariableDensityFlow1_pressure($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"
   "VariableDensityFlow1_rho($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"
   "VariableDensityFlow1_rho_velocity($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"
   "VariableDensityFlow1_rhsg($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"
   "VariableDensityFlow1_rhsf($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"
   "VariableDensityFlow1_velocity($DV_X,$DS_T,$DS_ALPHA,$DS_BETA,$DS_MU,$DS_X0,$DS_X1,$DS_C,$DS_L)"
   )
"PELICANS special insert (tests cases)"
)

(defvar pel-keywords 
  '(
    "module"
    "end"
    ) "*List of PELICANS data file keywords."
)

(defconst pel-sep "[ \t\n]*")
(defconst pel-symbol-re "\\([a-z_][a-z_0-9#]*\\)")
(defconst pel-vars-re "\\<\\([$][BDIS][SV]_[a-z_][a-z_0-9-]*\\)\\>")
(defconst pel-comment-re "\\(//\\)")
;; begin and end of major blocks
(defconst pel-begin-re "\\(module\\)\\>")
(defconst pel-end-re "\\(end module\\)\\>")
;; any block match (all the previous ones)
(defconst pel-block-re "\\<\\(module\\|end module\\)\\>")

;; Highlighting patterns
(defvar pel-intrinsic-procedures nil "*List of PELICANS data file intrinsic operators")

(defun pel-font-lock-keywords-1 ()
  (list
   ;; special declarations
   '("\\<module[ \t]*\\(\\sw+\\)" 1 font-lock-type-face)
   )
  )

(defun pel-font-lock-keywords-2 ()
  (list
   '("\\<module[ \t]*\\(\\sw+\\)" 1 font-lock-type-face)
   (list
    (concat "\\<\\(" (mapconcat 'identity pel-intrinsic-procedures "\\|") "\\)\\>") 1 font-lock-function-name-face
    )
   (list
    (concat "\\<\\(" (mapconcat 'identity pel-intrinsic-constants "\\|") "\\)\\>") 1 font-lock-reference-face
    )
   (list pel-vars-re 1 'font-lock-variable-name-face)
   (concat "\\<\\(" (mapconcat 'identity pel-keywords "\\|") "\\)\\>")
   )
  )

(defvar pel-font-lock-keywords   (purecopy (pel-font-lock-keywords-2))
  "*Additional expressions to highlight in PEL mode.")

;; syntax table
(defvar pel-mode-syntax-table nil
  "Syntax table in use in PEL mode buffers.")

(if pel-mode-syntax-table
    ()
  (setq pel-mode-syntax-table (make-syntax-table))
  (modify-syntax-entry ?/ ". 12" pel-mode-syntax-table)  ; beg C++ comment style
  (modify-syntax-entry ?\n "> a" pel-mode-syntax-table)  ; end comment C++ comment style
  (modify-syntax-entry ?_ "w" pel-mode-syntax-table)   ; underscore in names
  (modify-syntax-entry ?- "." pel-mode-syntax-table)   ; minus sign in names
  (modify-syntax-entry ?. "." pel-mode-syntax-table)   ; minus sign in names
  (modify-syntax-entry ?# "w" pel-mode-syntax-table)   ; diesis sign in names
  (modify-syntax-entry ?\' "\"" pel-mode-syntax-table) ; string quote
  (modify-syntax-entry ?\" "\"" pel-mode-syntax-table) ; string quote
  (modify-syntax-entry ?\r " " pel-mode-syntax-table)  ; return is whitespace
  (modify-syntax-entry ?+ "." pel-mode-syntax-table)   ; punctuation
  (modify-syntax-entry ?= "." pel-mode-syntax-table)
  (modify-syntax-entry ?* "." pel-mode-syntax-table)
;  (modify-syntax-entry ?/ "." pel-mode-syntax-table)
  (modify-syntax-entry ?\\ "/" pel-mode-syntax-table)
)

;; keys
(defvar pel-mode-map ()
  "Keymap used in PEL mode.")
(if pel-mode-map
    ()
  (setq pel-mode-map (make-sparse-keymap))
  (define-key pel-mode-map "\C-c\C-c"    'pel-comment-region)
  (define-key pel-mode-map "\C-c\C-u"    'pel-uncomment-region)
  (define-key pel-mode-map "\C-\M-a"  'pel-beginning-of-module)
  (define-key pel-mode-map "\C-\M-e"  'pel-end-of-module)
  (define-key pel-mode-map "\C-\M-h"  'pel-mark-module)
  (define-key pel-mode-map "\C-\M-q"  'pel-indent-module)
  (define-key pel-mode-map "\C-j"     'pel-indent-new-line) ; LFD equals C-j
  (define-key pel-mode-map "\r"       'newline-and-indent)
  (define-key pel-mode-map "\C-c\C-p" 'pel-previous-statement)
  (define-key pel-mode-map "\C-c\C-n" 'pel-next-statement)
  (define-key pel-mode-map "\t"       'pel-indent-line)
  (define-key pel-mode-map "\C-c#" 'goto-char)
)
;; menus
(defvar pel-main-menu
  '("PEL" 
    ["Beautify module"     pel-beautify-module t]
    ["Indent module"       pel-indent-module t]
    ["Mark module"         pel-mark-module t]
    ["Beginning of module" pel-beginning-of-module t]
    ["End of module"       pel-end-of-module t]
    "-----"
    ["Comment region"      pel-comment-region t]
    ["Uncomment region"    pel-uncomment-region t]
    ["Indent region"       indent-region t]
    "-----"
    ["Customize"       (customize-group "pel") t]
    )
  "The PEL main menu")
;; menus
(defvar pel-insert-menu
  '("PEL-Insert" 
    ["Help"     pel-switch-to-doc-buffer t]
    "-----"
    )
  "The PEL insert menu")


(defun pel-insert (string) 
  (let ((lstring (replace-in-string string "\\(\\sw\\|[<|> ]\\|[\[]\\|[\]]\\)+\\([,\)]\\)" "\\2")))
    (insert lstring)
    (if (string= lstring string)
	nil
      (search-backward "(")
      (forward-char)
      )
    )
  )

(defvar pel-short-docs "Single line help for PEL operators" "")

(defun pel-install-menubar ()
  "Installs the PEL menu at the menubar."
  (let ((short-docs ""))
    (if (and (featurep 'menubar) current-menubar)
	(progn
	  (let (
		(buffer-menubar (copy-sequence current-menubar))
		)
	    (set-buffer-menubar buffer-menubar)
	    )
	  (add-submenu nil pel-main-menu)
	  (add-submenu nil pel-insert-menu)
	  ; ajout des operateurs pelicans
	  (mapcar
	   (function
	    (lambda (entry)
	      (let (
		    (name (car entry))
		    (desc (nth 1 entry))
		    (group (nth 2 entry))
		    sname
		    )
		(if (not desc) (setq desc ""))
		(if (not group) (setq group "*"))
		(if (string-match (concat "^\\(" pel-symbol-re "\\)") name)
		    (progn
		      (setq sname (substring name 0 (match-end 1)))
		      (setq pel-intrinsic-procedures (cons sname pel-intrinsic-procedures))
		      (add-menu-button
		       (list (car pel-insert-menu) group)
		       (vector sname (list 'pel-insert name) t))
		      )
		  (add-menu-button
		   (list (car pel-insert-menu) group)
		   (vector name (list 'insert name) t))
		  )
		(if desc
		    (setq short-docs (concat short-docs "\n" name ": " desc))
		  )
		)
	      )
	    )
	   pel-operators)

	  (add-menu-button (list (car pel-insert-menu)) "--")

	  ; ajout des cas tests pelicans
	  (let (
		name
		group
		sname
		(remain pel-special)
		)
	    (while remain
	      (setq name (pop remain))
	      (string-match "^\\(^[A-Za-z0-9]+\\)\\(_\\|\(\\)" name)
	      (setq group (substring name 0 (match-end 1)))
	      (string-match (concat "^\\(" pel-symbol-re "\\)") name)
	      (setq sname (substring name 0 (match-end 1)))
	      (setq pel-intrinsic-procedures (cons sname pel-intrinsic-procedures))
	      (add-menu-button
	       (list (car pel-insert-menu) group)
	       (vector sname (list 'insert name) t))
	      )
	    )
	  )
      )
    (setq pel-font-lock-keywords (purecopy (pel-font-lock-keywords-2)))
    (put 'pel-short-docs 'variable-documentation short-docs)
    )
  )

(add-hook 'post-command-hook 'pel-get-help)
(add-hook 'pre-command-hook 'pel-refresh-echo-area)
(add-hook 'pel-mode-hook 'pel-install-menubar)


;; A temporary position to make region operators faster
(defvar pel-cache-position nil)
(make-variable-buffer-local 'pel-cache-position)

;;;###autoload
(defun pel-mode ()
"Major mode for editing PELICANS data files.

\\[pel-indent-new-line] corrects current indentation and creates new indented line.
\\[pel-indent-line] indents the current line correctly. 
\\[pel-indent-module] indents the current module. 

Key definitions:
\\{pel-mode-map}

Variables controlling indentation style and extra features:

 pel-indentation
    Extra indentation within MODULES blocks.  (default 3)
 pel-comment-region
    String inserted by \\[pel-comment-region] at start of each line in 
    region.  (default \"//\")
 pel-indented-comment
    String holding the type of comment to be intended like code.
    (default \"//\")
 pel-directive-comment
    String of comment-like directive not to be indented (column 0).
    (default \"///\")
 pel-smart-end 
    From an END statement, check and fill the end using matching block start.
    Allowed values are 'fill and nil, which determine
    whether to fill the matching beginning.) (default 'fill)
 pel-auto-keyword-case
    Automatic change of case of keywords. (default nil)
    The possibilities are 'downcase-word, 'upcase-word, 'capitalize-word.
 pel-startup-message
    Set to nil to inhibit message first time PEL mode is used. (default t)
 pel-keywords
    List of keywords used for highlighting/upcase-keywords etc.

Turning on PEL mode calls the value of the variable `pel-mode-hook'
with no args, if that value is non-nil."
  (interactive)
  (kill-all-local-variables)
  (setq major-mode 'pel-mode)
  (setq mode-name "PEL")
  (set-syntax-table pel-mode-syntax-table)
  (use-local-map pel-mode-map)
  (make-local-variable 'indent-line-function)
  (setq indent-line-function 'pel-indent-line)
  (make-local-variable 'indent-region-function)
  (setq indent-region-function 'pel-indent-region)
  (make-local-variable 'require-final-newline)
  (setq require-final-newline t)
  (make-local-variable 'comment-start)
  (setq comment-start "//")
  (make-local-variable 'comment-start-skip)
  (setq comment-start-skip "// *")
  (make-local-variable 'comment-indent-function)
  (setq comment-indent-function 'pel-comment-indent)
  (setq indent-tabs-mode nil)
  ;; Setting up things for font-lock
  (if (string-match "Lucid" emacs-version)
      (put 'pel-mode 'font-lock-keywords-case-fold-search t)
    (make-local-variable 'font-lock-keywords) ; for Emacs version <= 19.28 
    (setq font-lock-keywords pel-font-lock-keywords)
    ;; (make-local-variable 'font-lock-defaults) ; for Emacs version > 19.28
    ;; (setq font-lock-defaults '(pel-font-lock-keywords t))
    )
  (make-local-variable 'font-lock-keywords-case-fold-search)
  (setq font-lock-keywords-case-fold-search t)
  (run-hooks 'pel-mode-hook)
  (if pel-startup-message
      (message "Emacs mode for PELICANS Hierarchical Data Structures.")
    )
  (setq pel-startup-message nil)
)

;; inline-functions
(defsubst pel-get-beg-of-line ()
  (save-excursion (beginning-of-line) (point))
)

(defsubst pel-get-end-of-line ()
  (save-excursion (end-of-line) (point))
)

(defsubst pel-in-string ()
  (let (
	(beg-pnt
	 (if (and pel-cache-position (> (point) pel-cache-position))
	     pel-cache-position
	   (point-min)
	   )
	 )
	)
    (nth 3 (parse-partial-sexp beg-pnt (point)))
  )
)
	    
(defsubst pel-in-comment ()
  (let (
	(beg-pnt
	 (if (and pel-cache-position (> (point) pel-cache-position))
	     pel-cache-position
	   (point-min))
	)
       )
    (nth 4 (parse-partial-sexp beg-pnt (point)))
  )
)


(defsubst pel-match-piece (arg)
  (if (match-beginning arg)
      (buffer-substring (match-beginning arg) (match-end arg))))


(defsubst pel-looking-at-module ()
  "Return (\"module\" name) if a module starts after point."
  (let (
	struct
	(label nil)
	)
    (if (looking-at (concat pel-begin-re pel-sep pel-symbol-re))
	(list (pel-match-piece 1) (pel-match-piece 2))
      )
    )
  )

(defsubst pel-looking-at-module-end ()
  "Return (\"module\" name) if a at module end."
  (let (
	struct
	(label nil)
	)
  (if (looking-at (concat "end[ \t]*" pel-begin-re "?\\([ \t]+\\("
			  pel-symbol-re  "\\)\\)?\\>"))
      (list (pel-match-piece 1) (pel-match-piece 3))
      )
    )
  )


(defsubst pel-line-continued ()
  (save-excursion
    (let (
	  (bol (pel-get-beg-of-line))
	  )
      (end-of-line)
      (while (pel-in-comment)
	(re-search-backward pel-comment-re bol)
	(skip-chars-backward "/")
	)
      (skip-chars-backward " \t")
      (= (preceding-char) ?&)
      )
    )
  )

(defsubst pel-present-statement-cont ()
  "Return continuation properties of present statement."
  (let (
	pcont
	cont
	)
    (save-excursion
      (setq pcont (if (pel-previous-statement) (pel-line-continued) nil))
      )
    (setq cont (pel-line-continued))
    (cond ((and (not pcont) (not cont)) 'single)
 	  ((and (not pcont) cont)       'begin)
 	  ((and pcont       (not cont)) 'end)
 	  ((and pcont       cont)       'middle)
 	  (t (error))
	  )
    )
  )

(defsubst pel-update-line ()
  (let (
	bol
	eol
	)
    (if pel-auto-keyword-case
	(progn
	  (setq bol (pel-get-beg-of-line) eol (pel-get-end-of-line))
	  (pel-change-keywords pel-auto-keyword-case bol eol)
	  )
      )
    )
  )
	  

;; Statement = statement line, a line which is neither blank, nor a comment.
(defun pel-previous-statement ()
  "Move point to beginning of the previous PEL statement.
Return nil if no previous statement is found."
  (interactive)
  (let (
	not-first-statement
	)
    (beginning-of-line)
    (while (and
	    (setq not-first-statement (zerop (forward-line -1)))
	    (looking-at "\\([ \t]*$\\|[ \t]*//\\)")
	    )
      )
    not-first-statement
    )
  )

(defun pel-next-statement ()
  "Move point to beginning of the next PEL statement.
Return nil if no later statement is found."
  (interactive)
  (let (
	not-last-statement
	)
    (beginning-of-line)
    (while (and
	    (setq not-last-statement (and (zerop (forward-line 1)) (not (eobp))))
	    (looking-at "\\([ \t]*\\($|//\\)\\)")
	    )
      )
    not-last-statement
    )
  )

(defun pel-beginning-of-module ()
  "Move point to the beginning of module. Return (type name) or nil if not found."
  (interactive)
  (let (
	(count 1)
	(case-fold-search t)
	pos
	)
    (while (and (not (zerop count)) (re-search-backward pel-block-re nil 'move))
      (beginning-of-line) (skip-chars-forward " \t")
      (setq pos (point))
      (cond
       ((looking-at (concat pel-begin-re pel-sep pel-symbol-re))
	(setq count (- count 1))
	)
       ((looking-at (concat pel-end-re pel-sep pel-symbol-re))
	(setq count (+ count 1))
	)
       )
      )
    (if (zerop count)
	(progn (goto-char pos) (pel-looking-at-module))
      (message "No beginning-found.")
      nil
      )
    )
  )


(defun pel-end-of-module ()
  "Move point to the end of module. Return (type name) or nil if not found."
  (interactive)
  (let (
	(count 1)
	(case-fold-search t)
	pos
	)
    (while (and 
	    (not (zerop count)) 
	    (re-search-forward pel-block-re nil 'move) )
      (beginning-of-line) (skip-chars-forward " \t")
      (cond
       ((looking-at (concat pel-begin-re pel-sep pel-symbol-re))
	(setq count (+ count 1))
	)
       ((looking-at (concat pel-end-re pel-sep pel-symbol-re))
	(setq count (- count 1))
	)
       )
      (end-of-line)
      (setq pos (point))
      )
    (if (zerop count)
	(progn (goto-char pos) pos)
      (message "No end found." )
      nil
      )
    )
  )

(defun pel-mark-module ()
  "Put mark at end of PEL module, point at beginning.
Marks are pushed and highlight (grey shadow) is turned on."
  (interactive)
  (let (
	(pos (point))
	program
	)
    (pel-end-of-module)
    (push-mark (point) t)
    (goto-char pos)
    (setq program (pel-beginning-of-module))
;   (push-mark (point))
    ;; The keywords in the preceding lists assume case-insensitivity.
    (if (string-match "Lucid" emacs-version)
	(zmacs-activate-region)
      (setq mark-active t)
      (setq deactivate-mark nil)
      )
    program
    )
  )

(defun pel-comment-line ()
  "comment current line"
  (beginning-of-line) (skip-chars-forward " \t")
  (cond
   ((looking-at (regexp-quote pel-directive-comment))
    )
   (t
    (beginning-of-line)
    (insert pel-comment-region)
    )
   )
  )

(defun pel-uncomment-line ()
  "uncomment current line"
  (beginning-of-line) (skip-chars-forward " \t")
  (cond
   ((looking-at (regexp-quote pel-directive-comment))
    )
   ((looking-at (regexp-quote pel-comment-region))
    (delete-region (point) (match-end 0))
    )
   )
  )

(defun pel-comment-region (beg-region end-region)
  "Comment/uncomment every line in the region.
Insert pel-comment-region at the beginning of every line in the region."
  (interactive "*r")
  (let (
	(end (make-marker))
	)
    (set-marker end end-region)
    (goto-char beg-region)
    (pel-comment-line)
    (while (and  (zerop (forward-line 1)) (< (point) (marker-position end)))
      (pel-comment-line)
      )
    (set-marker end nil)
    )
  )

(defun pel-uncomment-region (beg-region end-region)
  "Uuncomment every line in the region.
Remove pel-comment-region at the beginning of every line in the region
or."
  (interactive "*r")
  (let (
	(end (make-marker))
	)
    (set-marker end end-region)
    (goto-char beg-region)
    (pel-uncomment-line)
    (while (and  (zerop (forward-line 1)) (< (point) (marker-position end)))
      (pel-uncomment-line)
      )
    (set-marker end nil)
    )
  )

(defsubst pel-current-indentation ()
  "Return indentation of current line."
  (save-excursion
    (beginning-of-line)
    (skip-chars-forward " \t")
    (current-column)
    )
  )


(defsubst pel-indent-to (col &optional)
  "Indent current line to column COL.
If no-line-number nil, jump over a possible line-number."
  (beginning-of-line)
  (delete-horizontal-space)
  (indent-to col col)
)

(defun pel-calculate-indent ()
  "Calculate the indent column based on previous statements."
  (interactive)
  (let (
	icol
	cont
	(case-fold-search t)
	pnt
	motif
	beg
	)
    (save-excursion
      (beginning-of-line)
      (setq pnt (point))
      (if (not (pel-previous-statement))
	  (setq icol 0)
	(setq cont (pel-present-statement-cont))
	(if (eq cont 'end)
	    (while (not (eq 'begin (pel-present-statement-cont)))
	      (pel-previous-statement)
	      )
	  )
	(cond ((eq cont 'begin)
	       (setq icol (+ (pel-current-indentation) pel-continuation-indent)))
	      ((eq cont 'middle)
	       (setq icol (current-indentation)))
	      (t
	       (setq icol (pel-current-indentation))
	       (setq beg (point))
	       (skip-chars-forward " \t")
	       (if (looking-at "end module\\>")
		   (setq icol (+ icol pel-indentation))
		 )
	       (while (and (<= beg pnt) (re-search-forward pel-block-re nil 'move))
		 (setq beg (match-end 0))
		 (if (<= beg pnt)
		     (progn
		       (setq motif (match-string 0))
;		       (message "motif %s" motif)
		       (cond
			((string-match pel-end-re motif) (setq icol (- icol pel-indentation)))
			((string-match pel-begin-re motif) (setq icol (+ icol pel-indentation)))
;			((string-match pel-sbegin-re motif) (setq icol (+ icol pel-indentation)))
;			((string-match pel-send-re motif) (setq icol (- icol pel-indentation)))
			)
		   )
		   )
		 )
		   (goto-char pnt)
		   (beginning-of-line)
		   (skip-chars-forward " \t")
		   (if (looking-at "end\\>")
		       (setq icol (- icol pel-indentation))
		     )

	       )
	      )
	)
      )
    icol
    )
  )

(defsubst pel-comment-indent ()
  (cond ;;((looking-at "\\(!!!\\|^[$]\\)") 0)
	((looking-at pel-directive-comment) 0)
	((looking-at pel-indented-comment) (pel-calculate-indent))
	(t
	 (skip-chars-backward " \t")
	 (max (if (bolp) 0 (1+ (current-column))) comment-column)
	 )
	)
  )
(defsubst pel-equal-symbols (a b)
  "Compare strings neglecting case and allowing for nil value."
  (let ((a-local (if a (downcase a) nil))
	(b-local (if b (downcase b) nil)))
    (equal a-local b-local)))

(defun pel-block-match (beg-block beg-name end-block end-name)
  "Match end-module with beg-module and complete end-block if possible.
Leave point at the end of line."
  (search-forward "end" (pel-get-end-of-line))
  (catch 'no-match
    (if (not (pel-equal-symbols beg-block end-block))
	(if end-block
	    (progn
	      (message "END %s does not match %s." end-block beg-block)
	      (end-of-line) 
	      (throw 'no-match nil))
	  (message "Inserting %s." beg-block)
	  (insert (concat " " beg-block)))
      (search-forward end-block))
    (if (not (pel-equal-symbols beg-name end-name))
	(cond ((and beg-name (not end-name)) 
	       (message "Inserting %s." beg-name)
	       (insert (concat " " beg-name)))
	      ((and beg-name end-name) 
	       (message "Replacing %s with %s." end-name beg-name)
	       (search-forward end-name)
	       (replace-match beg-name t))
	      ((and (not beg-name) end-name) 
	       (message "Deleting %s." end-name)
	       (search-forward end-name)
	       (replace-match "")))
      (if end-name (search-forward end-name)))
    (if (not (looking-at "[ \t]+//")) (delete-region (point) (point-at-eol)))
    )
  )


(defun pel-match-end ()
  "From an end foo statement, find the corresponding foo including name."
  (interactive)
  (let ((count 1)
	(matching-beg nil)
	(end-point (point))
	(case-fold-search t)
	beg-name end-name beg-block end-block end-struct)

    (if (save-excursion
	  (beginning-of-line)
	  (skip-chars-forward " \t")
	  (setq end-struct (pel-looking-at-module-end))
	  )
	(progn
	  (setq end-block (car end-struct))
	  (setq end-name  (car (cdr end-struct)))
	  ; (message "end_block=%s end-name=%s" end-block end-name)
	  (save-excursion
	    (beginning-of-line)
	    (while (and (not (zerop count))
			(re-search-backward 
			 (concat "\\(" pel-begin-re "\\)") nil t))
	      (beginning-of-line) (skip-chars-forward " ")
	      (cond ((setq matching-beg
			   (cond
			    ((pel-looking-at-module))))
		     (setq count (- count 1)))
		    ((looking-at (concat "end[ \t]*" pel-begin-re "\\b"))
		     (setq count (+ count 1)))
		    )
	      )
	    (if (not (zerop count))
		(message "No starting block match found.")
	      (pel-update-line)

	      (setq beg-block (car matching-beg))
	      (setq beg-name (car (cdr matching-beg)))
	      (goto-char end-point)
	      (beginning-of-line)
	      (pel-block-match beg-block beg-name end-block end-name)
	      )
	    )
	  )
      ; (message "end_struct" end-struct)
      )
    )
  )


(defun pel-indent-line (&optional no-update)
  "Indent current line as PEL code."
  (interactive)
  (let (
	(indent 0)
	(no-line-number nil)
	(pos (make-marker))
	(case-fold-search t)
	)
    (set-marker pos (point))
    (beginning-of-line)
    (skip-chars-forward " \t")
    (if (looking-at pel-comment-re)
	(setq indent (pel-comment-indent))
      (if (and (looking-at "end") pel-smart-end) (pel-match-end))
      (setq indent (pel-calculate-indent))
    )
    ;; (message "indent=%d col=%d" indent (current-column))
      (pel-indent-to indent)
    ;; If initial point was within line's indentation,
    ;; position after the indentation.  Else stay at same point in text.
    (if (< (point) (marker-position pos))
	(goto-char (marker-position pos))
    )
    (if (not no-update) (pel-update-line))
    (set-marker pos nil)
  )
)


(defun pel-indent-new-line ()
  "Reindent the current PEL line, insert a newline and indent the newline.
If run in the middle of a line, the line is not broken."
  (interactive)
  (let (
	string
	cont
	(case-fold-search t)
	)
    (beginning-of-line)			; Reindent where likely to be needed.
    (if (or (looking-at "\\(end module[ \t]*\\|!\\)"))
	(pel-indent-line 'no-update))
    (end-of-line)
    (delete-horizontal-space)		;Destroy trailing whitespace
    (setq string (pel-in-string))
    (pel-update-line)
    (newline)
    (pel-indent-line 'no-update)
    )
  )


(defun pel-indent-region (beg-region end-region)
  "Indent every line in region by forward parsing."
  (interactive "*r")
  (let (
	(end-region-mark (make-marker))
	(save-point (point-marker))
	ind-lev
	ind-curr
	ind-b
	)
    (set-marker end-region-mark end-region)
    (goto-char beg-region)
    ;; first find a line which is not a continuation line or comment
    (beginning-of-line)
    (while (and ;(looking-at "\\([ \t]*$\\)")
		(progn (pel-indent-line );'no-update)
		       (zerop (forward-line 1))
		)
		(< (point) end-region-mark)
           )
    )
    (if (string-match "Lucid" emacs-version)
	(zmacs-deactivate-region)
      (deactivate-mark)
    )
  ) ;let
)

(defun pel-indent-module ()
  "Properly indent the module which contains point."
  (interactive)
  (save-excursion
    (let (program)
      (setq program (pel-mark-module))
      (if program
	  (progn
	    (message (concat "Indenting " (car program) " " (car (cdr program))"."))
	    (pel-indent-region (point) (mark))
	    (message (concat "Indenting " (car program) " " (car (cdr program)) "...done."))
          )
	(message "Indenting the whole file.")
	(pel-indent-region (point) (mark))
	(message (concat "Indenting the whole file...done."))
      )
    )
  )
)

(defun pel-beautify-module ()
  "Properly indent & downcase the module wich contains point."
  (interactive)
  (pel-indent-module)
  (save-excursion
    (let (program beg end)
      (setq program (pel-mark-module))
      (setq beg (point))
      (setq end (mark))
      (if program
	  (progn
	    (message "Beautify %s %s..." (car program) (car (cdr program)))
	    (pel-beautify-region beg end)
	    (message "Beautify %s %s...done" (car program) (car (cdr program)))
	  )
	(message "Beautify the whole file...")
	(message "Beautify the whole file...done")
      )
    )
  )
)



(defun pel-upcase-keywords ()
  "Upcase all PEL keywords in the buffer."
  (interactive)
  (pel-change-keywords 'upcase-word)
)


(defun pel-upcase-region-keywords (beg end)
  "Upcase all PEL keywords in the region."
  (interactive "*r")
  (pel-change-keywords 'upcase-word beg end)
)


;; Change the keywords according to argument.
(defun pel-change-keywords (change-word &optional beg end)
  (save-excursion
    (setq beg (if beg beg (point-min)))
    (setq end (if end end (point-max)))
    (let (
	  (keyword-re 
	   (concat "\\<\\(" (mapconcat 'identity pel-keywords "\\|")
		   "\\)\\>")
	   )
	  (ref-point (point-min))
	  state
	  )
      (goto-char beg)
      (while (re-search-forward keyword-re end t)
	(if
	    (progn
	      (setq state (parse-partial-sexp ref-point (point)))
	      (or (nth 3 state) (nth 4 state)
		  (save-excursion	; Check for cpp directive.
		    (beginning-of-line)
		    (skip-chars-forward " \t")
		    (looking-at "#")))
	      )
	    ()
	  (setq ref-point (point))
	  (funcall change-word -1)
	  )
	)
      )
    )
  )
(defun pel-beautify-region (beg end)
  (save-excursion
    (let (
	  (tbeg beg)
	  (tend end)
	  p
	  )
      (goto-char end) (backward-char 1) (insert " _TAGEND_")
      (goto-char beg)
      (while (re-search-forward "\\([!\"']\\|^\\$\\|[ \t]\\)" end t)
	(backward-char 1)
	(cond
	 ((looking-at "[ \t]_TAGEND_") (replace-match  "") (goto-char end))
	 ((looking-at "//") (end-of-line))
	 ((looking-at "\"") (forward-char 1) (re-search-forward "[\"]" end t))
	 ((looking-at "'")  (forward-char 1) (re-search-forward "[']"  end t))
	 ((looking-at "[ \t]+")
	  (if (eq (current-column) 0)
	      (skip-chars-forward " \t")
	    (setq p (point))
	    (replace-match  " ")
	    (setq tend (+ tend (- p (point))))
	    )
	  )
	 )
	)
      )
    )
  )

;;; Split LIST into sublists of max length N
;;; Example (fume-split '(1 2 3 4 5 6 7 8) 3)-> '((1 2 3) (4 5 6) (7 8))
;;;

(defun pel-split (list n)
  (let ((i 0)
        result
        sublist
        (remain list))
    (while remain
      (if (= n (setq sublist (cons (car remain) sublist)
                     remain (cdr remain)
                     i (1+ i)))
          ;; We have finished a sublist
          (setq result (cons (nreverse sublist) result)
                sublist nil
                i 0)))
    ;; There might be a sublist (if the length of LIST mod n is != 0)
    ;; that has to be added to the result list.
    (if sublist
        (setq result (cons (nreverse sublist) result)))
    (nreverse result)))

(defvar pel-doc-buffer " *pel-doc*"
  "Where the documentation can be found.")

(defvar pel-last-help nil
  "The last help message, for echo area refresh.")

(make-variable-buffer-local 'pel-last-help)

(defun pel-get-help ()
  "Get one-line docs on the symbol at the point.
The data for these docs is a little bit obsolete and may be in fact longer
than a line. Your contribution to update/shorten it is appreciated."
  (interactive)
  (save-match-data			; May be called "inside" query-replace
    (save-excursion
      (let ((word (pel-word-at-point-hard)))
;	(message "get help on %s" word)
	(if word
	      (pel-describe-symbol word))
))))

(defun pel-word-at-point-hard ()
  (current-word t)
)


(defun pel-describe-symbol (val)
  "Display the documentation of symbol at point."
  (let ((enable-recursive-minibuffers t)
	args-file regexp)

    (setq regexp (concat "^" 
			 (regexp-quote val) 
			 "\\([ \t\(]\\|$\\)"
			 )
	  )

    ;; get the buffer with the documentation text
    (pel-switch-to-doc-buffer)

    ;; lookup in the doc
    (goto-char (point-min))
    ;(message "searching %s" regexp)
    (let ((case-fold-search nil))
      (list 
       (if (re-search-forward regexp (point-max) t)
	   (save-excursion
	     (beginning-of-line 1)
	     (let ((lnstart (point)))
	       (end-of-line)
	       (setq pel-last-help
		     (pel-message "%s" (buffer-substring lnstart (point))))))
)))))

(defun pel-switch-to-doc-buffer ()
  "Go to the perl documentation buffer and insert the documentation."
  (interactive)
  (let ((buf (get-buffer-create pel-doc-buffer)))
    (if (interactive-p)
	(switch-to-buffer-other-window buf)
      (set-buffer buf))
    (if (= (buffer-size) 0)
	(progn
	  (insert (documentation-property 'pel-short-docs
					  'variable-documentation))
	  (setq buffer-read-only t)))))

(defun pel-refresh-echo-area ()
  (and pel-last-help
       (if (and ;pel-mode
		(not executing-kbd-macro)
		(not cursor-in-echo-area)
		(not (eq (selected-window) (minibuffer-window))))
	   (pel-message pel-last-help)
	 (setq pel-last-help nil))))

(defun pel-message (&rest args)
  (let ((omessage pel-last-help))
    (cond ((eq (car args) pel-last-help))
          ((or (null args)
               (null (car args)))
           (setq pel-last-help nil))
          (t
           (setq pel-last-help (apply 'format args))))
    ;; Do not put pel-help messages in the log
    (if pel-last-help
	(display-message 'no-log pel-last-help)
      (and omessage
	   (clear-message 'no-log))))
  pel-last-help)

(provide 'pel)

;;; pel.el ends here

