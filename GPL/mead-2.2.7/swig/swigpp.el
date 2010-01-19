;;; swigpp.el --- make swig interface files from annotated C++ headers
;;; ID: $Id: swigpp.el,v 1.61 2001/06/26 18:32:24 bergsma Exp $
;;; SOURCE: $Source: /cvs-repository/bashford/cvsroot/mead/swig/swigpp.el,v $

;;; Commentary:

;; This file provides the command `swigpp-process-file' and its batch variant,
;; `batch-swigpp-process-file'.  Args are IN-FILE-NAME (usually "../Foo.h"),
;; OUT-FILE-NAME (usually "Foo.j") and MODULE-NAME (usually "Foo").  See the
;; GNUmakefile in this directory for `batch-swigpp-process-file' usage.
;;
;; The input file is searched for class declarations prefixed with "//!wrap!"
;; directives.  For each such class, the public members are processed to form
;; a swig-compatible class declaration in the output file.  Class members
;; prefixed w/ "//!nowrap!" are ignored, as are all text between "//!nowrap!+"
;; and "//!nowrap!-".
;;
;; In the output, ctors and dtors are grouped first, followed by operators,
;; followed by other class members.
;;
;; Additionally, lines between "#if SWIGPP_LITERAL_INCLUDE" and
;; "#endif // SWIGPP_LITTERAL_INCLUDE" are included directly in the output
;; without any processing.

;;; Code:

(require 'cl)				; sanity

(defvar swigpp-debug (getenv "SWIGPP_DEBUG")
  "Non-nil means swigpp.el writes debugging info to stdout.")

(defun swigpp-message-buffer (&optional buf)
  (with-current-buffer (or buf (current-buffer))
    (message "[%s] %d lines<<<\n%s\n>>>"
             (buffer-name)
             (count-lines (point-min) (point-max))
             (buffer-string))))

;;;---------------------------------------------------------------------------
;;; Textual ignorance

(defvar swigpp-wrap-one "//!wrap!"
  "Directive to wrap the next class or member.")
(defvar swigpp-nowrap-one "//!nowrap!"
  "Directive to not wrap the next class or member.")
(defvar swigpp-nowrap-start "//!nowrap!+"
  "Directive to start not wrapping successive following members.")
(defvar swigpp-nowrap-end "//!nowrap!-"
  "Directive to stop not wrapping successive members.")

(defun swigpp-next-member ()
  "Advance point past newline following end of member.  Return region.
If the ending regexp is not found, return nil."
  (let ((beg (progn (forward-comment (buffer-size))
		    (beginning-of-line)
		    (point))))
    (when (re-search-forward ";.*\n" (point-max) t)
      (when swigpp-debug
	(message "Next member: <<<%s>>>" (buffer-substring beg (point))))
      (cons beg (point)))))

(defun swigpp-subvert-definition-blocks ()
  ;; Replace "{...}" with ";", stashing text in `block' text property.
  (goto-char (point-min))
  (let (block)
    (while (search-forward "{" (point-max) t)
      (forward-char -1)
      (let ((beg (point)))
        (forward-sexp 1)
        (end-of-line)
        (setq block (buffer-substring-no-properties beg (point)))
        (delete-region beg (point)))
    (insert ";")
    (put-text-property (1- (point)) (point) 'block block))))

(defun swigpp-mask-non-public ()
  (let (end)
    (goto-char (point-max))
    (while (progn (setq end (point))
		  (re-search-backward
		   "^\\(\\(public\\)\\|\\(private\\)\\|\\(protected\\)\\):"
		   (point-min) t))
      (let ((md (match-data)))
	(when (or (nth 6 md) (nth 8 md)) ; private or protected
	  (goto-char (match-beginning 0))
	  (when swigpp-debug
	    (message "Masking non-public: <<<%s>>>"
		     (buffer-substring (point) end)))
	  (delete-region (point) end))))
    (while (re-search-forward "^public:\\s-*\n" (point-max) t)
      (delete-region (match-beginning 0) (match-end 0)))))

(defun swigpp-mask-nowrap-region ()
  (goto-char (point-min))
  (while (search-forward swigpp-nowrap-start (point-max) t)
    (let ((beg (progn (beginning-of-line) (point))))
      (search-forward swigpp-nowrap-end)
      (end-of-line) (forward-char 1)
      (when swigpp-debug
	(message "Masking nowrap region: <<<%s>>>"
		 (buffer-substring beg (point))))
      (delete-region beg (point)))))

(defun swigpp-mask-single-nowrap ()
  (goto-char (point-min))
  (let ((single-re (concat (regexp-quote swigpp-nowrap-one) "\\s-*$")))
    (while (re-search-forward single-re (point-max) t)
      (let ((beg (progn (beginning-of-line) (point)))
	    (end (cdr (swigpp-next-member))))
	(when swigpp-debug
	  (message "Masking single nowrap: <<<%s>>>"
		   (buffer-substring beg end)))
	(delete-region beg end)))))

(defun swigpp-mask-arg-default-values ()
  (goto-char (point-min))
  (while (and (search-forward "=" (point-max) t)
              (not (= ?= (char-before (1- (point)))))
              (not (string= "operator" (buffer-substring-no-properties
                                        (- (point) 10) (- (point) 2))))
              (not (string= "operator" (buffer-substring-no-properties
                                        (- (point) 9) (- (point) 1)))))
    (let ((beg (1- (point)))
          (end (save-excursion
                 (search-forward ";")
                 (1- (point)))))
      (cond
       ;; first look for constants
       ((re-search-forward "[^,()]+[,);]" end t)
        (delete-region beg (1- (point))))
       ;; then, simple (zero- or one-arg) function calls
       ((re-search-forward "[^,()]+([^,()]*)\\s-*[,);]" end t)
        (delete-region beg (1- (point))))
       ;; bail on complex function calls for now
       (t (error "swigpp-mask-arg-default-values: cannot handle <<<%S>>>"
                 (buffer-substring beg end)))))))

(defun swigpp-mask ()
  (let ((procs '(swigpp-subvert-definition-blocks
		 swigpp-mask-non-public
		 swigpp-mask-nowrap-region
		 swigpp-mask-single-nowrap
                 swigpp-mask-arg-default-values
		 ;; Add other masking functions here.
		 )))
    (mapc #'(lambda (proc)
              (funcall proc)
              (when swigpp-debug
                (message "[after %s] <<<%s>>>\n"
                         (symbol-name proc)
                         (buffer-string))))
          procs)))

;;;---------------------------------------------------------------------------
;;; Signature analysis

(defun swigpp-tokenize (b e)
  (save-excursion
    (goto-char b)
    (let ((me (cons 'w (intern (get 'swigpp 'name))))
	  state tokens snekot)
      ;; first go one direction
      (while (< (point) e)
	(if (or (looking-at "//.*") (looking-at "public:"))
	    (goto-char (match-end 0))
	  (setq state (char-to-string (char-syntax (char-after (point)))))
	  (let ((token (cons state
			     (buffer-substring
			      (point)
			      (+ (point) (skip-syntax-forward state e))))))
	    (unless (string-match "[-> ]" state)	;;; ignore whitespace
	      (push (if (string= "w" (car token))
			(cons (intern (car token)) (intern (cdr token)))
		      token)
		    snekot)))))
      ;; remove cruft
      (when (equal '("." . ";") (car snekot))
        (setq snekot (cdr snekot)))
      (when (equal '(w . const) (car snekot))
	(setq snekot (cdr snekot)))
      ;; now go the other direction
      (while snekot
	(let ((one (car snekot))
	      (two (cadr snekot)))
	  (cond ((and (string= "." (car one))
		      (equal '(w . operator) two))
		 (push (cons 'operator (cdr one)) tokens)
		 (setq snekot (cddr snekot)))
		((equal '("(" . "(") one)
		 (push 'L-paren tokens)
		 (setq snekot (cdr snekot)))
		((equal '(")" . ")") one)
		 (push 'R-paren tokens)
		 (setq snekot (cdr snekot)))
		((equal '("." . ",") one)
		 (push 'comma tokens)
		 (setq snekot (cdr snekot)))
		((and (equal '("." . "&") one)
		      (equal me two))
		 (push 'me-ref tokens)
		 (setq snekot (cddr snekot)))
                ((and (equal '("." . "::") one)
                      (equal me two))
                 (push 'me-class tokens)
                 (setq snekot (cddr snekot)))
		((equal me one)
		 (push 'me tokens)
		 (setq snekot (cdr snekot)))
		((equal '(w . const) one)
		 (push 'const tokens)
		 (setq snekot (cdr snekot)))
                ((equal '(w . inline) one)
                 (push 'inline tokens)
                 (setq snekot (cdr snekot)))
		(t
		 (push one tokens)
		 (setq snekot (cdr snekot))))))
      tokens)))

(defun swigpp-mangle-template-arg (ls)
  ;; This only works for simple foo<bar> templates for now.
  (if (and (listp ls)
           (listp (nth 0 ls))
           (eq 'w (car (nth 0 ls)))
           (equal '("." . "<") (nth 1 ls))
           (listp (nth 2 ls))
           (eq 'w (car (nth 2 ls)))
           (equal '("." . ">") (nth 3 ls)))
      (cons (cons 'w (make-symbol
                      (concat
                       (symbol-name (cdr (nth 0 ls)))
                       "_"
                       (symbol-name (cdr (nth 2 ls)))
		       "&")))
	    (nthcdr 4 ls))
    ls))

(defun swigpp-parse-list (tokens)
  (when swigpp-debug
    (message "swigpp-parse-list: %S" tokens))
  (mapcar 'swigpp-mangle-template-arg
          (loop with acc = nil
                for tok in tokens
                do (push tok acc)
                when (memq tok '(comma L-paren R-paren))
                collect (prog1
                            ;; If we wanted to make this system properly
                            ;; recursive, this value should be passed through
                            ;; `swigpp-parse-tree'.
                            (reverse (cdr acc))
                          (setq acc nil)))))

(defun swigpp-parse-c++-type (tokens)	; call cparse.el here?
  tokens)

(defun swigpp-parse-tree (tokens)
  (cond
   ;; simple data members
   ;; (data TOKEN...)
   ((not (memq 'L-paren tokens))
    (cons 'data tokens))
   ;; ctors (what about dtors?)
   ;; (ctor (ARG1-TOKENS...) [(ARG2-TOKENS...)...])
   ((and (eq 'me (car tokens))
	 (eq 'L-paren (cadr tokens)))
    (cons 'ctor (mapcar 'swigpp-parse-c++-type
			(swigpp-parse-list (cddr tokens)))))
   ;; operators (need to further distinguish)
   ;; (operator "TYPE" (ARG1-TOKENS...) [(ARG2-TOKENS...)...])
   ((and (memq (car tokens) '(me me-ref))
	 (listp (cadr tokens))
	 (eq 'operator (caadr tokens)))
    (append (list 'operator (cdadr tokens))
	    (mapcar 'swigpp-parse-c++-type
		    (swigpp-parse-list (cdddr tokens)))))
   ;; everything else is a function
   ;; (func PRE-L-PAREN-TOKENS... nil (ARG1-TOKENS...) [(ARG2-TOKENS...)...])
   (t
    (let (ret final-ret)
      (while tokens
        (if (eq 'L-paren (car tokens))
            (setq final-ret (append (reverse ret)
                                    '(nil) ; the nil neck!
                                    (mapcar 'swigpp-parse-c++-type
                                            (swigpp-parse-list
                                             (cdr tokens))))
                  tokens nil)
          (setq ret (cons (car tokens) ret)
                tokens (cdr tokens))))
      ;(message "[[[final-ret]]] %S" final-ret)
      (cons 'func final-ret)))))

(defun swigpp-note-signatures ()
  (goto-char (point-min))
  (loop for member = (swigpp-next-member)
	while member
	collect (let ((res (swigpp-parse-tree
			    (swigpp-tokenize (car member) (cdr member)))))
		  (when swigpp-debug
		    (message "Parse tree: %S" res))
		  res)))

;;;---------------------------------------------------------------------------
;;; Mapping
;;;
;;; TODO: Factor `swigpp-map-class-operator' and
;;; `swigpp-map-extra-class-operator'.

(defun swigpp-infixify (ls)
  (let ((res '("(")))
    (while ls
      (let ((elem (car ls)))
	(if (listp elem)		; flatten only one level
	    (mapc (lambda (sub-elems)
		    (push sub-elems res))
		  elem)
	  (push elem res)))
      (setq ls (cdr ls))
      (when ls (push "," res)))
    (reverse (cons ")" res))))

;; TODO: Redesign.
(defun swigpp-map-ctor (ctors)
  ;; Map the first non-copying non-default constructor as *the* ctor.
  ;; Ignore subsequent ones.  If there are none, use default.
  (let (default copying the-ctor)
    (when (member '(nil) ctors)
      (when swigpp-debug (message "Found a default ctor"))
      (setq ; ctors (remove '(nil) ctors)
	    default '(nil)))
    (let ((try (find-if (lambda (formals)
			  (and (= 1 (length formals))
			       (memq 'me-ref (car formals))))
			ctors)))
      (when try
	(when swigpp-debug (message "Found a copying ctor"))
	(setq ctors (remove try ctors)
	      copying try)))
    (when (and swigpp-debug (> (length ctors) 0))
      (message "Ignoring %d ctors besides first" (1- (length ctors))))
;    (setq the-ctor (if (null ctors)
;		       (or default copying)
;		     (car ctors)))
;    `(sig me ,@(swigpp-infixify the-ctor) ";")
    (when ctors `(sig me ,@(swigpp-infixify (car ctors)) ";"))))

(defun swigpp-map-class-operator (operator)
  (when swigpp-debug (message "swigpp-map-class-operator: %S" operator))
  (let ((op (car operator))
	(formals (cdr operator)))
    (flet ((var-of (type-spec) (car (reverse type-spec)))
	   (binary-op (type func)
                      (and (= 1 (length formals))
                           (string= type op)
                           (not (string= "=" (substring op -1)))
                           `(sig me ,func ,@(swigpp-infixify formals) const
                                 "{return " ,(concat "self->operator" type)
                                 "(" ,(var-of (car formals)) ");}"))))
      (cond
       ;; unary minus
       ((and (equal '(nil) formals) (string= "-" op))
	`(sig me "__neg__ () const {return self->operator-();}"))
       ((binary-op "+" "__add__"))
       ((binary-op "-" "__sub__"))
       ((binary-op "*" "__mul__"))
       ((binary-op "/" "__div__"))
       ;; Add new operator rules here.
       ))))

(defun swigpp-map-extra-class-operator (operator)
  (when swigpp-debug (message "swigpp-map-extra-class-operator: %S" operator))
  (let ((op (car operator))
        (formals (cdr operator)))
    (when (= 2 (length formals))
      (setq formals (cdr formals))
      (flet ((var-of (type-spec) (car (reverse type-spec)))
             (binary-op (type func)
                        (and (string= type op)
                             (not (string= "=" (substring op -1)))
                             `(sig me ,func ,@(swigpp-infixify formals) const
                                   "{" me " tmp = (*self) " ,type
                                   ;; ,(var-of (car formals))
                                   ;;,(concat "(*self).operator" type)
                                   ,(var-of (car formals)) ";"
                                   " return tmp;}"))))
        (cond
         ((binary-op "+" "__add__"))
         ((binary-op "-" "__sub__"))
         ((binary-op "*" "__mul__"))
         ((binary-op "/" "__div__"))
         ;; Add new operator rules here.
         )))))

(defun swigpp-map-members (sigs)
  (when swigpp-debug
    (message "swigpp-map-members: %S" sigs))
  (flet ((filter (type) (mapcar 'cdr
				(remove-if-not (lambda (x)
						 (eq x type))
					       sigs
					       :key 'car))))
    (let ((data  (filter 'data))
	  (ctors (filter 'ctor))
	  (ops   (filter 'operator))
	  (funcs (filter 'func))
          class except add-methods inline top-level)

      (when swigpp-debug
	(message "data:  %S" data)
	(message "ctors: %S" ctors)
	(message "ops:   %S" ops)
	(message "funcs: %S" funcs))

      ;; Constructors
      (let ((bunch-of-ctors (swigpp-map-ctor ctors)))
        (when bunch-of-ctors (push bunch-of-ctors class)))

      ;; Operators
      (let (am)
	(mapc (lambda (op)
		(push (or (swigpp-map-class-operator op)
			  `(note ignoring operator: ,@op))
		      am))
	      ops)
	(push `(add-method ,@am) class))

      ;; Functions
      (mapc (lambda (func)
              (push (cons 'func func) class))
            funcs)

      ;; Data members
      ;;
      ;; These are passed through unmodifed.
      (mapc (lambda (datum)
              (push `(data ,@datum ";") class))
            data)

      ;; Return things in the proper order.
      (remove nil
              (mapcar (lambda (var)
                        (let ((val (eval var)))
                          (and val (cons var (reverse val)))))
                      '(class except add-methods inline top-level))))))

(defun swigpp-map-extra-class-methods (assumed-funcs)
  (when swigpp-debug
    (message "ASSUMED-FUNCS: %S" assumed-funcs))
  (let (operator)
    (loop for func in assumed-funcs
          ;; We only consider simple class-specific operators at this time.
          when (and (eq 'func (car func))
                    (let* ((fguts (cdr func))
                           (ret-name (car fguts))
                           (formals (cdr fguts)))
                      (when (eq 'inline (car-safe ret-name))
                        (setq ret-name (cdr ret-name)))
                      (and (memq (car-safe ret-name) '(me me-ref))
                           (listp (cadr ret-name))
                           (eq 'operator (caadr ret-name))
                           ;; this assignment is non-optimal
                           (setq operator (cons (cdadr ret-name) formals)))))
          collect ;(cons 'func
                        (swigpp-map-extra-class-operator operator)
                        ;)
          )))

;;;---------------------------------------------------------------------------
;;; Output

(defun swigpp-obuf () (get 'swigpp 'obuf))

(defun swigpp-insert (&rest s)
  (with-current-buffer (swigpp-obuf)
    (apply 'insert s)))

(defun swigpp-untokenize (token)
  (flet ((tagged (tag) (and (listp token) (eq tag (car token)))))
    (let ((name (get 'swigpp 'name)))
       (cond ((stringp token) token)
	    ((tagged 'w) (symbol-name (cdr token)))
            ((tagged 'operator) (concat "operator" (cdr token)))
            ((and (listp token) (equal "." (car token))) (cdr token))
            ((case token
               ((me) name)
               ((me-ref) (concat name "&"))
	       ((const) "const")
	       ((comma) ",")))		; todo: delete
            (t (error "swigpp-untokenize: Unrecognized token: %S" token))))))

(defun swigpp-render-note (note)
  (swigpp-insert (format "// %S\n" note)))

(defun swigpp-render-signature (sig)
  (when swigpp-debug (message "swigpp-render-signature: %S" sig))
  (swigpp-insert (mapconcat 'swigpp-untokenize sig " ") "\n"))

(defun swigpp-render-data (data)
  (when swigpp-debug (message "swigpp-render-data: %S" data))
  (swigpp-insert (mapconcat 'swigpp-untokenize data " ") "\n"))

(defun swigpp-render-class (class-tree)
  (when swigpp-debug (message "swigpp-render-class: %S" class-tree))
  (let ((name (get 'swigpp 'name))
        (inher (get 'swigpp 'inher)))
    (swigpp-insert "class " name inher "{\npublic:\n")
    (mapc 'swigpp-render class-tree)
    (swigpp-insert "}; // " name "\n\n")))

(defun swigpp-render-add-method (add-method &optional name)
  (when swigpp-debug (message "swigpp-render-add-method: %S" add-method))
  (swigpp-insert "%addmethods " (or name "") " {\n")
  (mapc 'swigpp-render add-method)
  (swigpp-insert "}\n"))

(defun swigpp-render-func (func)
  (when swigpp-debug (message "swigpp-render-func: %S" func))
  (while (car func)
    (swigpp-insert (swigpp-untokenize (car func)) " ")
    (setq func (cdr func)))
  (setq func (cdr func))                ; ignore neck
  (swigpp-insert "(")
  (swigpp-insert (mapconcat #'(lambda (formal-arg)
                                (mapconcat 'swigpp-untokenize
                                           formal-arg
                                           " "))
                            func
                            ", "))
  (swigpp-insert ");\n"))

(defun swigpp-render (tree)
  (when swigpp-debug (message "swigpp-render: %S" tree))
  (if (listp (car tree))
      (mapc 'swigpp-render tree)
    (funcall (case (car tree)
               ((class)       'swigpp-render-class)
               ((sig)         'swigpp-render-signature)
               ((data)        'swigpp-render-data)
               ((add-method)  'swigpp-render-add-method)
	       ((note)	      'swigpp-render-note)
               ((func)        'swigpp-render-func)
               (t (error "swigpp-render: Unknown tree type: %S" (car tree))))
             (cdr tree))))

;; TODO: Redesign.
(defun swigpp-process-class (workbuf name inher b e)
  (save-excursion
    (let ((srcbuf (current-buffer)))
      (append-to-buffer workbuf b e)
      (put 'swigpp 'name name)
      (put 'swigpp 'inher inher)
      ;; first do members
      (set-buffer workbuf)
      (swigpp-mask)
      (swigpp-insert "class " name " " inher " {\npublic:\n")
      (swigpp-insert (buffer-string))
      (swigpp-insert "}; // " name "\n\n")
      (erase-buffer)
;      ;; now look at rest of file
;      (set-buffer srcbuf)
;      (append-to-buffer workbuf (+ 2 e) (point-max)) ; hmm 2
;      (set-buffer workbuf)
;      (swigpp-mask)
;      (when swigpp-debug (swigpp-message-buffer))
;      (swigpp-render-add-method
;       (swigpp-map-extra-class-methods (swigpp-note-signatures))
;       (get 'swigpp 'name))
;      (erase-buffer)
      )))

;;;---------------------------------------------------------------------------
;;; Entry points

(defun swigpp-scrub-input-file (file)
  (set-buffer (generate-new-buffer "*swig ibuf*"))
  (setq buffer-file-name file)          ; ugh
  (insert-file-contents-literally file)
  (goto-char (point-min))
  (flush-lines "^#\\s-*include")
  (let (swigpp-directives)
    (save-excursion
      (while (re-search-forward "//!.*" (point-max) t)
        (push (cons (count-lines (point-min) (point)) (match-string 0))
              swigpp-directives)))
    (goto-char (point-min))
    (put 'swigpp 'lit
         (if (re-search-forward "^#if SWIGPP_LITERAL_INCLUDE.*\n"
                                (point-max) t)
             (let ((p (point)))
               (search-forward "#endif // SWIGPP_LITERAL_INCLUDE")
               (concat "/* SWIGPP_LITERAL_INCLUDE starts here */\n"
                       (buffer-substring p (match-beginning 0))
                       "/* SWIGPP_LITERAL_INCLUDE ends here */\n\n"))
           ""))
    (shell-command-on-region (point-min) (point-max) "g++ -E -" nil t)
    (goto-char (point-min))
    (flush-lines "^#")
    (mapc (lambda (note)
            (goto-line (car note))
            (insert (cdr note)))
          swigpp-directives)
    (when swigpp-debug (swigpp-message-buffer))
    (current-buffer)))

(defun swigpp-prettify-output-buffer (buf)
  (set-buffer buf)
  (goto-char (point-min))
  (c++-mode)
  (while (re-search-forward "^class \\S-+\\s-*" (point-max) 1)
    (c-indent-exp)
    (forward-sexp 1)))

(defun path-to-include-name (opath)
  ;; added by bashford
  ;; turn a string like "../../mead/libmead/PhysCond.h"
  ;; into "MEAD/PhysCond.h"
  ;; (Assumes that any  "*/libmead/ part will be accessible as MEAD
  ;; due to suitable -I flags to the compiler)
  (if (string-match "^.*/libmead/\\([^/]+$\\)" opath)
      (concat "MEAD/" 
	      (substring opath (match-beginning 1) (match-end 1)))
    opath))

      
(defun* swigpp-process-file (inf outf module-name)
  (let (c++-mode-hook find-file-hooks hs-minor-mode-hook ibuf classes)
    (put 'swigpp 'obuf (find-file outf))
    (erase-buffer)
    (insert "// Do not edit -- autogenerated file.\n"
	    "// Written " (format-time-string "%c")
	    " by swigpp.el $Revision: 1.61 $\n\n"
            "#ifndef _SWIGPP_" module-name "_j\n"
            "#define _SWIGPP_" module-name "_j\n\n"
	    "%module " module-name "\n"
;;            "%include ./meadtypes.i\n\n"
	    "%{\n#include \"" (path-to-include-name inf) "\"\n%}\n\n")
    (setq ibuf (swigpp-scrub-input-file inf))

    ;; process classes
    (let ((class-id-re (concat (regexp-quote swigpp-wrap-one)
			       ".*\n\\(class\\|struct\\)"
                               "\\s-+\\([a-zA-Z0-9_]+\\)"
                               "\\(\\s-*\\([^{]\\|\n\\)*\\){")))
      (goto-char (point-min))
      (while (re-search-forward class-id-re (point-max) t)
        (forward-char -1)
	(push (list (match-string 2)
                    (let ((inher (match-string 3)))
                      (if (string-match "<" inher) ; ignore templates for now
                          ""
                        inher))
                    (1+ (point))
		    (progn (forward-sexp 1) (1- (point))))
	      classes))
      (unless classes
	(message "No classes found, not writing %s" outf)
	(return-from swigpp-process-file))
      (setq classes (nreverse classes))
      (message "Found classes: %s" (mapconcat 'car classes " "))
      (let ((workbuf (generate-new-buffer (concat ".swigpp-" module-name))))
        (save-excursion
          (set-buffer workbuf)
          (let (c++-mode-hook)          ; ugh
            (c++-mode)
            (modify-syntax-entry ?_ "w")))
        (mapc (lambda (class)
                (when swigpp-debug (message "***\n*** %S\n***" class))
                (apply 'swigpp-process-class workbuf class))
              classes)
        (when swigpp-debug
          (set-buffer workbuf)
          (write-file (buffer-name)))
        (kill-buffer workbuf)))

    ;; literal text
    (swigpp-insert (get 'swigpp 'lit))

    ;; cleanup
    (swigpp-prettify-output-buffer (swigpp-obuf))
    (with-current-buffer (swigpp-obuf)
      (goto-char (point-max))
      (insert "#endif /* !_SWIGPP_" module-name "_j */\n")
      (save-buffer))
    (kill-buffer (swigpp-obuf))
    (kill-buffer ibuf)
    (setplist 'swigpp nil)))

(defun batch-swigpp-process-file ()
  (apply 'swigpp-process-file command-line-args-left)
  (message "batch-swigpp-process-file: Done %s" (car command-line-args-left))
  ;; This should not be required, but is included here to workaround
  ;; a bug in GNU Emacs 20.3.1 (i386-redhat-linux-gnu, X toolkit)
  ;; of Mon Apr 19 1999 on porky.devel.redhat.com.
  (kill-emacs))

;;;---------------------------------------------------------------------------
;;; Load-time actions

(provide 'swigpp)

;;; swigpp.el ends here
