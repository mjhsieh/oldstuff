;;; gen.el --- various text-instantiation commands
;;; Author: Thien-Thi Nguyen

;;; Commentary:

;; Following the code are the data sections, one for each template-type.
;; Currently there are three supported: list, vector, map.  Each data section
;; begins with the comment ";;; Data: TEMPLATE-TYPE" followed by a blank line
;; followed by lines beginning with ";" followed by another blank line.

;;; Code:

(defun gen-subst (data-type &rest subs)
  (save-excursion
    (set-buffer (generate-new-buffer "*gen*"))
    (erase-buffer)
    (insert-file-contents "./gen.el")
    (re-search-forward (concat "^;;; Data: " data-type "\n\n"))
    (delete-region (point-min) (point))
    (while (looking-at "^;")
      (delete-char 1) (end-of-line) (forward-char 1))
    (delete-region (point) (point-max))
    (let (c-types)
      (while subs
        (goto-char (point-min))
        (let ((c-type (cdar subs)))
          (setq c-types (cons c-type c-types))
          (while (search-forward (concat "@" (caar subs) "@") (point-max) t)
            (replace-match c-type t)))
        (setq subs (cdr subs)))
      (goto-char (point-min))
      (insert "%{\n#include <" data-type ">\n")
      (while c-types
        (let ((c-type (car c-types)))
          (unless (memq c-type (mapcar 'symbol-name
                                       '(int float double char)))
            (insert "#include \""
                    ;; if c-type starts w/ a capital letter, prefix "MEAD/"
                    (let ((first-char (substring c-type 0 1)))
                      (if (string= first-char (downcase first-char))
                          ""
                        "MEAD/"))
                    c-type
                    ".h\"\n")))
        (setq c-types (cdr c-types)))
      (insert "%}\n\n"))
    (prog1
        (buffer-string)
      (kill-buffer (current-buffer)))))

(defun gen-from-list (type)
  (when (symbolp type) (setq type (symbol-name type)))
  (insert (gen-subst "list" (cons "TYPE" type))))

(defun gen-from-vector (type)
  (when (symbolp type) (setq type (symbol-name type)))
  (insert (gen-subst "vector" (cons "TYPE" type))))

(defun gen-from-map (type-a type-b)
  (when (symbolp type-a) (setq type-a (symbol-name type-a)))
  (when (symbolp type-b) (setq type-b (symbol-name type-b)))
  (insert (gen-subst "map" (cons "TYPE-A" type-a) (cons "TYPE-B" type-b))))

(defun gen-meadtypes (from-list from-vector from-map)
  (mapcar #'gen-from-list from-list)
  (mapcar #'gen-from-vector from-vector)
  (mapcar #'(lambda (args) (apply #'gen-from-map args)) from-map))

;;; Data: list

;%{
;typedef std::list<@TYPE@> list_@TYPE@;
;%}
;
;%typemap(python,in) list_@TYPE@& ($basetype vi) {
;  if (PyList_Check($source)) {
;    int size = PyList_Size($source);
;    for (int i=0; i<size; ++i) {
;      PyObject *o = PyList_GetItem($source, i);
;      @TYPE@* a=0;
;      if ((SWIG_ConvertPtr(o,(void **) &a,SWIGTYPE_@TYPE@_p,1)) == -1) {
;	PyErr_SetString(PyExc_TypeError,
;			"Conversion of list item to @TYPE@ failed");
;	return NULL;
;      }
;      vi.push_back(*a);
;    }
;    $target = &vi;
;  }
;  else {
;    PyErr_SetString(PyExc_TypeError, "not a list");
;    return NULL;
;  }
;}
;
;%typemap(python,argout) list_@TYPE@& {
;  static PyObject *@TYPE@Ptr_class = 0;
;  if (!@TYPE@Ptr_class) {
;    //    cerr << "Setting up the @TYPE@Ptr_class pointer" << endl;
;    // The call below returns only a "borrowed" ref, which is good (I think)
;    if(PyObject *module = PyImport_ImportModule("Triangle"))
;      // Note: the call below returns only a "borrowed" ref,
;      // but I think that's okay, since of the "holding" module
;      // of @TYPE@Ptr deletes it, we're screwed anayway.
;      @TYPE@Ptr_class = PyDict_GetItemString(PyModule_GetDict(module),
;					"@TYPE@Ptr");
;    else
;      return NULL;
;    if (!PyCallable_Check(@TYPE@Ptr_class)) {
;      PyErr_SetString(PyExc_TypeError, "@TYPE@Ptr_class is not callable");
;      return NULL;
;    }
;  }
;  PyObject *newlist = PyList_New(0);
;  if (!newlist) return NULL;
;  for (list_@TYPE@::const_iterator p = $source->begin(); p!=$source->end(); ++p) {
;    @TYPE@ * fp = new @TYPE@(*p);
;    PyObject *o = SWIG_NewPointerObj((void *)fp, SWIGTYPE_@TYPE@_p);
;    if (@TYPE@Ptr_class) {
;      //      cerr << "building the arglist" << endl;
;      PyObject *arglist = Py_BuildValue("(O)", o);
;      PyObject *optr = PyEval_CallObject(@TYPE@Ptr_class, arglist);
;      Py_DECREF(arglist);
;      PyList_Append(newlist, optr);
;      Py_DECREF(optr);
;    }
;    else {
;      PyList_Append(newlist, o);
;    }
;    Py_XDECREF(o);
;  }
;  PyList_SetSlice($arg, 0, PyList_Size($arg), newlist);
;  Py_DECREF(newlist);
;}
;
;%typemap(python,out) list_@TYPE@* {
;  static PyObject *@TYPE@Ptr_class = 0;
;  if (!@TYPE@Ptr_class) {
;    //    cerr << "Setting up the @TYPE@Ptr_class pointer" << endl;
;    // The call below returns only a "borrowed" ref, which is good (I think)
;    if(PyObject *module = PyImport_ImportModule("Triangle"))
;      // Note: the call below returns only a "borrowed" ref,
;      // but I think that's okay, since of the "holding" module
;      // of @TYPE@Ptr deletes it, we're screwed anayway.
;      @TYPE@Ptr_class = PyDict_GetItemString(PyModule_GetDict(module),
;					"@TYPE@Ptr");
;    else
;      return NULL;
;    if (!PyCallable_Check(@TYPE@Ptr_class)) {
;      PyErr_SetString(PyExc_TypeError, "@TYPE@Ptr_class is not callable");
;      return NULL;
;    }
;  }
;  $target = PyList_New(0);
;  if (!$target) return NULL;
;  for (list_@TYPE@::const_iterator p = $source->begin(); p!=$source->end(); ++p) {
;    @TYPE@ * fp = new @TYPE@(*p);
;    PyObject *o = SWIG_NewPointerObj((void *)fp, SWIGTYPE_@TYPE@_p);
;    if (@TYPE@Ptr_class) {
;      //      cerr << "building the arglist" << endl;
;      PyObject *arglist = Py_BuildValue("(O)", o);
;      PyObject *optr = PyEval_CallObject(@TYPE@Ptr_class, arglist);
;      Py_DECREF(arglist);
;      PyList_Append($target, optr);
;      Py_DECREF(optr);
;    }
;    else {
;      PyList_Append($target, o);
;    }
;    Py_XDECREF(o);
;  }
;  delete $source;  // Or should it be in a "clean up" typemap? FIXME?
;}
;

;;; Data: vector

;%{
;typedef std::vector<@TYPE@> vector_@TYPE@;
;%}
;
;%typemap(python,in) vector_@TYPE@& ($basetype vpe) {
;  if (PyList_Check($source)) {
;    int size = PyList_Size($source);
;    for (int i=0; i<size; ++i) {
;      PyObject *o = PyList_GetItem($source, i);
;      @TYPE@ *peptr;
;      if ((SWIG_ConvertPtr(o,(void **) & peptr,
;			   SWIGTYPE_@TYPE@_p,1)) != -1) {
;	// if here, the list item, o was succesfully translated
;	vpe.push_back(*peptr);
;      }
;      else {
;	PyErr_SetString(PyExc_TypeError, "list must contain @TYPE@ items");
;	return NULL;
;      }
;    }
;    $target = &vpe;
;  }
;  else {
;    PyErr_SetString(PyExc_TypeError, "not a list");
;    return NULL;
;  }
;}
;

;;; Data: map

;/* TODO: map<@TYPE-A@,@TYPE-B@> */
;

;;; gen.el ends here
