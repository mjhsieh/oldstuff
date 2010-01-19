*** MEAD_shadow.cc	Thu Jun 21 11:09:44 2001
--- MEAD_shadow.cc.orig	Thu Jun 21 11:25:42 2001
*************** SWIG_TypeCheck(char *c, swig_type_info *
*** 129,134 ****
--- 129,135 ----
        ty->next = s;
        return s;
      }
+     if (s == s->next) break;
      s = s->next;
    }
    return 0;
