#!/usr/bin/env perl
if(scalar @ARGV==0 || $ARGV[0] =~ /^-h/oi){
   print STDERR "Usage: perl $0 your_interface_file > new_interface_file\n";
   exit;
}
$refgone="&";
$topointer="*";
if($ARGV[0] eq "-ref"){
   $refgone="";
   $topointer="";
   shift @ARGV;
}
%binopmap=(
'<'=>'__lt__', '>'=>'__gt__', '=='=>'__eq__', '!='=>'__ne__',
'<='=>'__le__', '>='=>'__ge__',
'+'=>'__add__', '-'=>'__sub__','*'=>'__mul__', '/'=>'__div__',
'&'=>'__and__', '^'=>'__xor__','|'=>'__or__', '%'=>'__mod__',
'+='=>'__iadd__','-='=>'__isub__','*='=>'__imul__','/='=>'__idiv__',
'&='=>'__iand__','^='=>'__ixor__','|='=>'__ior__','%='=>'__imod__',
'<<'=>'__lshift__','>>'=>'__rshift__',
'[]'=>'__getitem__','()'=>'__call__'
);
%unopmap=( '-' =>'__neg__','+' =>'__pos__');
# These can be rop's
%ropmap=('+'=>'__radd__', '-'=>'__rsub__','*'=>'__rmul__', '/'=>'__rdiv__');
# Count the occurrences of the binary operators
%binopfreq=(
'__lt__', 0, '__gt__', 0, '__eq__', 0, '__ne__', 0, '__le__', 0, '__ge__', 0,
'__add__', 0, '__sub__', 0, '__mul__', 0, '__div__', 0,
'__radd__', 0, '__rsub__', 0, '__rmul__', 0, '__rdiv__', 0,
'__and__', 0, '__xor__', 0, '__or__', 0, '__mod__', 0,
'__iadd__', 0, '__isub__', 0, '__imul__', 0, '__idiv__', 0,
'__iand__', 0, '__ixor__', 0, '__ior__', 0, '__imod__', 0,
'__lshift__', 0, '__rshift__', 0,
'__call__', 0);
# Count the occurrences of the unary operators
%unopfreq=('-', 0, '+', 0);
# Count the occurrences of the operators necessary to implement __cmp__
%cmpfreq=('==', 0, '<', 0, '>', 0);

$lastline="";
$foundendif=0;

$ifile=shift @ARGV;

# First count the occurrences of the binary operators for each class
$class="";
$lastclass="";
$num_classes=0;
open(IFILE,$ifile);
while(<IFILE>){
   if(m|^\s*//|o){next;}
   $_ =~ s/\/\/.*$//o;
   $foundkey=0;
   for $c(keys %binopfreq){
      if(/\b$c\s*\(/o){$binopfreq{$c}++; $foundkey=1; last;} }
   if($foundkey){next;}
   elsif(/class\s+(\w+)/o){
      if($class ne ""){
         for $c(keys %binopfreq){$binopcount{$class}{$c} = $binopfreq{$c};}}
      $lastclass=$class;
      $class=$1;
      $classnames[$num_classes]=$class;
      $num_classes++;
      for $c(keys %binopfreq){$binopfreq{$c}=0;}
   }elsif(/%addmethods\s+(\w+)/o){
      if($1 ne $class){
         $lastclass=$class;
         $class=$1;
         for $c(keys %binopfreq){$binopcount{$lastclass}{$c} += $binopfreq{$c};}
         %binopfreq=%{$binopcount{$class}};}
   }elsif(/\boperator\b\s*[^\w\s]+\s*\(/o){
      $s=$_; $h=$_;
      if($_ !~ /^(.+?)operator\s*(\W+?)\s*\(/o || !exists $binopmap{$2}){
         next;
      }
      $rettype=$1; $op=$2;

      if($_ !~ /;/o && $_ !~ /\{/o){
         while(<IFILE>){
            $s.=$_;
            if(/;/o){ last;}
         }
      }elsif($_ =~ /\{/o){
         $p=1;
         while(<IFILE>){
            $s.=$_;
            s/([\{\}])/if($1 eq "{"){$p++}else{$p--} $1/goe;
            if($p==0 && /\}/o){ last;}
         }
      }

# Handle the two parameter case
      $s=~ /\(\s*(.*,)?\s*(.*)\s*\)/mo;
      $param1=$1;
      $param2=$2;
      $param=$param2;
      $param1=~ /^(.*)\s*,/o;
      $param1=$1;
#      print "Operator: $op \n";
#      print "First parameter: $param1 \n";
#      print "Second parameter: $param2 \n";

      if($param =~ s/^(\w+)$/$1 oprnd2/o){
      }elsif($param =~ s/((const\s+)?\w+\s*[^\w\s]*)\s*\w+$/$1 oprnd2/o){
      }elsif($param =~ s/((const\s+)?\w+\s*[^\w\s]*)$/$1 oprnd2/o){ }
#      print "Param: $param \n";

      if($param ne ""){
# Handle the rop case where the first param is not type $class
#         print "param eq \"\" \n";
         if($param1 ne ""){
            if (exists $ropmap{$op} && index($param1, "$class") < 0){
               $__op__ = $ropmap{$op};
            }else{
               $__op__ = $binopmap{$op};
               if (exists $binopfreq{$__op__}){$binopfreq{$__op__}++;}
            }
         }else{
            $__op__ = $binopmap{$op};
            if (exists $binopfreq{$__op__}){$binopfreq{$__op__}++;}
         }
      }else{
#         print "param ne \"\" \n";
         $__op__ = $unopmap{$op};
      }
   }
}
close(IFILE);

if($class ne ""){
   for $c(keys %binopfreq){$binopcount{$class}{$c} = $binopfreq{$c};}
}

#print "The following counts were found\n";
#foreach $elem (@classnames){
#   print "For class $elem \n";
#   print "Binaray Op's \n";
#   %binopfreq= %{$binopcount{$elem}};
#   foreach $op (keys %binopfreq){
#      print "$op $binopfreq{$op}\n";}
#}

# Now process the operators
$class="";
$lastclass="";
open(IFILE,$ifile);
$cmtout="//";
while(<IFILE>){
  if (/#endif/o){
      $foundendif = 1;
      $lastline = $_;
      next;
   }elsif($foundendif){
      print $lastline;
      $foundendif = 0;
   }
   if(m|^\s*//|o){ print; next;}
   elsif(/\%inline/o){ $cmtout="";}
   elsif(/\%\}/o){ $cmtout="//";}
   elsif(/class\s+(\w+)/o){
      if($class ne ""){
         for $c(keys %unopfreq) {$unopcount{$class}{$c} = $unopfreq{$c};}
         for $c(keys %cmpfreq) {$cmpcount{$class}{$c} = $cmpfreq{$c};}}
      $lastclass=$class;
      $class=$1;
      %binopfreq=%{$binopcount{$class}};
      for $c(keys %unopfreq){$unopfreq{$c}=0;}
      for $c(keys %cmpfreq){$cmpfreq{$c}=0;}
   }elsif(/%addmethods\s+(\w+)/o){
      if($1 ne $class){
         $lastclass=$class;
         $class=$1;
         for $c(keys %unopfreq) {$unopcount{$lastclass}{$c} += $unopfreq{$c};}
         for $c(keys %cmpfreq) {$cmpcount{$lastclass}{$c} += $cmpfreq{$c};}
         %binopfreq=%{$binopcount{$class}};
         %unopfreq=%{$unopcount{$class}};
         %cmpfreq=%{$cmpcount{$class}};}
   }elsif(/\boperator\b\s*[^\w\s]+\s*\(/o){
      if(($comment = index($_, "//")) > 0){
        $s = substr($_, 0, $comment); $h=$s;
      }else{
         $s=$_; $h=$_;
      }
      if($_ !~ /^(.+?)operator\s*(\W+?)\s*\(/o || !exists $binopmap{$2}){
         print STDERR "WARNING: Operator overloading at Line $. cannot be translated... ignored\n   =>$_";
         print $_; next;
      }
      $rettype=$1; $op=$2;
      print "$cmtout$_";
      if($_ !~ /;/o && $_ !~ /\{/o){
         while(<IFILE>){
            print "$cmtout$_";
            $s.=$_;
            if(/;/o){ last;}
         }
      }elsif($_ =~ /\{/o){
         $p=1;
         while(<IFILE>){
            print "$cmtout$_";
            $s.=$_;
            s/([\{\}])/if($1 eq "{"){$p++}else{$p--} $1/goe;
            if($p==0 && /\}/o){ last;}
         }
      }
#      $s=~ /\(\s*(.*)\s*\)/mo;
#      $param=$1;

# Handle the two parameter case
      $s=~ /\(\s*(.*,)?\s*(.*)\s*\)/mo;
      $param1=$1;
      $param2=$2;
      $param=$param2;
      $param1=~ /^(.*)\s*,/o;
      $param1=$1;
#      print "Operator: $op \n";
#      print "First parameter: $param1 \n";
#      print "Second parameter: $param2 \n";

      if($param =~ s/^(\w+)$/$1 oprnd2/o){
      }elsif($param =~ s/((const\s+)?\w+\s*[^\w\s]*)\s*\w+$/$1 oprnd2/o){
      }elsif($param =~ s/((const\s+)?\w+\s*[^\w\s]*)$/$1 oprnd2/o){
      }
#      print "Param: $param \n";

      if($param ne ""){
# Handle the rop case where the first param is not type $class
#         print "param ne \"\" \n";
         if($param1 ne ""){
            if(exists $ropmap{$op} && index($param1, "$class") < 0){
               $__op__ = $ropmap{$op};
#               print "is an rop $__op__ \n";
            }else{
               $__op__ = $binopmap{$op};
#               print "is an binop $__op__ \n";
            }
         }else{
            $__op__ = $binopmap{$op};
#            print "is an binop $__op__ \n";
         }
      }else{
#         print "param eq \"\" \n";
         $__op__ = $unopmap{$op};
#         print "is a unop $__op__ \n";
      }
# See if we are defining multiple operators of the same kind
      if (exists $binopfreq{$__op__} && $binopfreq{$__op__} > 1){
         $ismultiple=1;
         $ending=";";
         if(!exists $wrapper{$class}){
            $wrapper{$class}="%wrapper %{\n";
            $wrapper{$class}.="#ifdef __cplusplus\n}\n#endif\n\n";
         }
      }else{
         $ismultiple=0;
         $ending="{";
      }
      if(!exists $addmethods{$class}){
         $addmethods{$class}="%addmethods $class\{\n";
      }
      if(exists $cmpfreq{$op}){$cmpfreq{$op}++;}
      $w="";
      if($op eq "[]"){
         $h =~ s/^\s*(.+?)\s*$refgone\s*operator\s*(\W+?)\s*\((.*)\).*$/   $1 $topointer $binopmap{$op}\($param){/o;
         $addmethods{$class}.=$h."     return ${refgone}(*self)[oprnd2];\n   }\n";
         $rettype =~ s/$refgone//o;
         $h =~ s/^.*?__getitem__\s*\(.*?\)/   void __setitem__($param,$rettype _value)/o;
         $h .= "     (*self)[oprnd2]=_value;\n   }\n";
      }elsif($op eq "()"){
         $param=~s/\)\(//o;
         if($param ne ""){
            if($param1 ne ""){
               if($param1 =~ s/^(\w+)$/$1 oprnd2/o){
               }elsif($param1 =~ s/((const\s+)?\w+\s*[^\w\s]*)\s*\w+$/$1 oprnd2/o){
               }elsif($param1 =~ s/((const\s+)?\w+\s*[^\w\s]*)$/$1 oprnd2/o){ }
               $param1=~s/\)\(//o;
               $param=~s/oprnd2/oprnd3/o;
               $params="${param1}, ${param}";
               $oprnds="oprnd2, oprnd3";
            }else{
               $params=${param};
               $oprnds="oprnd2";
            }
         }else{
            $params="";
            $oprnds="";
         }
         if($ismultiple){
            $h =~ s/^\s*(.+?)\s*operator\s*(\W+?)\s*\((.*)\).*$/   $1 $binopmap{$op}\($params);/o;
            $rettype =~ s/^\s*//o;
            $w = "$rettype ${class}_$__op__(${class} *self, ${params}) {\n";
            $w .= "   return self->operator$op ($oprnds);\n}\n";
         }else{
            $h =~ s/^\s*(.+?)\s*operator\s*(\W+?)\s*\((.*)\).*$/   $1 $binopmap{$op}\($params){/o;
            $h .= "     return (*self)($oprnds);\n   }\n";
         }
      }elsif($param eq ""){
         $unopfreq{$op}++;
         $h =~ s/^\s*(.+?)operator\s*(\W+?)\s*\((.*)\).*$/   $1 $unopmap{$op}\(){/o;
         $h .= "     return $op(*self);\n   }\n";
      }else{
# Handle the rop case where the first param is not type $class
         if($param1 ne "" && exists $ropmap{$op} && index($param1, "$class") < 0){
            if($param1 =~ s/^(\w+)$/$1 oprnd2/o){
            }elsif($param1 =~ s/((const\s+)?\w+\s*[^\w\s]*)\s*\w+$/$1 oprnd2/o){
            }elsif($param1 =~ s/((const\s+)?\w+\s*[^\w\s]*)$/$1 oprnd2/o){ }
            $h =~ s/^\s*(.+?)operator\s*(\W+?)\s*\((.*)\).*$/   $1 $ropmap{$op}\($param1)$ending/o;
         }else{
            $h =~ s/^\s*(.+?)operator\s*(\W+?)\s*\((.*)\).*$/   $1 $binopmap{$op}\($param)$ending/o;
         }
         $rettype =~ s/^\s*//o;
         if($ismultiple){
            if ($param1 eq ""){
               $w = "$rettype ${class}_$__op__(${class} *self, ${param}) {\n";
               $w .= "   return self->operator$op (oprnd2);\n}\n";
            }else{
               $rettype =~ s/$refgone//o;
               if(exists $ropmap{$op} && $ropmap{$op} eq $__op__){
                  if($param2 =~ s/^(\w+)$/$1 *oprnd1/o){
                  }elsif($param2 =~ s/const\s+?(\w+)\s*[^\w\s]*\s*\w+$/$1 *oprnd1/o){
                  }elsif($param2 =~ s/const\s+?(\w+)\s*[^\w\s]*$/$1 *oprnd1/o){
                  }elsif($param2 =~ s/(\w+)\s*[^\w\s]*\s*\w+$/$1 *oprnd1/o){ }
                  $w = "$rettype ${class}_$__op__(${param2}, ${param1}) {\n";
               }else{
                  if($param1 =~ s/^(\w+)$/$1 *oprnd1/o){
                  }elsif($param1 =~ s/const\s+?(\w+)\s*[^\w\s]*\s*\w+$/$1 *oprnd1/o){
                  }elsif($param1 =~ s/const\s+?(\w+)\s*[^\w\s]*$/$1 *oprnd1/o){
                  }elsif($param1 =~ s/(\w+)\s*[^\w\s]*\s*\w+$/$1 *oprnd1/o){ }
                  $w = "$rettype ${class}_$__op__(${param1}, ${param}) {\n";
               }
               $w .= "   return operator$op (*oprnd1, oprnd2);\n}\n";
            }
         }else{
            $h .= "     return *self $op oprnd2;\n   }\n";
         }
      }
      $addmethods{$class}.=$h;
      if ($w ne ""){$wrapper{$class}.=$w;}
      next;
   }
   print;
}
close(IFILE);

if($class ne ""){
   for $c(keys %unopfreq) {$unopcount{$class}{$c} = $unopfreq{$c};}
   for $c(keys %cmpfreq) {$cmpcount{$class}{$c} = $cmpfreq{$c};}
}

print "\n// ----------- Operator overloading code ----------\n";
for $c(keys %addmethods){
   %unopfreq=%{$unopcount{$c}};
   %cmpfreq=%{$cmpcount{$c}};
# Add the corresponding unary __pos__ if __neg__ is defined
   if ($unopfreq{'-'} > 0 && $unopfreq{'+'} == 0){
      $h =  "   $c __pos__(){\n";
      $h .= "     return *self;\n   }\n";
      $addmethods{$c}.=$h;
   }

# Add the __cmp__ function if the appropriate operators exist
   if ($cmpfreq{'=='} > 0 && $cmpfreq{'<'} > 0 && $cmpfreq{'>'} > 0){
      $h =  "   int __cmp__(const ${c}& oprnd2){\n";
      $h .= "     int retval;\n";
      $h .= "     if (*self == oprnd2)\n";
      $h .= "       retval = 0;\n";
      $h .= "     else if (*self < oprnd2)\n";
      $h .= "       retval = -1;\n";
      $h .= "     else if (*self > oprnd2)\n";
      $h .= "       retval = 1;\n";
      $h .= "     return retval;\n   }\n";
      $addmethods{$c}.=$h;
   }

   print "$addmethods{$c}\n};\n";
}

print "\n// ----------- Wrapped operator overloading code ----------\n";
for $c(keys %wrapper){
   print "$wrapper{$c}\n";
   print "#ifdef __cplusplus\nextern \"C\" {\n#endif\n\n";
   print "%}\n";
}

if ($foundendif){
   print $lastline;
}
