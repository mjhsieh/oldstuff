#!/usr/bin/env perl
$numArgs = $#ARGV + 1;
$ifile=shift @ARGV;
if($numArgs > 1){
   $fwd_mod=shift @ARGV;
   $mod_defined=1;
}else{
   $fwd_mod="";
   $mod_defined=0;
}

use File::Basename;
($base,$path,$type) = fileparse($ifile);
$dot = rindex($base, '.');
if ($dot > 0){
   $basename = substr($base, 0, $dot);
}else{
   $basename = $base;
}

open(F,"swig -python -shadow -c++ $ifile 2>&1 |");
$errCount=0;
while(<F>){
   if(/multiply defined /o){
#      /Line (\d+)\s*\.\s*(\w+)\s+(\S+)/o;
      /Line (\d+)\s*\.\s*((static\s+)?\w+)\s+(\S+)/o;
      $n=$1;
      $t=$2;
      $f=$4;
      $t=~s/static\s*//o;
      $md{$n}=[$t,$f,0];
      if(/\(member /o){
         $md{$n}[2]=1;
      }
   }
}
close(F);

$firstclass="";
$lastline="";
$foundendif=0;
$inclass=0;

open(IFILE,$ifile);
unshift @line,"";
while(<IFILE>){
  if(m|^\s*//|o){print; next;}
  if (/#endif/o){
      $foundendif = 1;
      $lastline = $_;
      next;
   }elsif($foundendif){
      print $lastline;
      $foundendif = 0;
   }
   #$line[$num]="%name(${name}__L$num) ".$line[$num];
   if(/%module\s+(\S+)/o){
      $mod=$1;
      if (!$mod_defined){
         $fwd_mod=$mod;
         $mod_defined=1;
      }
   }
   if(/class\s+(\w+)/o){
      $class=$1;
      if($firstclass eq ""){
        $firstclass=$class;
        open(F,">${basename}_defs.defs");
        }
      $inclass=1;
   }elsif(/%addmethods\s+(\w+)/o){
      $class=$1;
      $inclass=0;
   }
   if(exists $md{$.}){
      $s=$_;
      if($_ !~ /\)\s*;/o && $_ =~ /\)\s*\{/o){
         while(<IFILE>){
            $s.=$_;
            if(/\)\s*;/o || /\)\s*\{/o){ last;}
         }
      }
      $c=0;
      $s =~ s/,/$c++,","/oge;
      if($s !~ /\(\s*(void)?\s*\)/m){
         $c++;
      }
      $numargs = $c;

# Is this a Constructor?
      $isconst = ($md{$.}[0] =~ /(Const)/o);
# Is this a member function?
      $ismember = ($md{$.}[2]);

# Determine the return type, the number and types of parameters
      $func = $md{$.}[1];
      $s=~ /^\s*(\w+\s*\W?)[\s\w]*\(/o;
      $rettype=$1;
      $rettype=~ s/&//o;
      $rettype=~ s/^\s*//o;
      $rettype=~ s/\s*$//o;
      $s=~ /\(\s*(.*,)*\s*(.*)\s*\)/o;
      $param1=$1;
      $param2=$2;
#      print "rettype = $rettype numargs = $c \n";
      @params=();
      while ($param1 =~ /\s*(.*?),/o){
         push (@params, $1);
         $param1=$';
      }
      push (@params, $param2);
      $numparams = scalar(@params);

      if($ismember || $isconst){

# Member functions actually have one more argument, self!
         if($ismember){
            $numargs++;
         }
# Is this a copy constructor?
         if($isconst && $numparams == 1){
            $param2 =~ /((const\s+)?\w+)/o;
            $type=$1;
            $type =~ s/const\s*//o;
            if($type eq $class){
# If so, define a deepcopy method using %name
# It will be called from a special __deepcopy__ in the derived class below
               $numargs = -$c;
               print "%name(${class}__deepcopy) $s";
               if($firstclass ne ""){
                  print F "def ${class}__deepcopy(*args):\n"; }
            }
         }

         push @{$redef{$class}{$md{$.}[1]}},[$.,$numargs,$rettype];

      }else{
         push @{$redef{""}{$md{$.}[1]}},[$.,$numargs,$rettype];
      }
      $md{$.}[1] =~ s/^new_//o;
      print "%name($md{$.}[1]__L$.) $s";

# Keep track of all the renamed external def's that we might want to
# clean up later
      if((($inclass && !$ismember) || $isconst) && $firstclass ne ""){
         print F "def $md{$.}[1]__L$.(*args):\n"; }

      next;
   }
   print;
}
close(IFILE);
if($firstclass ne ""){
   close(F);
}

#
# Clean up swig generated files
#
$wrapfile = "${basename}_wrap.c";
$pyfile = "${mod}.py";
#system("rm ${wrapfile} ${pyfile}");
@filelist = ($wrapfile, $pyfile);
$cnt = unlink @filelist;

print "%pragma(python) include=\"$basename\_overloaded.py\"\n";

open(F,">${basename}_overloaded.py");
print F "import types, re\n";
for $cl(keys %redef){
 $deepcopy=0;
 if($cl ne ""){
    if(!exists $defmethods{$cl}){
       $defmethods{$cl}="\n";
    }
#   print F "class $cl\($cl):\n";
   print F "class __Dummy:\n";
   $spaces = "   ";
 }else{
   $spaces = "";
 }
 $spacer1 = $spaces;
   for $f(keys %{$redef{$cl}}){
      $el=""; 
      if($f =~ /^new_/o){
         $clname="";
         $defmethods{$cl}.="${cl}.__init__ = __Dummy.__dict__['__init__']\n";
         print F "${spacer1}def __init__\(self,*args):\n";
      }else{
         $clname="${cl}_";
         $defmethods{$cl}.="${cl}.${f} = __Dummy.__dict__['${f}']\n";
         print F "${spacer1}def $f\(*args):\n";
      }
      print F "${spacer1}   trydefault = 0\n";

      @funcids = @{$redef{$cl}{$f}};
      for $i(@funcids){
         if ($i->[1] < 0){
            $i->[1] = -$i->[1];
            $deepcopy = 1;
         }
      }
      @sortedfuncids = sort {($a->[1] <=> $b->[1])} @funcids;

#      print "Sorted Funcs\n";
#      for $i(@sortedfuncids){
#         print "$i->[0], $i->[1], $i->[2]\n";
#      }

      $numfuncs = scalar(@sortedfuncids);
      $numargs = $sortedfuncids[0][1];
      if ($numfuncs > 1){
         $has_multiple = ($numargs == $sortedfuncids[1][1]);
      }else{
         $has_multiple = 0;
      }
      $islast = 0;
      $fcount = 0;
      $funcname = $f;
      if($funcname =~ s/^new_// || $cl eq ""){
         $defname=$funcname;
      }else{
         $defname="${cl}.${funcname}";
      }

      $ifunc = 1;
      for $i(@sortedfuncids){
         $rettype = $i->[2];
         $base_rettype = $rettype;
         $base_rettype =~ s/\*//o;
#         print "rettype = $rettype base_rettype = $base_rettype \n";
         if ($numargs != $i->[1]){
            $numargs = $i->[1];
            if ($has_multiple){
               $spacer1 = $spaces;
            }
            $fcount = 0;
            if ($ifunc < $numfuncs){
               $has_multiple = ($numargs == $sortedfuncids[$ifunc][1]);
            }else{
               $has_multiple = 0;
            }
            $islast = 0;
         }elsif($fcount > 0){
            $spacer1 = "$spacer1   ";
            if($ifunc < $numfuncs){
               $islast = ($numargs != $sortedfuncids[$ifunc][1]);
            }else{
               $islast = 1;
            }
         }

         if ($fcount == 0){
            if ($el eq ""){
               print F "${spacer1}   try:\n"
            }
            print F "$spacer1      ${el}if len(args)==$i->[1]:\n";
         }
         $fcount++;

         if ($has_multiple){
            if(!$islast){
               print F "$spacer1         try:\n";
               $spacer2 = "   ";
            }else{
               $spacer2 = "";
            }
         }else{
            $spacer2= "";
         }
         if($f =~ /^new_/o){
            print F "${spacer1}${spacer2}         self.this = apply(${fwd_mod}c.${f}__L$i->[0],args)\n";
         }else{
            if($f =~ /^__/o){
               $defmethods{$cl}.="del ${cl}.__dict__['_${cl}${f}__L$i->[0]']\n";
            }else{
               $defmethods{$cl}.="del ${cl}.__dict__['${f}__L$i->[0]']\n";
            }
            print F "${spacer1}${spacer2}         val = apply(${fwd_mod}c.${clname}${f}__L$i->[0],args)\n";
            if($rettype eq $cl){
               print F "${spacer1}${spacer2}         if val: val = ${cl}Ptr(val) ; val.thisown = 1\n";
            }else{
               print F "${spacer1}${spacer2}         if val and isinstance(val, types.StringType) and re.search('_p_${base_rettype}', val):\n";
	       print F "${spacer1}${spacer2}            val = ${base_rettype}Ptr(val) ; val.thisown = 1\n";
            }
         }
         if ($has_multiple && !$islast){
            print F "$spacer1         except:\n";
         }
         $el="el";
         $ifunc++;
      }
      if ($has_multiple){
         $spacer1 = $spaces;
      }
      print F "$spacer1      else:\n";
      print F "$spacer1         trydefault = 1\n";
      if($f =~ /^new_/o){
         print F "$spacer1         self.this = apply(${fwd_mod}c.${f},args)\n";
      }else{
         print F "$spacer1         val = apply(${fwd_mod}c.${clname}${f},args)\n";
         print F "$spacer1         if val and isinstance(val, types.StringType):\n";
         print F "$spacer1            mo = re.search(r'^_[\\da-fA-F]*_p_(.*)', val)\n";
         print F "$spacer1            if mo:\n";
         print F "$spacer1               rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'\n";
         print F "$spacer1               val = eval(rettype) ; val.thisown = 1\n";
      }
      print F "$spacer1   except StandardError, e:\n";
      print F "$spacer1      if not trydefault:\n";
      print F "$spacer1         laste = e\n";
      print F "$spacer1         try:\n";
      if($f =~ /^new_/o){
         print F "$spacer1            self.this = apply(${fwd_mod}c.${f},args)\n";
      }else{
         print F "$spacer1            val = apply(${fwd_mod}c.${clname}${f},args)\n";
         print F "$spacer1            if val and isinstance(val, types.StringType):\n";
         print F "$spacer1               mo = re.search(r'^_[\\da-fA-F]*_p_(.*)', val)\n";
         print F "$spacer1               if mo:\n";
         print F "$spacer1                  rettype = val[mo.start(1):mo.end(1)] + 'Ptr(val)'\n";
         print F "$spacer1                  val = eval(rettype) ; val.thisown = 1\n";
      }
      print F "$spacer1         except StandardError, e:\n";
      print F "$spacer1            error = 'No ${defname} methods apply using the given arguments.\\n'\n";
      print F "$spacer1            if re.search('arguments', e.args[0]):\n";
      print F "$spacer1               error = error + 'Possible error is ' + laste.args[0]\n";
      print F "$spacer1            else:\n";
      print F "$spacer1               error = 'One possible error is ' + e.args[0] + '\\n'\n";
      print F "$spacer1               error = error + 'Another possible error is ' + laste.args[0]\n";
      if($f =~ /^new_/o){
         print F "$spacer1            self.thisown = 0\n";
      }
      print F "$spacer1            raise Error, error\n";
      print F "$spacer1      else:\n";
      if($f =~ /^new_/o){
         print F "$spacer1         self.thisown = 0\n";
      }
      print F "$spacer1         if re.search('arguments', e.args[0]):\n";
      print F "$spacer1            if len(args)==1:\n";
      print F "$spacer1               nargs = '1 argument'\n";
      print F "$spacer1            else:\n";
      print F "$spacer1               nargs = str(len(args)) + ' arguments'\n";
      print F "$spacer1            error = 'No ${defname} methods take ' + nargs\n";
      print F "$spacer1            raise Error, error\n";
      print F "$spacer1         else:\n";
      print F "$spacer1            raise e\n";
      if($f =~ /^new_/o){
         print F "$spacer1   self.thisown = 1\n";
         print F "$spacer1   return\n";
      }else{
         print F "$spacer1   return val\n";
      }
   }
# This is the __deepcopy__ that gets called from copy.deepcopy().
   if($deepcopy){
      $defmethods{$cl}.="${cl}.__deepcopy__ = __Dummy.__dict__['__deepcopy__']\n";
      print F "   def __deepcopy__(self, memo = None):\n";
      print F "      val = ${fwd_mod}c.new_${cl}__deepcopy(self)\n";
      print F "      if val: val = ${cl}Ptr(val) ; val.thisown = 1\n";
      print F "      return val\n";
   }
# Now assign these override methods to the actual class,
# and delete the mangled name methods from the class
   print F "$defmethods{$cl}\n";
 }

if ($foundendif){
   print $lastline;
}
close(F);
