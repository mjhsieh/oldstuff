# Only tested on Mac OS X + git 1.6

# VCS
_bold=$(tput bold)
_normal=$(tput sgr0)
_vcs_readlink() {
  if [ -d $1 ]; then
    (cd $1; pwd)
  else
    return 1
  fi
}
__vcs_dir() {
  local vcs base_dir sub_dir ref
  sub_dir() {
    local sub_dir
    sub_dir=$(_vcs_readlink "${PWD}")
    sub_dir=${sub_dir#$1}
    echo ${sub_dir#/}
  }
  git_dir() {
    base_dir=$(git rev-parse --show-cdup 2>/dev/null) || return 1
    base_dir=$(_vcs_readlink "$base_dir/..")
    sub_dir=$(git rev-parse --show-prefix)
    sub_dir=${sub_dir%/}
    ref=$(git symbolic-ref -q HEAD || git name-rev --name-only HEAD 2>/dev/null)
    ref=${ref#refs/heads/}
    vcs="git"
  }
  cvs_dir() {
    [ -d "CVS" ] && [ -d "CVS/Repository" ] || return 1
    cvsbase="."
    while [ -d "${cvsbase}/../CVS" ]; do
      cvsbase="${cvsbase}/.."
    done
    cvsbase="$(_vcs_readlink ${cvsbase})"
    cvsbranch="$(< CVS/Repository)"
    rrn=${cvsbase##*/}
    cvsbranch=${cvsbranch##${rrn}/}
    [ -z ${cvsbranch} ] && cvsbranch=${rrn}
    ref="${cvsbranch}:"
    vcs="cvs"
  }
  svn_dir() {
    [ -d ".svn" ] || return 1
    base_dir="."
    while [ -d "$base_dir/../.svn" ]; do base_dir="$base_dir/.."; done
    base_dir=$(_vcs_readlink "$base_dir")
    sub_dir=$(sub_dir "${base_dir}")
    ref=$(svn info "$base_dir" | awk '/^URL/ { sub(".*/","",$0); r=$0 } /^Revision/ { sub("[^0-9]*","",$0); print r":"$0 }')
    vcs="svn"
  }
#  svk_dir() {
#    [ -f ~/.svk/config ] || return 1
#    base_dir=$(awk '/: *$/ { sub(/^ */,"",$0); sub(/: *$/,"",$0); if (match("'${PWD}'", $0"(/|$)")) { print $0; d=1; } } /depotpath/ && d == 1 { sub(".*/","",$0); r=$0 } /revision/ && d == 1 { print r ":" $2; exit 1 }' ~/.svk/config) && return 1
#    ref=${base_dir##*
#}
#    base_dir=${base_dir%%
#*}
#    sub_dir=$(sub_dir "${base_dir}")
#    vcs="svk"
#  }
  hg_dir() {
    base_dir="."
    while [ ! -d "$base_dir/.hg" ]; do base_dir="$base_dir/.."; [ $(_vcs_readlink "${base_dir}") = "/" ] && return 1; done
    base_dir=$(_vcs_readlink "$base_dir")
    sub_dir=$(sub_dir "${base_dir}")
    ref=$(< "${base_dir}/.hg/branch")
    vcs="hg"
  }
  git_dir ||
  svn_dir ||
#  svk_dir ||
  hg_dir ||
  cvs_dir ||
  base_dir="$PWD"
# echo "${vcs:+($vcs)}${_bold}${base_dir/$HOME/~}${_normal}${vcs:+\[$ref\]${_bold}${sub_dir}${_normal}}"
  if [ -n "$vcs" ]; then
    __vcs_prefix="($vcs)"
    __vcs_base_dir="${base_dir/$HOME/~}"
    __vcs_ref="[$ref]"
    __vcs_sub_dir="${sub_dir}"
  else
    __vcs_prefix=''
    __vcs_base_dir="${PWD/$HOME/~}"
    __vcs_ref=''
    __vcs_sub_dir=''
  fi
}
PROMPT_COMMAND=__vcs_dir
PS1='\u@\h:$__vcs_prefix\[${_bold}\]${__vcs_base_dir}\[${_normal}\]${__vcs_ref}\[${_bold}\]${__vcs_sub_dir}\[${_normal}\]\$ '
