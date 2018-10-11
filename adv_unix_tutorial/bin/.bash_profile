# mary lauren benton | 2015
# .bash_profile

# get any aliases and functions from bashrc
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# load shell dotfiles, and then some:
# * ~/.path can be used to extend `$PATH`.
for file in ~/.{path,bash_prompt,aliases,functions}; do
  [ -r "$file" ] && source "$file"
done
unset file

# set colorscheme
eval `dircolors -b /home/bentonml/.dircolors`

# user specific environment and startup programs
PATH=$PATH:$HOME/bin
export PATH

