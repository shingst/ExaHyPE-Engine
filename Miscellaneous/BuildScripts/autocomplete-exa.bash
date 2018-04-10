# Bash completion for exa
# First, we just list the commands.
# As a next step, all installed applications should be given.

# see also https://github.com/jarun/googler/blob/master/auto-completion/bash/googler-completion.bash
# for further examples.

_exa() {
  _exa_commands=$(exa help-shortlist)

  local cur
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  COMPREPLY=( $(compgen -W "${_exa_commands}" -- ${cur}) )

  return 0
}
complete -o nospace -F _exa  exa
