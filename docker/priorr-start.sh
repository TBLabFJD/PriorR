#!/bin/sh
# Starts PriorR as the user priorr (not root). The data volume may come from an older image that
# ran as root: it is given to priorr first (only when it isn't already its owner).
set -e
if [ "$(id -u)" = "0" ]; then
  mkdir -p /var/priorr
  if [ "$(stat -c %U /var/priorr)" != "priorr" ]; then chown -R priorr:priorr /var/priorr; fi
  exec setpriv --reuid=priorr --regid=priorr --init-groups "$@"
fi
exec "$@"
