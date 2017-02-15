#!/bin/sh
while true; do
  inotifywait -e modify -e create -e delete -f *.py .
  ps auxw | grep manage.py
  kill
done
