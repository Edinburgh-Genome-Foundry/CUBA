#!/bin/sh
./rqworker_start.sh & watchmedo shell-command \
    --patterns="*.py;*.txt" \
    --recursive \
    --command='kill $(cat rqworker_pid) ; ./rqworker_start.sh'
    .
