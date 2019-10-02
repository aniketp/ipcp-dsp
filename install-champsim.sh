#!/usr/bin/env bash

champdir="../champsim"
command="git clone https://github.com/ChampSim/ChampSim.git ${champdir}"
if [ ! -d ${champdir} ]; then
    eval ${command}
else
    echo "Champsim already installed."; exit 1
fi

