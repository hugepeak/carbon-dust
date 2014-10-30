#! /usr/bin/python

import os, commands

cmds = []
cmds.append( "../../nucnet-tools-code/examples/analysis/print_properties" )
cmds.append( "../../nucnet-tools-code/examples/analysis/print_mass_fractions_in_zones" )

length = 1
for cmd in cmds:
    tmp = len( commands.getoutput("ls %s" % cmd).split() )
    if tmp > length: length = tmp

if length != 1:
    os.system("./make_analysis.sh")
