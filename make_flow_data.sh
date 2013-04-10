#!/bin/bash

FILE=out

print_carbon_flows ${FILE} '' "[product = 'he2' and product = 'gamma']" > idl/c2_flows1
print_carbon_flows ${FILE} '' "[product = 'he2' and product = 'n']" > idl/c2_flows2
print_carbon_flows ${FILE} '' "[reactant = 'he2' and reactant = 'gamma']" > idl/c2_flows3
print_carbon_flows ${FILE} '' "[reactant = 'he2' and reactant = 'n']" > idl/c2_flows4

print_carbon_flows ${FILE} '' "[(product = 'h2' and product = 'gamma') or (reactant = 'h2' and reactant = 'gamma') or (reactant = 'h2' and reactant = 'electron')]" > idl/co_flows

print_carbon_flows ${FILE} '' "[reactant = 'li3' and reactant = 'h1']" > idl/c3_flows1
print_carbon_flows ${FILE} '' "[reactant = 'li3' and reactant = 'gamma']" > idl/c3_flows2
print_carbon_flows ${FILE} '' "[reactant = 'li3' and reactant = 'n']" > idl/c3_flows3

print_carbon_flows ${FILE} '' "[reactant = 'be4' and reactant = 'gamma']" > idl/c4_flows1
print_carbon_flows ${FILE} '' "[reactant = 'be4' and reactant = 'n']" > idl/c4_flows2
print_carbon_flows ${FILE} '' "[reactant = 'li3' and reactant = 'h1']" > idl/c4_flows3

