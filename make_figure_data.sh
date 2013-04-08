#!/bin/bash

FILE=out

print_properties ${FILE} time t9 rho > idl/prop
print_mass_fractions_in_zones ${FILE} n h1 h2 he2 li3 be4 b5 c6C n7C o8C o8R > idl/mass_fractions

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

compute_carbon_flows ${FILE} '' > flows
sed_flow_names.sh
