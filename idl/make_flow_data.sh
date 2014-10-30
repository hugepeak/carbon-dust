#!/bin/bash

FILE="../out"

print_carbon_flows ${FILE} '' "[product = 'he2' and product = 'gamma']" > c2_flows1
print_carbon_flows ${FILE} '' "[product = 'he2' and product = 'n']" > c2_flows2
print_carbon_flows ${FILE} '' "[reactant = 'he2' and reactant = 'gamma']" > c2_flows3
print_carbon_flows ${FILE} '' "[reactant = 'he2' and reactant = 'n']" > c2_flows4

print_carbon_flows ${FILE} '' "[(product = 'h2' and product = 'gamma') or (reactant = 'h2' and reactant = 'gamma') or (reactant = 'h2' and reactant = 'electron')]" > co_flows

print_carbon_flows ${FILE} '' "[reactant = 'li3' and reactant = 'h1']" > c3_flows1
print_carbon_flows ${FILE} '' "[reactant = 'li3' and reactant = 'gamma']" > c3_flows2
print_carbon_flows ${FILE} '' "[reactant = 'li3' and reactant = 'n']" > c3_flows3

print_carbon_flows ${FILE} '' "[reactant = 'be4' and reactant = 'gamma']" > c4_flows1
print_carbon_flows ${FILE} '' "[reactant = 'be4' and reactant = 'n']" > c4_flows2
print_carbon_flows ${FILE} '' "[reactant = 'li3' and reactant = 'h1']" > c4_flows3

