#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Chemical elements and compunds"""

__author__ = "shiva senthilkumar"

from collections import namedtuple

# TODO: define the Element datatype (using namedtuple)

Element = namedtuple("Element", ["number", "symbol", "name", "weight"])
# Element = ...

elements = []


def load_elements(filename):
    """Load text database of chemical elements to initialize the module.

    Add the Element objects to the global elements list"""
    with open(filename) as db:
        for line in db:
            line_parts = line.split()  # splitting line into 4 parts
            number = int(line_parts[0])
            symbol = str(line_parts[1])
            name = str(line_parts[2])
            weight = float(line_parts[3])
            correct_type = Element(number, symbol, name, weight)
            elements.append(correct_type)


def element_by_number(number):
    """Return the chemical element of a given atomic number

    Return None if element is not found.
    """
    for each_element in elements:
        if each_element[0] == number:
            return elements[each_element[0] - 1]
    else:
        return None


def element_by_symbol(symbol):
    """Return the chemical element of a given symbol

    Return None if element is not found.
    """
    for each_element in elements:
        if each_element[1] == symbol:
            return elements[each_element[0] - 1]
    else:
        return None


def element_by_name(name):
    """Return the chemical element of a given name.

    Name matching is case-insensitive.
    Return None if element is not found.
    """
    name = name.upper()
    for each_element in elements:
        if each_element[2].upper() == name:
            return elements[each_element[0] - 1]
    else:
        return None


def compound_elements(formula):
    """Return the list of elements in a given compound formula."""
    
    indexes = []
    for character in formula:
        if character.isupper():
            indexes.append(formula.index(character))
    
    next_indexes = [x + 1 for x in indexes]
    
    for each_index in next_indexes:
        if formula[each_index].isdigit():
            return (element_by_symbol(formula[each_index-1])
        else:
            return (element_by_symbol(formula[each_index-1:each_index+1]))
            
        
def compound_weight(formula):
    """Return the total weight of a given compound formula."""
    pass


if __name__ == "__main__":
    load_elements("elements.txt")

    # Feel free to use (umcomment) these examples to try out your code

    print("The 17. element:", element_by_number(17).name)
    print("The name of Ge:", element_by_symbol("Ge").name)
    print("The weight of Oxygen:", element_by_name("oxygen").weight)

    print("The elements in salt:", compound_elements("H2O"))
    # print("The molecular weight of Caffeine", compound_weight("C8H10N4O2"))
