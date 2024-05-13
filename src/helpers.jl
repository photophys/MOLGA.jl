module Helpers

using PeriodicTable
export element_symbol, atomic_number

"""
    element_symbol(atomic_number)

Get the symbol of an element of the periodic table as `String`.

# Arguments

  - `atomic_number::Int`

# Example

```jldoctest; setup=:(using MOLGA.Helpers)
julia> element_symbol(1)
"H"
```
"""
element_symbol(atomic_number::Int) = String(elements[atomic_number].symbol)

"""
    atomic_number(element_symbol)

Get the atomic number of an element of the periodic table as `Int`.

# Arguments

  - `element_symbol::AbstractString`: The case-sensitive element symbol.

# Example

```jldoctest; setup=:(using MOLGA.Helpers)
julia> atomic_number("Ar")
18
```
"""
atomic_number(element_symbol::AbstractString) = elements[Symbol(element_symbol)].number

end
