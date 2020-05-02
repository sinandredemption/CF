# CF
A C++ library for continued fraction arithmetic based on Gosper's algorithm

**Prerequisite**

CF requires MPIR library, including it's C++ wrapper

**Quickstart guide**

Include the file cf.h and cf.cpp in your compilation.
A class named "CF", which stands for continued fraction, is provided.
You can initialize a CF from a fraction (<code>CF pi(355, 113)</code>) or progressively add terms using the <code>add_term()</code> function.
See the function pi_simple() from demo.cpp for a quick idea of how to use the library.

**Features**

- Addition, subtraction, multiplication, and division of continued fractions
- Arithmetic between continued fractions and fractions
- To keep this library simple-to-use, arithmetic using operator overloading can be done
- Precision-capped arithmetic is supported
- Included demo.cpp, which uses the library to calculate continued fraction of pi

## Author and License

Syed Fahad, distributed under MIT License
