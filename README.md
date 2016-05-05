
New algorithms are based on formula: (1+f')·X·tan(E) - b·e'^2·sin(E) - Y = 0 . 

This formula is transformed afterwards to quartic equaton: t^4 + A·t^3 + B·t - 1 = 0 .

Four algorithms are created for solving this quartic equation: 1) direct solution (please, see my 'Quartic' repository) 2) optimized direct solution (taylored for this specific case) 3) Newton-Raphson based iterative solution 4) Solution based on so-called semiquadratic interpolation. The main idea behind all these methods is - speed. To simplify procedure and to reduce using of trancedential functions, in order to enhance the performancies.

New algorithms are tested against the most popular and well-known methods: Moritz-Heiskanen (implemented in proj), Bowring, Borkowsky (recommended by IERS). The tested methods are enumerated and explained in 'methods.txt'. The comparison results are presented in 'results.txt'. In a word - new methods are more than  promissing. They can cope with the old ones and they win. There is a clear logical explaination why they are so fast(er) and I believe that they will be used everywhere where performance matters. (especially - nr2) 

And, performances are everything! They save time and energy; they save the planet's resources.
