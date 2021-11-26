# sparse_redes_cpp
A C++ program to sparsify mobility networks as long as they come in the same format as the files in citiesMult.
Change the `state` string variable to the name of the state, or the directory where the network files are.

Build with the command:
`c++ .\sparse_redes.cpp -o sparse_redes.exe -Ofast` on Windows, `c++ .\sparse_redes.cpp -o sparse_redes -Ofast` on Linux

There needs to be 3 text files inside the state directory, called: `mobnetwork.txt`, `Poparea.txt`, `Citypatch.txt`.

The program will output a file named `results.txt` with the matrix of Resistances of that network, and `newmobnetwork.txt` with the sparsified matrix.

The resistances matrix takes a long time to compute, so after the first time, inside `main.cpp` comment the macro `#define ESCRIBIR_MATRIZ` out, and leave the macro `#define LEER_MATRIZ` uncommented.
This way you can sparsify the matrix to different amounts of links and not calculate the resistances every time.

## Input files format
### mobnetwork.txt
A 3 column, space separated values text file where:

1. The first column contains the origin patches id
2. The second column contains the destination patches id
3. The third column contains the link weights

The network must be layed out so that the first column is strictly sorted ascending by value, then the second column after that. There can't be any duplicated links
(eg. the first two columns are *key*). Otherwise the program will fail.

### Poparea.txt
A 3 column, space separated values text file where:

1. The first column contains the patch id
2. The second column contains the population of the patch
3. The third column contains the area of the patch

### Citypatch.txt
Characteristic of multiscale networks. A 2 column, space separated values text file where:

1. The first column contains the patch id
2. The second column contains the area of the patch

## Precautions
- Make sure the `Poparea.txt` and `Citypatch.txt` contain all the patches from 0 to N-1 that appear on the network.
- If a patch does not appear on the network, but it is described in the `Poparea.txt` and `Citypatch.txt` files, it will mean there is a disconnected node, and the Laplacian will
be more degenerated, thus making it possible for the pseudo-inverse calculation to fail.
- Do not skip numbers on the numeration of the patch id's. They must range from 0 to N-1.
